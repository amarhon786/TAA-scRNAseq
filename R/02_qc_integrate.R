suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(patchwork)
})

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)

RAW_DIR <- "data/raw/GSE155468"
PROC_DIR <- "data/processed"
FIG_DIR <- "results/figures"
META_FILE <- "data/raw/sample_metadata.csv"
dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR,  recursive = TRUE, showWarnings = FALSE)

meta <- read.csv(META_FILE, stringsAsFactors = FALSE)
stopifnot(!any(is.na(meta$sex)), !any(is.na(meta$disease)))

sample_files <- list.files(RAW_DIR, pattern = "^GSM\\d+_.*\\.txt\\.gz$", full.names = TRUE)
names(sample_files) <- str_extract(basename(sample_files), "^GSM\\d+")
stopifnot(length(sample_files) == nrow(meta))

# header row is 1 col short (no gene col) so read it separately
read_one_matrix <- function(path) {
  con <- gzfile(path, "rt")
  cell_ids <- strsplit(readLines(con, n = 1), "\t")[[1]]
  close(con)
  dt <- data.table::fread(path, sep = "\t", header = FALSE, skip = 1,
                          data.table = FALSE)
  mat <- as.matrix(dt[, -1, drop = FALSE])
  rownames(mat) <- dt[[1]]
  colnames(mat) <- cell_ids
  Matrix::Matrix(mat, sparse = TRUE)
}

objs <- imap(sample_files, function(path, gsm) {
  message("Reading ", gsm, " (", basename(path), ")")
  counts <- read_one_matrix(path)
  obj <- CreateSeuratObject(counts = counts, project = gsm, min.cells = 3, min.features = 200)
  m <- meta[meta$geo_accession == gsm, ]
  obj$sample  <- gsm
  obj$sex     <- m$sex[1]
  obj$disease <- m$disease[1]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})

pdf(file.path(FIG_DIR, "qc_violin_prefilter.pdf"), width = 12, height = 4)
for (g in names(objs)) {
  print(VlnPlot(objs[[g]],
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0) + plot_annotation(title = g))
}
dev.off()

# standard cutoffs, close to what the paper used
objs <- map(objs, ~ subset(.x, subset = nFeature_RNA > 500 &
                                        nFeature_RNA < 6000 &
                                        percent.mt   < 15))

merged <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, nfeatures = 3000)
merged <- ScaleData(merged, vars.to.regress = c("percent.mt", "nCount_RNA"))
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

merged <- IntegrateLayers(object = merged, method = HarmonyIntegration,
                          orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
merged <- JoinLayers(merged)

merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)

markers <- list(
  SMC_contractile = c("MYH11", "ACTA2", "CNN1", "MYL9", "TAGLN", "LMOD1"),
  SMC_modulated   = c("KLF4", "TNFRSF11B", "VCAN", "FN1", "MGP"),
  Fibroblast      = c("DCN", "LUM", "COL1A1", "COL1A2", "PDGFRA", "PI16"),
  Endothelial     = c("PECAM1", "VWF", "CDH5", "CLDN5", "KDR"),
  Macrophage      = c("CD68", "CD163", "MRC1", "LYZ", "C1QA", "AIF1"),
  T_cell          = c("CD3D", "CD3E", "CD8A", "CD4", "TRAC"),
  B_cell          = c("CD79A", "MS4A1", "CD19"),
  NK              = c("NKG7", "KLRD1", "GNLY", "NCAM1"),
  Mast            = c("KIT", "TPSAB1", "CPA3"),
  Plasma          = c("MZB1", "JCHAIN", "IGKC")
)

for (ct in names(markers)) {
  merged <- AddModuleScore(merged, features = list(markers[[ct]]),
                           name = paste0(ct, "_score"))
}
score_cols <- paste0(names(markers), "_score1")

cluster_scores <- merged@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), .groups = "drop")

cluster_label <- cluster_scores %>%
  pivot_longer(-seurat_clusters, names_to = "celltype", values_to = "score") %>%
  group_by(seurat_clusters) %>%
  slice_max(score, n = 1) %>%
  mutate(celltype = str_remove(celltype, "_score1$")) %>%
  select(seurat_clusters, celltype)

merged$celltype <- cluster_label$celltype[match(merged$seurat_clusters,
                                                cluster_label$seurat_clusters)]

saveRDS(merged, file.path(PROC_DIR, "seurat_annotated.rds"))
write.csv(cluster_label, file.path(PROC_DIR, "cluster_celltype_map.csv"), row.names = FALSE)

pdf(file.path(FIG_DIR, "umap_overview.pdf"), width = 14, height = 5)
print(
  DimPlot(merged, group.by = "celltype", label = TRUE, repel = TRUE) +
  DimPlot(merged, group.by = "disease") +
  DimPlot(merged, group.by = "sex")
)
dev.off()

message("\nDone. Saved ", file.path(PROC_DIR, "seurat_annotated.rds"))
message("Inspect ", file.path(PROC_DIR, "cluster_celltype_map.csv"),
        " - refine cluster labels if a cluster looks miscalled before DE.")
message("Next: source('R/03_celltype_abundance.R')")
