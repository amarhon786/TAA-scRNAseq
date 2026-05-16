# pseudobulk + DESeq2 per cell type

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(tidyverse)
  library(Matrix)
})

obj <- readRDS("data/processed/seurat_annotated.rds")

MIN_CELLS_PER_GROUP <- 10
MIN_SAMPLES_PER_CT <- 4

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
md <- obj@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(sample_ct = paste(sample, celltype, sep = "__"))

cells_per <- md %>% dplyr::count(sample, celltype, sex, disease, name = "n_cells")
write.csv(cells_per, "results/tables/cells_per_sample_celltype.csv", row.names = FALSE)

keep_groups <- cells_per %>% dplyr::filter(n_cells >= MIN_CELLS_PER_GROUP)

celltypes <- keep_groups %>%
  dplyr::count(celltype, name = "n_samples") %>%
  dplyr::filter(n_samples >= MIN_SAMPLES_PER_CT) %>%
  dplyr::pull(celltype)

message("Cell types with enough samples for DE: ", paste(celltypes, collapse = ", "))

run_de_one_celltype <- function(ct) {
  message("\n=== ", ct, " ===")
  cells_in <- md %>%
    dplyr::filter(celltype == ct,
           sample %in% keep_groups$sample[keep_groups$celltype == ct]) %>%
    dplyr::pull(cell)

  sample_vec <- md$sample[match(cells_in, md$cell)]
  pb <- sapply(unique(sample_vec), function(s) {
    Matrix::rowSums(counts[, cells_in[sample_vec == s], drop = FALSE])
  })
  pb <- as.matrix(pb)
  storage.mode(pb) <- "integer"

  cd <- md %>%
    dplyr::filter(sample %in% colnames(pb)) %>%
    dplyr::distinct(sample, sex, disease) %>%
    tibble::column_to_rownames("sample")
  cd <- cd[colnames(pb), , drop = FALSE]
  # F + Control are the ref levels
  cd$sex <- factor(cd$sex, levels = c("F", "M"))
  cd$disease <- factor(cd$disease, levels = c("Control", "TAA"))

  # interaction needs all 4 sex x disease groups, otherwise fall back to additive
  use_interaction <- length(unique(cd$sex)) == 2 &&
                     length(unique(cd$disease)) == 2 &&
                     all(table(cd$sex, cd$disease) >= 1)
  design <- if (use_interaction) ~ sex + disease + sex:disease
            else                  ~ sex + disease

  dds <- DESeqDataSetFromMatrix(countData = pb, colData = cd, design = design)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, quiet = TRUE)

  write_contrast <- function(res, label) {
    out <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene") %>%
      dplyr::arrange(padj) %>%
      dplyr::mutate(celltype = ct, contrast = label, .before = 1)
    write.csv(out, sprintf("results/tables/de_%s_%s.csv", gsub("[^A-Za-z0-9]+", "_", ct), label),
              row.names = FALSE)
    out
  }

  results_list <- list()

  rn <- resultsNames(dds)

  # avg over sex = disease + 0.5*interaction
  if (use_interaction && all(c("disease_TAA_vs_Control", "sexM.diseaseTAA") %in% rn)) {
    cv <- rep(0, length(rn)); names(cv) <- rn
    cv["disease_TAA_vs_Control"] <- 1
    cv["sexM.diseaseTAA"] <- 0.5
    results_list$TAA_vs_Control <- write_contrast(results(dds, contrast = cv), "TAA_vs_Control")
  } else if ("disease_TAA_vs_Control" %in% rn) {
    results_list$TAA_vs_Control <- write_contrast(
      results(dds, name = "disease_TAA_vs_Control"),
      "TAA_vs_Control"
    )
  }

  if (use_interaction && all(c("disease_TAA_vs_Control", "sexM.diseaseTAA") %in% rn)) {
    results_list$Male_TAA_vs_Male_Control <- write_contrast(
      results(dds, contrast = list(c("disease_TAA_vs_Control", "sexM.diseaseTAA"))),
      "Male_TAA_vs_Male_Control"
    )
  }

  if (use_interaction && all(c("sex_M_vs_F", "sexM.diseaseTAA") %in% rn)) {
    results_list$Male_TAA_vs_Female_TAA <- write_contrast(
      results(dds, contrast = list(c("sex_M_vs_F", "sexM.diseaseTAA"))),
      "Male_TAA_vs_Female_TAA"
    )
  }

  if ("disease_TAA_vs_Control" %in% rn) {
    results_list$Female_TAA_vs_Female_Control <- write_contrast(
      results(dds, name = "disease_TAA_vs_Control"),
      "Female_TAA_vs_Female_Control"
    )
  }

  if ("sex_M_vs_F" %in% rn) {
    results_list$Male_Control_vs_Female_Control <- write_contrast(
      results(dds, name = "sex_M_vs_F"),
      "Male_Control_vs_Female_Control"
    )
  }

  results_list
}

de_all <- lapply(celltypes, function(ct) {
  tryCatch(run_de_one_celltype(ct),
           error = function(e) { message("skipped ", ct, ": ", e$message); NULL })
})
names(de_all) <- celltypes

saveRDS(de_all, "data/processed/de_results.rds")
message("\nAll DE done. Tables in results/tables/de_*.csv")
message("Next: source('R/05_enrichment.R')")
