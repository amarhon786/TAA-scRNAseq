suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})
source("R/theme_nature.R")

obj <- readRDS("data/processed/seurat_annotated.rds")
de_all <- readRDS("data/processed/de_results.rds")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

umap_pt <- function(group, palette = NULL, title = group, label = FALSE, override_size = 3) {
  p <- DimPlot(obj, group.by = group, label = label, repel = TRUE, label.size = 2.6,
               label.box = FALSE, pt.size = 0.6, raster = TRUE, raster.dpi = c(1200, 1200)) +
       ggtitle(title) +
       guides(color = guide_legend(override.aes = list(size = override_size, alpha = 1),
                                   ncol = 1)) +
       theme(legend.position   = "bottom",
             legend.direction  = "vertical",
             legend.key.size   = unit(2.5, "mm"),
             legend.title      = element_blank(),
             axis.line         = element_line(linewidth = 0.3),
             plot.title        = element_text(hjust = 0.5))
  if (!is.null(palette)) p <- p + scale_color_manual(values = palette)
  p
}
p_umap <- umap_pt("celltype", label = TRUE, title = "Cell type") +
  umap_pt("disease", DISEASE_COLORS, title = "Disease") +
  umap_pt("sex", SEX_COLORS, title = "Sex") +
  plot_layout(nrow = 1)
save_nature(p_umap, "results/figures/umap_panel.pdf",
            width_mm = NATURE_DOUBLE_COL, height_mm = 110)

# pb over all cells per sample, only used for the sample PCA
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
sample_vec <- obj$sample
pb_all <- sapply(unique(sample_vec), function(s) {
  Matrix::rowSums(counts[, sample_vec == s, drop = FALSE])
})
storage.mode(pb_all) <- "integer"
cd_all <- obj@meta.data %>%
  tibble::as_tibble() %>%
  dplyr::distinct(sample, sex, disease) %>%
  tibble::column_to_rownames("sample")
cd_all <- cd_all[colnames(pb_all), , drop = FALSE]

dds_all <- DESeqDataSetFromMatrix(pb_all, cd_all, design = ~ disease)
vsd <- vst(dds_all, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("disease", "sex"), returnData = TRUE)

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = disease, shape = sex)) +
  geom_point(size = 2.4, stroke = 0.4) +
  geom_text_repel(aes(label = name), size = 2.2, color = "grey20",
                  segment.size = 0.2, min.segment.length = 0, box.padding = 0.3) +
  scale_color_manual(values = DISEASE_COLORS) +
  scale_shape_manual(values = c(F = 16, M = 17)) +
  labs(title = "Pseudobulk PCA - disease and sex",
       x = sprintf("PC1 (%.0f%%)", attr(pca_data, "percentVar")[1] * 100),
       y = sprintf("PC2 (%.0f%%)", attr(pca_data, "percentVar")[2] * 100))

save_nature(p_pca, "results/figures/pseudobulk_pca.pdf",
            width_mm = NATURE_SINGLE_COL, height_mm = 85)

volcano <- function(df, title) {
  df <- df %>% dplyr::filter(!is.na(padj))
  df$sig <- df$padj < 0.05 & abs(df$log2FoldChange) > 0.5
  df$dir <- ifelse(df$log2FoldChange > 0, "Up", "Down")
  df$col <- with(df, ifelse(!sig, "ns", ifelse(dir == "Up", "Up", "Down")))
  top_lab <- df %>% dplyr::filter(sig) %>% dplyr::arrange(padj) %>% head(15)

  ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = col), size = 0.6, alpha = 0.7, stroke = 0) +
    geom_text_repel(data = top_lab, aes(label = gene),
                    size = 2.1, color = "black",
                    segment.size = 0.2, min.segment.length = 0,
                    box.padding = 0.25, max.overlaps = 20) +
    scale_color_manual(values = c(ns = "grey80", Up = "#D55E00", Down = "#0072B2"),
                       guide = "none") +
    labs(title = title,
         x = expression(log[2]~fold~change),
         y = expression(-log[10]~adj.~italic(P)))
}

pdf("results/figures/volcano_panel.pdf", width = NATURE_ONE_HALF / 25.4, height = 95 / 25.4)
for (ct in names(de_all)) {
  if (is.null(de_all[[ct]])) next
  for (contrast in names(de_all[[ct]])) {
    print(volcano(de_all[[ct]][[contrast]], paste0(ct, " - ", contrast)))
  }
}
dev.off()

ma_plot <- function(df, title) {
  df <- df %>% dplyr::filter(!is.na(padj), baseMean > 0)
  df$sig <- df$padj < 0.05 & abs(df$log2FoldChange) > 0.5
  df$col <- with(df, ifelse(!sig, "ns", ifelse(log2FoldChange > 0, "Up", "Down")))
  top_lab <- df %>% dplyr::filter(sig) %>%
    dplyr::arrange(padj) %>% head(12)

  ggplot(df, aes(baseMean, log2FoldChange)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "black") +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(aes(color = col), size = 0.6, alpha = 0.7, stroke = 0) +
    geom_text_repel(data = top_lab, aes(label = gene),
                    size = 2.1, color = "black",
                    segment.size = 0.2, min.segment.length = 0,
                    box.padding = 0.25, max.overlaps = 20) +
    scale_x_log10(labels = scales::label_log()) +
    scale_color_manual(values = c(ns = "grey80", Up = "#D55E00", Down = "#0072B2"),
                       guide = "none") +
    labs(title = title,
         x = "Mean of normalized counts (log scale)",
         y = expression(log[2]~fold~change))
}

pdf("results/figures/ma_panel.pdf",
    width  = NATURE_ONE_HALF / 25.4,
    height = 95 / 25.4)
for (ct in names(de_all)) {
  if (is.null(de_all[[ct]])) next
  for (contrast in names(de_all[[ct]])) {
    print(ma_plot(de_all[[ct]][[contrast]],
                  paste0(ct, " - ", contrast)))
  }
}
dev.off()

smc_key <- intersect(c("SMC_contractile", "SMC_modulated"), names(de_all))
if (length(smc_key) > 0) {
  smc_key <- smc_key[1]
  sex_de <- de_all[[smc_key]]$M_vs_F
  if (!is.null(sex_de)) {
    top_genes <- sex_de %>%
      filter(!is.na(padj)) %>% arrange(padj) %>% head(30) %>% pull(gene)

    smc_obj <- subset(obj, celltype == smc_key)
    avg <- AverageExpression(smc_obj, features = top_genes, group.by = c("disease", "sex"),
                             slot = "data")$RNA

    ann_col <- data.frame(
      disease = sub("_.*",  "", colnames(avg)),
      sex     = sub(".*_",  "", colnames(avg)),
      row.names = colnames(avg)
    )
    pheatmap(
      avg,
      scale = "row",
      cluster_cols = FALSE,
      annotation_col = ann_col,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      main = paste0("Top sex-DE genes in ", smc_key),
      filename = "results/figures/heatmap_smc_sex_de.pdf",
      width = 6, height = 8
    )
  }
}

sanity <- intersect(c("RPS4Y1", "DDX3Y", "EIF1AY", "XIST", "CARTPT"), rownames(obj))
# drop all-zero genes, they just give empty panels
sanity <- sanity[vapply(sanity, function(g) {
  max(Seurat::GetAssayData(obj, assay = "RNA", layer = "data")[g, ]) > 0
}, logical(1))]
message("Sanity features kept (non-zero): ", paste(sanity, collapse = ", "))

if (length(sanity) > 0) {
  n_features <- length(sanity)

  feat_h_mm <- 45 * n_features
  cairo_pdf("results/figures/sanity_featureplot.pdf",
            width  = NATURE_DOUBLE_COL / 25.4,
            height = feat_h_mm / 25.4,
            onefile = TRUE)
  print(FeaturePlot(obj, features = sanity, split.by = "sex",
                    raster = TRUE, raster.dpi = c(1200, 1200),
                    pt.size = 0.7, order = TRUE) &
        theme(plot.title  = element_text(size = 8, face = "bold"),
              axis.title  = element_text(size = 6),
              axis.text   = element_text(size = 5),
              axis.line   = element_line(linewidth = 0.25),
              axis.ticks  = element_line(linewidth = 0.25),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 6),
              legend.key.size = unit(2.5, "mm"),
              strip.text  = element_text(size = 7)))
  dev.off()

  cairo_pdf("results/figures/sanity_violin.pdf",
            width  = NATURE_DOUBLE_COL / 25.4,
            height = (22 * n_features + 25) / 25.4,
            onefile = TRUE)
  print(VlnPlot(obj, features = sanity, group.by = "celltype",
                split.by = "sex", pt.size = 0,
                stack = TRUE, flip = TRUE, cols = SEX_COLORS) &
        theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 6),
              axis.text.y  = element_text(size = 5),
              axis.title   = element_text(size = 7),
              strip.text.y = element_text(size = 7, face = "bold", angle = 0),
              strip.background = element_blank(),
              legend.position = "bottom",
              legend.text = element_text(size = 6),
              legend.key.size = unit(2.5, "mm"),
              plot.title = element_blank()))
  dev.off()

  old_file <- "results/figures/sanity_y_chrom_cartpt.pdf"
  if (file.exists(old_file)) file.remove(old_file)
}

message("Figures done. See results/figures/")
