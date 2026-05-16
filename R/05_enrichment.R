# ORA + GSEA per contrast

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(fgsea)
  library(msigdbr)
  library(org.Hs.eg.db)
  library(patchwork)
})
source("R/theme_nature.R")

de_all <- readRDS("data/processed/de_results.rds")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables",  recursive = TRUE, showWarnings = FALSE)

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:KEGG_LEGACY")
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcollection = "GO:BP")

to_named_list <- function(df) split(df$gene_symbol, df$gs_name)
sets <- list(
  Hallmark = to_named_list(hallmark),
  KEGG     = to_named_list(kegg),
  GO_BP    = to_named_list(gobp)
)

PADJ_CUT <- 0.1
LFC_CUT <- 0.5

run_enrichment <- function(de_df) {
  sig_up   <- de_df %>% dplyr::filter(padj < PADJ_CUT, log2FoldChange >  LFC_CUT) %>% dplyr::pull(gene)
  sig_down <- de_df %>% dplyr::filter(padj < PADJ_CUT, log2FoldChange < -LFC_CUT) %>% dplyr::pull(gene)
  universe <- de_df$gene

  ora_objs <- list()
  ora_df   <- tibble()
  for (sn in names(sets)) {
    term2gene <- stack(sets[[sn]])[, c("ind", "values")]
    names(term2gene) <- c("term", "gene")
    for (dir in c("up", "down")) {
      gene_set <- if (dir == "up") sig_up else sig_down
      if (length(gene_set) < 5) next
      # cutoffs = 1 so no terms get dropped on p/q
      e <- enricher(gene_set, TERM2GENE = term2gene, universe = universe,
                    pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
      if (is.null(e) || nrow(e@result) == 0) next
      key <- paste(sn, dir, sep = "_")
      ora_objs[[key]] <- e
      ora_df <- bind_rows(ora_df, as.data.frame(e) %>% mutate(geneset = sn, direction = dir))
    }
  }

  rank_vec <- de_df %>%
    dplyr::filter(!is.na(stat)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    { setNames(.$stat, .$gene) } %>%
    sort(decreasing = TRUE)

  gsea_df <- map_dfr(names(sets), function(sn) {
    res <- fgsea(pathways = sets[[sn]], stats = rank_vec, minSize = 10, maxSize = 500)
    as.data.frame(res) %>%
      mutate(geneset = sn,
             leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
  })

  list(ora_objs = ora_objs, ora_df = ora_df,
       gsea_df = gsea_df, rank_vec = rank_vec)
}

enrich_all <- list()
for (ct in names(de_all)) {
  if (is.null(de_all[[ct]])) next
  for (contrast in names(de_all[[ct]])) {
    label <- paste(ct, contrast, sep = "__")
    message("Enriching: ", label)
    enrich_all[[label]] <- run_enrichment(de_all[[ct]][[contrast]])
    safe <- gsub("[^A-Za-z0-9]+", "_", label)
    write.csv(enrich_all[[label]]$ora_df, sprintf("results/tables/enrichment_ORA_%s.csv", safe),
              row.names = FALSE)
    write.csv(enrich_all[[label]]$gsea_df, sprintf("results/tables/enrichment_GSEA_%s.csv", safe),
              row.names = FALSE)
  }
}
saveRDS(enrich_all, "data/processed/enrichment_results.rds")

dotplot_kegg <- function(gsea_df, label, top_n = 20) {
  df <- gsea_df %>%
    dplyr::filter(geneset == "KEGG", padj < 0.25) %>%
    dplyr::slice_min(padj, n = top_n) %>%
    dplyr::mutate(
      pathway_short = stringr::str_remove(pathway, "^(HALLMARK|KEGG_LEGACY|KEGG|GOBP)_"),
      pathway_short = stringr::str_replace_all(pathway_short, "_", " "),
      pathway_short = stringr::str_to_sentence(pathway_short),
      pathway_short = stringr::str_wrap(pathway_short, width = 42),
      direction     = ifelse(NES > 0, "Up", "Down")
    )
  if (nrow(df) == 0) return(NULL)
  ggplot(df, aes(NES, reorder(pathway_short, NES))) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "black") +
    geom_segment(aes(xend = 0, yend = pathway_short),
                 color = "grey85", linewidth = 0.3) +
    geom_point(aes(size = size, fill = -log10(padj)),
               shape = 21, color = "black", stroke = 0.3) +
    scale_fill_viridis_c(option = "magma", end = 0.9, name = "-log10\n(padj)") +
    scale_size_continuous(range = c(2.2, 6), name = "Set\nsize",
                          breaks = scales::breaks_pretty(4)) +
    labs(x = "Normalized enrichment score (NES)",
         y = NULL,
         title = paste0("KEGG enrichment - ", label),
         subtitle = paste0("Top ", nrow(df), " pathways (padj < 0.25)")) +
    theme(axis.text.y = element_text(size = 7.5, lineheight = 0.85),
          panel.grid.major.y = element_line(color = "grey95", linewidth = 0.2),
          legend.box = "vertical",
          legend.spacing.y = unit(1, "mm"))
}

# sizes are in mm, /25.4 for inches
cairo_pdf("results/figures/enrichment_dotplots.pdf",
          width  = NATURE_DOUBLE_COL / 25.4,
          height = 150 / 25.4,
          onefile = TRUE)
for (label in names(enrich_all)) {
  p <- dotplot_kegg(enrich_all[[label]]$gsea_df, gsub("__", " - ", label))
  if (!is.null(p)) print(p)
}
dev.off()

plot_gsea_running <- function(en, label) {
  top <- en$gsea_df %>%
    dplyr::filter(padj < 0.25) %>%
    dplyr::slice_min(padj, n = 1)
  if (nrow(top) == 0) return(NULL)
  pathway_genes <- sets[[top$geneset]][[top$pathway]]
  if (is.null(pathway_genes)) return(NULL)
  plotEnrichment(pathway_genes, en$rank_vec) +
    labs(title = stringr::str_replace_all(
                   stringr::str_remove(top$pathway, "^(HALLMARK|KEGG|GOBP)_"),
                   "_", " "),
         subtitle = sprintf("%s   NES=%.2f   padj=%.2g",
                            label, top$NES, top$padj),
         x = "Gene rank", y = "Enrichment score")
}

pdf("results/figures/enrichment_gsea_running.pdf",
    width = NATURE_ONE_HALF / 25.4, height = 70 / 25.4)
for (label in names(enrich_all)) {
  p <- plot_gsea_running(enrich_all[[label]],
                         gsub("__", " - ", label))
  if (!is.null(p)) print(p)
}
dev.off()

message("Enrichment done.")
message("  Tables : results/tables/enrichment_*.csv")
message("  Figures: results/figures/enrichment_dotplots.pdf")
message("           results/figures/enrichment_gsea_running.pdf")
message("Next: source('R/06_figures.R')")
