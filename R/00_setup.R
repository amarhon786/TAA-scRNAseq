if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::init(bare = TRUE, restart = FALSE)

cran_pkgs <- c("BiocManager", "Seurat", "tidyverse", "Matrix", "patchwork", "ggrepel", "pheatmap",
               "RColorBrewer", "data.table")

bioc_pkgs <- c(
  "GEOquery",
  "DESeq2",
  "clusterProfiler",
  "fgsea",
  "msigdbr",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "SingleCellExperiment",
  "scran",
  "scater",
  "DropletUtils"
)

to_install <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")

# no bioc version pin, BiocManager picks the one that matches this R
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc)) BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)

renv::snapshot(prompt = FALSE)

invisible(lapply(c(cran_pkgs, bioc_pkgs), function(p) {
  ok <- requireNamespace(p, quietly = TRUE)
  message(sprintf("  %s %s", if (ok) "[ok]" else "[MISSING]", p))
}))

message("\nSetup complete. Next: source('R/01_download.R')")
