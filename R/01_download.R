# GSE155468 (Li 2020): 8 TAA (4M/4F) + 3 ctrl (1M/2F)

suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

GSE_ID <- "GSE155468"
RAW_DIR <- "data/raw"
META_FILE <- "data/raw/sample_metadata.csv"

dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)

message("Fetching ", GSE_ID, " metadata from GEO...")
gse <- getGEO(GSE_ID, GSEMatrix = TRUE, getGPL = FALSE)
pheno <- pData(gse[[1]])

write.csv(pheno, file.path(RAW_DIR, "GSE155468_phenoData_full.csv"), row.names = FALSE)

# col names are specific to this series, check the pheno csv if reusing
sex_col <- "gender:ch1"
disease_col <- "disease state:ch1"
stopifnot(all(c(sex_col, disease_col) %in% names(pheno)))

meta <- tibble(
  geo_accession = pheno$geo_accession,
  title         = as.character(pheno$title),
  sex = case_when(
    str_to_lower(pheno[[sex_col]]) == "female" ~ "F",
    str_to_lower(pheno[[sex_col]]) == "male"   ~ "M",
    TRUE ~ NA_character_
  ),
  disease = case_when(
    str_detect(str_to_lower(pheno[[disease_col]]), "control|normal|non-aneurysmal") ~ "Control",
    str_detect(str_to_lower(pheno[[disease_col]]), "ataa|aneurysm")                 ~ "TAA",
    TRUE ~ NA_character_
  )
)

write.csv(meta, META_FILE, row.names = FALSE)
message("Sample metadata written to ", META_FILE)
print(meta)

message("Downloading supplementary files (this can take a few minutes)...")
supp <- getGEOSuppFiles(GSE_ID, baseDir = RAW_DIR, makeDirectory = TRUE, fetch_files = TRUE)
print(supp)

tar_files <- list.files(file.path(RAW_DIR, GSE_ID),
                        pattern = "\\.tar$", full.names = TRUE)
for (tf in tar_files) {
  message("Extracting ", basename(tf))
  untar(tf, exdir = file.path(RAW_DIR, GSE_ID))
}

message("\nDownload complete. Next: source('R/02_qc_integrate.R')")
message("Inspect ", META_FILE, " - if any sex/disease is NA, edit by hand before continuing.")

