#!/usr/bin/env Rscript

## grab_geo.R
## Run this from the 01_proc_data/ directory

library(GEOquery)

acc <- "GSE171894"

# directories relative to 01_proc_data/
base_dir    <- "../data"
supp_dir    <- file.path(base_dir, acc)                 # ../data/GSE171894
raw_out_dir <- file.path(base_dir, paste0(acc, "_raw")) # ../data/GSE171894_raw

dir.create(base_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(raw_out_dir, recursive = TRUE, showWarnings = FALSE)

message("Downloading supplementary files for ", acc, " into: ", supp_dir)

getGEOSuppFiles(
  acc,
  baseDir       = base_dir,   # puts files under ../data/GSE171894
  makeDirectory = TRUE
)

raw_tar <- file.path(supp_dir, paste0(acc, "_RAW.tar"))
if (!file.exists(raw_tar)) {
  stop("Could not find RAW tar file at: ", raw_tar)
}

message("Found RAW tar: ", raw_tar)
message("Untarring RAW data into: ", raw_out_dir)

untar(raw_tar, exdir = raw_out_dir)

message("Untar complete. RAW files are now in: ", raw_out_dir)

## Series matrix: use only for phenotype / sample metadata
message("Downloading series matrix and extracting sample metadata...")

gse_list <- getGEO(acc, GSEMatrix = TRUE)
gse      <- gse_list[[1]]

pheno <- pData(gse)
saveRDS(pheno, file.path(base_dir, paste0(acc, "_pheno.rds")))

message("Saved sample metadata: ", file.path(base_dir, paste0(acc, "_pheno.rds")))
message("Done.")

