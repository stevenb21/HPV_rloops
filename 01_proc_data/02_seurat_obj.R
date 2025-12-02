#!/usr/bin/env Rscript

library(Seurat)
library(Matrix)

# directory that contains GSM5236544_HPV+.1.txt.gz etc.
raw_dir <- "../data/GSE171894_raw"

# list all .txt.gz files
files <- list.files(raw_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)

if (length(files) == 0) {
  stop("No .txt.gz files found in: ", raw_dir)
}

# make nice sample IDs from filenames
sample_ids <- sub("\\.txt\\.gz$", "", basename(files))
names(files) <- sample_ids

make_seurat_from_file <- function(sample_id, f) {
  message("Reading ", f)

  # genes in rows, cells in columns
  mat <- read.table(
    f,
    header      = TRUE,
    row.names   = 1,
    sep         = "\t",
    check.names = FALSE
  )

  # convert to matrix / sparse matrix
  mat <- as.matrix(mat)
  mat <- Matrix(mat, sparse = TRUE)

  # build Seurat object
  so <- CreateSeuratObject(
    counts       = mat,            # normalized UMI counts, but usable
    project      = "GSE171894",
    min.cells    = 3,
    min.features = 200
  )

  # annotate with sample ID
  so$sample_id <- sample_id

  # prefix cell names with sample_id to keep them unique
  colnames(so) <- paste(sample_id, colnames(so), sep = "_")

  return(so)
}

# build one Seurat object per file
objs_list <- mapply(
  make_seurat_from_file,
  sample_id = sample_ids,
  f         = files,
  SIMPLIFY  = FALSE
)

# name the list by sample id
names(objs_list) <- sample_ids

# optional: save the list
saveRDS(objs_list, "../data/GSE171894_seurat_list.rds")

# optional: merge all samples into one Seurat object
if (length(objs_list) > 1) {
  merged <- Reduce(function(x, y) merge(x, y),
                   objs_list)
} else {
  merged <- objs_list[[1]]
}

saveRDS(merged, "../data/GSE171894_seurat_merged.rds")

message("Built Seurat objects for samples: ",
        paste(sample_ids, collapse = ", "))
message("Saved:")
message("  - per-sample objects: ../data/GSE171894_seurat_list.rds")
message("  - merged object:      ../data/GSE171894_seurat_merged.rds")

