#!/usr/bin/env Rscript

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(ragg)
set.seed(123)

#-----------------------------------
# 0. Input / output paths
#-----------------------------------
data_dir <- "../data"
res_dir  <- "../res"

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

postqc_path <- file.path(data_dir, "GSE171894_seurat_postqc.rds")
harmony_rds <- file.path(data_dir, "GSE171894_seurat_harmony.rds")

#-----------------------------------
# 1. Load post-QC Seurat object
#-----------------------------------
seu <- readRDS(postqc_path)

message("Loaded post-QC object with ", ncol(seu), " cells and ", nrow(seu), " genes.")

#-----------------------------------
# 2. Set batch variable for Harmony
#-----------------------------------
# Change this if your batch variable has a different name
batch_var <- "sample_id"

if (!batch_var %in% colnames(seu@meta.data)) {
  stop(
    "Batch variable '", batch_var,
    "' not found in seu@meta.data. Available columns: ",
    paste(colnames(seu@meta.data), collapse = ", ")
  )
}

message("Using '", batch_var, "' as batch variable for Harmony.")

#-----------------------------------
# 3. Ensure PCA exists
#-----------------------------------
# Your previous script already ran PCA, but this makes the script robust
if (!"pca" %in% names(seu@reductions)) {
  message("No PCA reduction found; running PCA on variable features.")
  if (length(VariableFeatures(seu)) == 0) {
    stop("No variable features found. Make sure FindVariableFeatures was run before Harmony.")
  }
  seu <- RunPCA(seu, features = VariableFeatures(seu))
} else {
  message("PCA reduction found; reusing existing PCA.")
}

# Number of PCs to use for Harmony / neighbors / UMAP
n_pcs <- 30

#-----------------------------------
# 4. Run Harmony integration
#-----------------------------------
message("Running Harmony integration on ", n_pcs, " PCs...")

seu <- RunHarmony(
  object        = seu,
  group.by.vars = batch_var,
  dims.use      = 1:n_pcs
)

message("Harmony integration complete. Harmony reduction stored in seu@reductions$harmony.")

#-----------------------------------
# 5. Neighbors, clustering, UMAP on Harmony
#-----------------------------------
message("Running neighbors, clustering, and UMAP on Harmony embedding...")

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_pcs)
seu <- FindClusters(seu, resolution = 0.6)

seu <- RunUMAP(
  seu,
  reduction = "harmony",
  dims = 1:n_pcs
)

#-----------------------------------
# 6. UMAP plots (Harmony)
#-----------------------------------
umap_harmony_clusters_pdf <- file.path(res_dir, "GSE171894_umap_harmony_clusters.pdf")
pdf(umap_harmony_clusters_pdf, width = 6, height = 5)
p_umap_harmony_clusters <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("GSE171894 Harmony UMAP - clusters")
print(p_umap_harmony_clusters)
dev.off()

umap_harmony_batch_pdf <- file.path(res_dir, "GSE171894_umap_harmony_batch.pdf")
pdf(umap_harmony_batch_pdf, width = 6, height = 5)
p_umap_harmony_batch <- DimPlot(
  seu,
  reduction = "umap",
  group.by = batch_var
) +
  ggtitle(paste0("GSE171894 Harmony UMAP - ", batch_var)) +
  theme(legend.position="bottom")
print(p_umap_harmony_batch)
dev.off()


umap_harmony_batch_png <- file.path(res_dir, "GSE171894_umap_harmony_batch.png")

ragg::agg_png(
  filename = umap_harmony_batch_png,
  width = 6, height = 5, units = "in", res = 300
)

p_umap_harmony_batch <- DimPlot(
  seu,
  reduction = "umap",
  group.by = batch_var
) +
  ggtitle(paste0("GSE171894 Harmony UMAP - ", batch_var)) +
  theme(
    legend.position = "bottom"
  ) + guides(color = guide_legend(nrow=2,byrow=TRUE))

print(p_umap_harmony_batch)
dev.off() 


#-----------------------------------
# 7. Save Harmony-integrated object
#-----------------------------------
saveRDS(seu, harmony_rds)

message("Harmony pipeline complete.")
message("Outputs:")
message("  - Harmony UMAP (clusters): ", umap_harmony_clusters_pdf)
message("  - Harmony UMAP (", batch_var, "): ", umap_harmony_batch_pdf)
message("  - Harmony-integrated Seurat RDS: ", harmony_rds)


#-----------------------------------
# 8. Cell counts per cluster per sample
#-----------------------------------

# make sure required columns exist
if (!all(c("seurat_clusters", "sample_id") %in% colnames(seu@meta.data))) {
  stop("Expected 'seurat_clusters' and 'sample_id' in meta.data but did not find them.")
}

# contingency table: clusters x samples
cluster_sample_tab <- table(
  seurat_cluster = seu$seurat_clusters,
  sample_id      = seu$sample_id
)

# convert to long-format data.frame (one row per cluster-sample)
cluster_sample_df <- as.data.frame(cluster_sample_tab)
colnames(cluster_sample_df) <- c("seurat_cluster", "sample_id", "n_cells")

# write to file
cellcount_path <- file.path(res_dir, "GSE171894_cell_counts_per_cluster_per_sample.tsv")
write.table(
  cluster_sample_df,
  file      = cellcount_path,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

message("  - Cell count table: ", cellcount_path)

