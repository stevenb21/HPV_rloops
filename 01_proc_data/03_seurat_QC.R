#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)
library(patchwork)

#-----------------------------------
# 0. Output dirs
#-----------------------------------
data_dir    <- "../data"
res_dir <- "../res"

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

#-----------------------------------
# 1. Load merged object
#-----------------------------------
seu_path <- file.path(data_dir, "GSE171894_seurat_merged.rds")
seu <- readRDS(seu_path)

message("Loaded object with ", ncol(seu), " cells and ", nrow(seu), " genes.")

#-----------------------------------
# 2. Basic QC metrics
#-----------------------------------
# % mitochondrial genes (human: MT-)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# If you want ribosomal too, uncomment:
# seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RPL|^RPS")
message("0")
# Violin plots (multi-feature => patchwork object; use print())
qc_violin_pdf <- file.path(res_dir, "GSE171894_qc_violin.pdf")
pdf(qc_violin_pdf, width = 10, height = 4)
p_vln <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
print(p_vln)
dev.off()
message("1")
# Scatter QC plot (simple ggplot; + theme_bw() is fine here)
qc_scatter_pdf <- file.path(res_dir, "GSE171894_qc_scatter.pdf")
pdf(qc_scatter_pdf, width = 8, height = 4)
p_scatter <- FeatureScatter(
  seu,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) 
print(p_scatter)
dev.off()
message("2")
#-----------------------------------
# 3. Filter low-quality cells
#-----------------------------------
# tweak thresholds after inspecting QC plots
min_features <- 200
max_features <- 6000
max_mt       <- 20   # percent

seu <- subset(
  seu,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt   < max_mt
)

message("After QC filtering: ", ncol(seu), " cells remain.")

#-----------------------------------
# 4. Normalization, variable features, scaling
#-----------------------------------
# Counts are already “normalized UMIs” in GEO, but we still log-normalize
seu <- NormalizeData(
  seu,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

seu <- FindVariableFeatures(
  seu,
  selection.method = "vst",
  nfeatures = 2000
)

top10 <- head(VariableFeatures(seu), 10)
message("Top 10 variable features: ", paste(top10, collapse = ", "))

seu <- ScaleData(seu, features = VariableFeatures(seu))

#-----------------------------------
# 5. PCA
#-----------------------------------
seu <- RunPCA(seu, features = VariableFeatures(seu))

elbow_pdf <- file.path(res_dir, "GSE171894_elbowplot.pdf")
pdf(elbow_pdf, width = 6, height = 4)
p_elbow <- ElbowPlot(seu, ndims = 50)
print(p_elbow)
dev.off()

# choose number of PCs (30 is a reasonable default)
n_pcs <- 30

#-----------------------------------
# 6. Neighbors, clustering, UMAP
#-----------------------------------
seu <- FindNeighbors(seu, dims = 1:n_pcs)
seu <- FindClusters(seu, resolution = 0.6)

seu <- RunUMAP(seu, dims = 1:n_pcs)

# UMAP by clusters
umap_clusters_pdf <- file.path(res_dir, "GSE171894_umap_clusters.pdf")
pdf(umap_clusters_pdf, width = 6, height = 5)
p_umap_clusters <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) + ggtitle("GSE171894 UMAP - clusters")
print(p_umap_clusters)
dev.off()

# UMAP by sample_id if available
if ("sample_id" %in% colnames(seu@meta.data)) {
  umap_sample_pdf <- file.path(res_dir, "GSE171894_umap_sample.pdf")
  pdf(umap_sample_pdf, width = 6, height = 5)
  p_umap_sample <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "sample_id"
  ) + ggtitle("GSE171894 UMAP - sample_id")
  print(p_umap_sample)
  dev.off()
}


#-----------------------------------
# 8. Save processed object
#-----------------------------------
postqc_rds <- file.path(data_dir, "GSE171894_seurat_postqc.rds")
saveRDS(seu, postqc_rds)

message("Pipeline complete.")
message("Outputs:")
message("  - QC violin:   ", qc_violin_pdf)
message("  - QC scatter:  ", qc_scatter_pdf)
message("  - Elbow plot:  ", elbow_pdf)
message("  - UMAP plots:  ", umap_clusters_pdf,
        if ("sample_id" %in% colnames(seu@meta.data)) paste0(", ", umap_sample_pdf) else "")
message("  - Seurat RDS:  ", postqc_rds)

