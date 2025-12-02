#devtools::install_github("guokai8/rcellmarker")
library(Seurat)
library(rcellmarker)
library(dplyr)
library(ggplot2)

seu <- readRDS("../data/GSE171894_seurat_harmony.rds")

seu <- JoinLayers(seu)

# Make sure identities are the Harmony-based clusters
Idents(seu) <- "seurat_clusters"


markers <- FindAllMarkers(
  seu,
  only.pos        = TRUE,   # only positive markers
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Keep only the columns rcellmarker needs
markers_for_cm <- markers[, c("cluster", "gene", "avg_log2FC", "p_val_adj")]
head(markers_for_cm)



# Tissue markers ----------------------------------------------------------
dat <- rcellmarker:::.getdata(species = "human", db = "default")
tissues_human <- sort(unique(dat$tissueType))


# pick relevant tissues (only those that actually exist in tissues_human)
relevant_tissues <- c(
  "Uterine cervix"
)

relevant_tissues <- intersect(relevant_tissues, tissues_human)
relevant_tissues


res_cm_tissue <- cellMarker(
  x      = markers_for_cm,
  type    = "seurat",
  species = "human",
  keytype = "SYMBOL",
  tissue  = relevant_tissues,  # <— key change
  weight  = 0.25,              # modest weight on logFC
  padj    = 0.05,              # only significant markers
  minSize = 5                  # require ≥5 genes per cell-type call
)


cm_df_tissue <- marker(res_cm_tissue)   #


cluster2ct_tissue <- cm_df_tissue %>%
  count(Cluster, cellType, name = "n_genes") %>%
  group_by(Cluster) %>%
  slice_max(order_by = n_genes, n = 1, with_ties = FALSE) %>%
  ungroup()

cluster2ct_tissue # microglial cells in cervical tissue... suspicious

#  Build a mapping from cluster -> cell type
cluster2ct_tissue$Cluster <- as.character(cluster2ct_tissue$Cluster)

cluster_map <- setNames(cluster2ct_tissue$cellType,
                        cluster2ct_tissue$Cluster)   # names = cluster IDs

#  Use that mapping to get a per-cell vector, then DROP names
celltype_vec <- cluster_map[as.character(seu$seurat_clusters)]
celltype_vec <- unname(celltype_vec)  # or as.vector(celltype_vec)

#  Assign to Seurat metadata
seu$celltype_cm <- celltype_vec

# sanity check
table(seu$seurat_clusters, seu$celltype_cm)


# weird -------------------------------------------------------------------

weird_clusters <- c(7, 9, 11, 15)

for (cl in weird_clusters) {
  cat("\n==== Cluster", cl, "====\n")
  print(
    FindMarkers(seu, ident.1 = cl, only.pos = TRUE, logfc.threshold = 0.25) %>%
      head(20)
  )
}


# manual annotations ------------------------------------------------------

# # start from the CellMarker label
# seu$celltype_cm_raw <- seu$celltype_cm
# 
# # initialize a new, coarse annotation
# seu$celltype_coarse <- NA_character_
# 
# # T / NK clusters
# seu$celltype_coarse[seu$seurat_clusters %in% c(0, 2, 4, 6, 10, 13, 16)] <- "T/NK"
# 
# # monocytes / myeloid
# seu$celltype_coarse[seu$seurat_clusters %in% c(1, 3, 8, 7)] <- "Myeloid-like"
# 
# # epithelial / tumor
# seu$celltype_coarse[seu$seurat_clusters %in% c(5, 12, 9, 11, 15)] <- "Epithelial/tumor-like"
# 
# # plasmacytoid DC
# seu$celltype_coarse[seu$seurat_clusters %in% c(14)] <- "pDC"
# 
# # sanity check
# table(seu$seurat_clusters, seu$celltype_coarse)


# final -------------------------------------------------------------------

# define a manual mapping from cluster ID -> final cell type
final_map <- c(
  "0"  = "T/NK",
  "1"  = "Monocyte/Myeloid",
  "2"  = "T/NK",
  "3"  = "Monocyte/Myeloid",
  "4"  = "T/NK",                    # could refine to Treg-like if markers support it
  "5"  = "Epithelial",
  "6"  = "T/NK",
  "7"  = "B cell",                 # was 'Microglial cell'
  "8"  = "Monocyte/Myeloid",
  "9"  = "Basal epithelial/tumor", # was 'Migration phase fetal germ cell'
  "10" = "T/NK",
  "11" = "Cycling tumor/epithelial", # was 'Neural progenitor cell'
  "12" = "Epithelial",
  "13" = "T/NK",
  "14" = "pDC",
  "15" = "CAF / stromal",          # was 'Leydig precursor cell'
  "16" = "T/NK"
)

# turn cluster IDs into a character vector
clusters_chr <- as.character(seu$seurat_clusters)

# sanity check: all clusters have a mapping?
setdiff(unique(clusters_chr), names(final_map))
# should print character(0); if not, you forgot a cluster in final_map

# build per-cell label vector, drop names so Seurat doesn't get confused
celltype_vec <- unname(final_map[clusters_chr])

seu@meta.data$celltype_final <- celltype_vec

p_cm <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "celltype_final",
  label    = TRUE,
  repel = TRUE
) + ggtitle("GSE171894 Harmony UMAP\ncell types") +
  theme(legend.position="bottom")

res_dir <- "../res"

# Save as PNG
ggsave(
  filename = file.path(res_dir, "GSE171894_umap_harmony_cellmarker.png"),
  plot     = p_cm,
  width    = 10,
  height   = 6,
  dpi      = 300
)

saveRDS(seu, file.path("..", "data", "GSE171894_seurat_harmony_annotated.rds"))


table(seu$seurat_clusters, seu$celltype_final)
