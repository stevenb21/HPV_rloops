library(data.table)
library(patchwork)

seu <- readRDS(file.path("..", "data", "GSE171894_seurat_harmony_annotated.rds"))

seu$HPV_status <- ifelse(grepl("HPV\\+", seu$sample_id), "HPV+", "HPV-")
seu$HPV_status <- factor(seu$HPV_status, levels = c("HPV+", "HPV-"))



# Check:
#table(seu$sample_id, seu$HPV_status)


# Rloop regulators --------------------------------------------------------


# grabbed "Gene list of R-loop regulators" under 2. R-loop regulators:
# https://rloopbase.nju.edu.cn/download.jsp
rloop_path <- "../data/gene_list_of_R-loop_regulators.gz"

rloop_genes <- fread(rloop_path)
rloop_genes <- rloop_genes[rloop_genes$Species == "Homo sapiens", 1:8]


rlrg_symbols <- unique(rloop_genes$Regulator)

# Intersect with genes present in your Seurat object
rlrg_symbols <- intersect(rlrg_symbols, rownames(seu))
length(rlrg_symbols)
head(rlrg_symbols)

DefaultAssay(seu) <- "RNA"  # or whatever your raw/normalized assay is called



# Rloop score per-cell ----------------------------------------------------

seu <- AddModuleScore(
  object   = seu,
  features = list(rlrg_symbols),
  name     = "RloopScore",
  assay    = DefaultAssay(seu)
)

# wip below ---------------------------------------------------------------

# ------------------------------------------------------------------
# Choose malignant cells only
# ------------------------------------------------------------------
# Adjust this to match your malignant label(s):
malignant_labels <- c("Epithelial",
                      "Cycling tumor/epithelial",
                      "Basal epithelial/tumor")  # e.g. your celltype_final label for tumor cells

seu_malignant <- subset(
  seu,
  subset = celltype_final %in% malignant_labels
)

# Extract per-cell R-loop scores + metadata
df_mal <- FetchData(
  seu_malignant,
  vars = c("RloopScore1", "sample_id", "HPV_status", "celltype_final")
) %>%
  rename(
    celltype = celltype_final,
    RloopScore = RloopScore1
  )

# quick sanity check
table(df_mal$HPV_status)
table(df_mal$sample_id, df_mal$HPV_status)


# ------------------------------------------------------------------
# Tumor-level R-loop score for malignant cells
# ------------------------------------------------------------------

sample_scores <- df_mal %>%
  group_by(sample_id, HPV_status) %>%
  summarise(
    mean_Rloop = mean(RloopScore, na.rm = TRUE),
    median_Rloop = median(RloopScore, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

sample_scores




# average the normalized expression of your RLRG genes
#
# subtract the average of a control gene set matched by expression level
#
# So you get a relative score (higher = relatively enriched for RLRG expression compared to background, in that cell).



# Viz ---------------------------------------------------------------------



p_rl <- FeaturePlot(
  seu,
  features  = "RloopScore1",
  reduction = "umap",
  order     = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95"
)

res_dir <- "../res"

# Save as PNG
ggsave(
  filename = file.path(res_dir, "GSE171894_umap_rloop_cells.png"),
  plot     = p_rl,
  width    = 10,
  height   = 6,
  dpi      = 300
)


# HPV status --------------------------------------------------------------

# Derive HPV status from sample_id
# sample_id examples:
# "GSM5236544_HPV+.1" "GSM5236545_HPV+.2" "GSM5236546_HPV-.1" "GSM5236547_HPV-.2"



# Set identities to cell type (replace "celltype_col" with your column)
Idents(seu) <- "celltype_final"

df_rl <- FetchData(
  seu,
  vars = c("RloopScore1", "sample_id", "HPV_status", "celltype_final")
) %>%
  rename(
    celltype = celltype_final
  ) %>%
  filter(!is.na(celltype))

p_violin_box <- ggplot(df_rl, aes(x = celltype, y = RloopScore1, fill = HPV_status)) +
  geom_violin(position = position_dodge(width = 0.8), trim = FALSE) +
  geom_boxplot(width = 0.15, position = position_dodge(width = 0.8), outlier.size = 0.3) +
  labs(
    x = "Cell type",
    y = "R-loop module score",
    fill = "HPV status",
    title = "R-loop scores across cell types by HPV status"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(
  filename = file.path(res_dir, "GSE171894_violin_box_RloopScore_by_celltype_HPV_status.png"),
  plot     = p_violin_box,
  width    = 14,
  height   = 6,
  dpi      = 300
)


# If not already created:
df_rl <- FetchData(
  seu,
  vars = c("RloopScore1", "sample_id", "HPV_status", "celltype_final")
) %>%
  rename(celltype = celltype_final)

# Global violin + boxplot (all cells)
p_global_violin_box <- ggplot(df_rl, aes(x = HPV_status, y = RloopScore1, fill = HPV_status)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.size = 0.4) +
  labs(
    x = "HPV status",
    y = "R-loop module score",
    title = "Global R-loop scores by HPV status"
  ) +
  theme_bw()



# immune_plots ------------------------------------------------------------

# Pull metadata
meta <- seu@meta.data

# If your sample column is called something else, e.g. orig.ident, do:
# meta$sample <- meta$orig.ident

# (Optional) keep only major immune cell types
immune_types <- c("B cell", "T/NK", "Monocyte/Myeloid", "pDC")
meta_immune <- meta %>%
  filter(celltype_final %in% immune_types)

frac_by_status <- meta_immune %>%
  group_by(HPV_status, celltype_final) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

immune_frac <- ggplot(frac_by_status, aes(x = HPV_status, y = frac, fill = celltype_final)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "HPV status",
    y = "Fraction of cells",
    fill = "Cell type",
    title = "Fraction of major immune cell types by HPV status"
  ) +
  theme_bw()

ggsave("../res/immune_cell_frac.png", plot = immune_frac)


# exhaustion --------------------------------------------------------------
# Exhaustion score per-cell ---------------------------------------------

exhaustion_genes <- c(
  "PDCD1",   # PD-1
  "CTLA4",
  "LAG3",
  "HAVCR2",  # TIM-3
  "TIGIT",
  "TOX",
  "EOMES",
  "CXCL13",
  "ENTPD1",  # CD39
  "TNFRSF9"  # 4-1BB
)

# Keep only genes present in the object
exhaustion_genes <- intersect(exhaustion_genes, rownames(seu))
length(exhaustion_genes)
exhaustion_genes

# Compute exhaustion module score
seu <- AddModuleScore(
  object   = seu,
  features = list(exhaustion_genes),
  name     = "ExhaustionScore",
  assay    = DefaultAssay(seu)
)

# This creates a column "ExhaustionScore1" in meta.data
head(seu@meta.data[, c("RloopScore1", "ExhaustionScore1")])

DefaultAssay(seu) <- "RNA"

df_tnk <- FetchData(
  seu,
  vars = c("CD8A", "CD8B", "RloopScore1", "ExhaustionScore1",
           "sample_id", "HPV_status", "celltype_final")
)

# Keep only T/NK cluster
df_tnk <- df_tnk %>%
  filter(celltype_final == "T/NK") 

# Define CD8 T cells
cd8_threshold <- 0.5  # tweak if needed
df_cd8 <- df_tnk %>%
  mutate(
    is_CD8 = (CD8A > cd8_threshold | CD8B > cd8_threshold)
  ) %>%
  filter(is_CD8)

# Quick check
table(df_cd8$HPV_status)
nrow(df_cd8)


# scatter exhaustion ------------------------------------------------------


# Optional: correlation
cor_cd8_spearman <- cor(
  df_cd8$RloopScore1,
  df_cd8$ExhaustionScore1,
  method = "spearman",
  use    = "complete.obs"
)
cor_cd8_spearman

# Scatter plot
p_cd8_scatter <- ggplot(
  df_cd8,
  aes(x = RloopScore1,
      y = ExhaustionScore1,
      color = HPV_status)
) +
  geom_point(alpha = 0.4, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  labs(
    x = "R-loop module score",
    y = "Exhaustion module score",
    color = "HPV status",
    title = "CD8 T cells (T/NK cluster): R-loop vs exhaustion score"
  ) +
  theme_bw()

p_cd8_scatter

ggsave(
  filename = file.path(res_dir, "GSE171894_CD8_scatter_Rloop_vs_Exhaustion.png"),
  plot     = p_cd8_scatter,
  width    = 6,
  height   = 5,
  dpi      = 300
)


# violin exhaustion -------------------------------------------------------

# Median split on R-loop score
med_Rloop_cd8 <- median(df_cd8$RloopScore1, na.rm = TRUE)

df_cd8 <- df_cd8 %>%
  mutate(
    Rloop_group = if_else(
      RloopScore1 >= med_Rloop_cd8,
      "Rloop-high",
      "Rloop-low"
    ),
    Rloop_group = factor(Rloop_group, levels = c("Rloop-low", "Rloop-high"))
  )

# Violin of exhaustion score by R-loop group
p_cd8_violin <- ggplot(
  df_cd8,
  aes(x = Rloop_group,
      y = ExhaustionScore1,
      fill = Rloop_group)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.size = 0.3) +
  labs(
    x = "R-loop group (CD8 T cells)",
    y = "Exhaustion module score",
    title = "CD8 T cells: exhaustion in R-loop-high vs R-loop-low"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

p_cd8_violin

ggsave(
  filename = file.path(res_dir, "GSE171894_CD8_violin_Exhaustion_by_Rloop_group.png"),
  plot     = p_cd8_violin,
  width    = 4,
  height   = 5,
  dpi      = 300
)



# anitgens ----------------------------------------------------------------


# Core MHC I genes (edit as needed)
mhc1_genes <- c(
  "HLA-A", "HLA-B", "HLA-C",
  "B2M", "TAP1", "TAP2",
  "PSMB8", "PSMB9"
)

# Keep only genes actually present in your object
mhc1_genes <- mhc1_genes[mhc1_genes %in% rownames(seu)]
mhc1_genes

malignant_labels <- c(
  "Epithelial",
  "Cycling tumor/epithelial",
  "Basal epithelial/tumor"
)

seu_epi <- subset(
  seu,
  subset = celltype_final %in% malignant_labels
)

# Sanity check:
table(seu_epi$celltype_final, seu_epi$HPV_status)


seu_epi <- AddModuleScore(
  object   = seu_epi,
  features = list(mhc1_genes),
  name     = "MHC1"
)

Idents(seu_epi) <- "HPV_status"

p_mhc1_hpv_simple <- VlnPlot(
  seu_epi,
  features = "MHC11",   # the module score
  pt.size  = 0
) +
  theme_bw(base_size = 10) +
  labs(
    x = "HPV status",
    y = "MHC I module score",
    title = "MHC I expression by HPV status"
  ) +
  NoLegend()

ggsave(
  filename = file.path(res_dir, "GSE171894_MHC1_module_by_HPV.png"),
  plot     = p_mhc1_hpv_simple,
  width    = 4,   # much smaller
  height   = 3,
  dpi      = 300
)


Idents(seu_epi) <- "Rloop_group"

p_mhc1_rloop_simple <- VlnPlot(
  seu_epi,
  features = "MHC11",
  pt.size  = 0
) +
  theme_bw(base_size = 10) +
  labs(
    x = "R-loop group",
    y = "MHC I module score",
    title = "MHC I expression by R-loop group"
  ) +
  NoLegend()

ggsave(
  filename = file.path(res_dir, "GSE171894_MHC1_module_by_RloopGroup.png"),
  plot     = p_mhc1_rloop_simple,
  width    = 4,
  height   = 3,
  dpi      = 300
)

