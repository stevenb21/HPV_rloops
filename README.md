# HPV R-loops (scRNA-seq analysis)

Research pipeline to test whether **HPV status** is associated with **R-loop–linked transcriptional programs** in single-cell data.

**Scope note:** the research question and the repo’s implementation are orthogonal—this is see-through, stepwise analysis code (scripts + saved Seurat objects), not a polished package. The dataset used here ultimately lacked key metadata/structure needed to cleanly answer the central question, but the analysis design and scoring framework are reusable.

An initial attempt to extend Chen et al.’s HBV-focused approach to HPV.

![fig](./Rationale.png)

Includes QC → integration (Harmony) → clustering/UMAP → marker-based annotation → R-loop gene-set scoring/plots.

**Related prior work:** earlier R-loop analyses supporting *Moffitt et al.* are at **[ZPR1_TCGA](https://github.com/stevenb21/ZPR1_TCGA)**.

**Future work:** incorporate **viral integration status** into the R-loop score/models using a head & neck cancer dataset.
