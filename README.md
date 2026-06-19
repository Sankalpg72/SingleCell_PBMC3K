# PBMC3K Single-Cell RNA-seq: Preprocessing Perturbation Study

## Intent

This repository documents a systematic investigation into how small, deliberate changes in preprocessing decisions affect clustering outcomes and marker gene identification in the **PBMC3K dataset** (Peripheral Blood Mononuclear Cells, 3000 cells, 10x Genomics).

The core question driving this work is deceptively simple:

> *How much does a single preprocessing choice change the biology you recover?*

Each run in this repository changes **one variable at a time** relative to the previous run. The goal is not to produce a final, definitive analysis — it is to make the sensitivity of single-cell pipelines visible and traceable. Because every step in the pipeline feeds into the next (HVG selection → PCA → neighbourhood graph → Leiden clustering → marker genes), even a small upstream change can compound into dramatically different cell type recovery downstream.

---

## Dataset

**PBMC3K** — 2,700 peripheral blood mononuclear cells from a healthy donor, sequenced using the 10x Genomics Chromium platform. Freely available from [10x Genomics](https://www.10xgenomics.com/datasets).

Raw data lives in `data/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/`.

Expected cell types in this dataset: CD4+ T cells, CD8+ T cells, NK cells, B cells, CD14+ Classical Monocytes, FCGR3A+ Non-Classical Monocytes, Dendritic Cells, and Megakaryocytes (Platelets).

---

## Repository Structure

```
SingleCell_github/
├── data/
│   └── pbmc3k_filtered_gene_bc_matrices/
│       └── filtered_gene_bc_matrices/
│           └── hg19/               # Raw 10x MTX files
├── Figures/
│   ├── Run_1/                      # QC plots, UMAP, clustering outputs for Run 1
│   ├── Run_2/                      # QC plots, UMAP, clustering outputs for Run 2
│   └── Run_3/                      # QC plots, UMAP, clustering outputs for Run 3
├── scripts/                        # Analysis notebooks
├── .gitignore
└── README.md
```

> `Reports/` and `.ipynb_checkpoints/` are excluded from version control via `.gitignore`.

---

## Runs

All three runs use the same base workflow built on [Scanpy](https://scanpy.readthedocs.io/en/stable/). The pipeline follows the general structure of the Scanpy PBMC3K legacy tutorial with deliberate modifications at specific steps. The HVG method is `seurat_v3` with `layer="counts"` (raw counts) across all runs. Leiden clustering uses `resolution=0.7` throughout.

---

### Run 1 — Full Preprocessing (Doublet Removal + Ribosomal Filtering)

**Script:** `scripts/Run_1_Single_cell_github.ipynb`

**What was done:**
- Doublet detection and removal using scVI
- Ribosomal gene filtering (genes starting with `RPS`, `RPL` removed)
- Mitochondrial gene percentage filtering (`pct_MT < 5`)
- Basic cell and gene filtering (`min_genes=200`, `max_genes=2500`, `min_cells=3`)
- Normalisation, log1p, HVG (`seurat_v3`, n=2000, `layer="counts"`)
- Regression on `total_counts`, `pct_MT`, `pct_Ribo`
- Scaling, PCA, neighbours, UMAP, Leiden (res=0.7), Wilcoxon DE

**Clusters recovered:** 8

**Cell types identified:**

| Cluster | Cell Type |
|---------|-----------|
| 0 | CD8+ T cells |
| 1 | CD4+ T cells |
| 2 | Low-quality / Doublets (spurious) |
| 3 | B cells |
| 4 | Classical Monocytes |
| 5 | NK cells |
| 6 | Non-Classical Monocytes |
| 7 | Doublets (29 cells, residual) |

**What was missing:** Platelets, Dendritic cells

**Key observation:** Doublet removal and ribosomal filtering together suppress the platelet transcriptional signal below the clustering threshold. The platelet population (PPBP, PF4, SDPR) does not resolve as its own cluster — these cells are likely absorbed into monocyte-adjacent clusters. A small residual doublet cluster (Cluster 7, n=29) persists, indicating doublet detection was imperfect rather than complete.

---

### Run 2 — Doublet Removal Skipped

**Script:** `scripts/Run_2_Single_cell_github.ipynb`

**Change from Run 1:** Doublet detection and removal step skipped. Everything else identical.

**Clusters recovered:** 8

**Cell types identified:**

| Cluster | Cell Type |
|---------|-----------|
| 0 | CD8+ T cells |
| 1 | B cells |
| 2 | CD4+ T cells |
| 3 | Classical Monocytes |
| 4 | NK cells |
| 5 | Low-quality / Doublets (spurious) |
| 6 | Non-Classical Monocytes |
| 7 | Doublets |

**What was missing:** Platelets, Dendritic cells

**Key observation:** Skipping doublet removal had minimal impact on the major immune populations — CD4+ T cells, CD8+ T cells, NK cells, B cells, and both monocyte subsets all resolved correctly. The doublet cluster persists and is slightly larger than in Run 1. Platelet and Dendritic cell populations remain unrecovered, confirming that doublet removal alone was not responsible for their absence in Run 1. This points to ribosomal gene filtering as the more likely suppressing factor.

---

### Run 3 — Doublet Removal and Ribosomal Filtering Both Skipped

**Script:** `scripts/Run_3_Single_cell_github.ipynb`

**Change from Run 2:** Ribosomal gene filtering also skipped. No doublet removal, no ribo filtering. Followed the Scanpy tutorial structure more closely.

**Clusters recovered:** 7

**Cell types identified:**

"CD8 T CELLS",
    "",
    "",
    "",
    "",
    "",
    "DENDRITIC CELLS",
    "MEGAKARYOCYTES",

| Cluster | Cell Type |
|---------|-----------|
| 0 | CD4 T CELLS |
| 1 | NK CELLS |
| 2 | B CELLS |
| 3 | CD14 + MONOCYTES |
| 4 | F3GR3A+ MONOCYTES |
| 5 | Dendritic Cells |
| 6 | Megakaryocytes |

**What was missing:** CD4+ T cells, CD8+ T cells (merged with NK), Dendritic cells

=

## Core Insight

There is no single universally correct preprocessing pipeline. Every choice involves a trade-off:

- **Ribosomal filtering** helps T cell subsets separate cleanly but suppresses small populations like platelets whose transcriptomes overlap with ribosomal-high signatures
- **Doublet removal** reduces noise but imperfect detection means residual doublet clusters still appear
- **Ribosomal regression** (without filtering) is a middle ground — it dampens ribo variance in PCA without removing platelet-adjacent genes entirely, but may be insufficient alone at lower resolutions

The sensitivity observed here — where a single filtering decision changes which cell types you recover — is inherent to the sequential, compounding structure of scRNA-seq pipelines. A change in gene filtering changes HVG selection, which changes PCA space, which changes the neighbourhood graph, which changes Leiden cluster boundaries, which changes marker gene rankings.

---

## Dependencies

```
scanpy
scvi-tools
pandas
numpy
matplotlib
seaborn
leidenalg
```

---

## Reference

Scanpy PBMC3K legacy tutorial: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html

Wolf et al. (2018) — Scanpy: large-scale single-cell gene expression data analysis. *Genome Biology*.
