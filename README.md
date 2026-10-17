# EBV_NPC

Analysis code for the manuscript **"Plasma EBV DNA Grouping Uncovers a Stage-Independent Myeloid-Inflammatory Signature in Nasopharyngeal Carcinoma"** (submitted to the *Journal of Translational Medicine*).

## Overview

This repository contains the analysis pipelines underlying the manuscript, which integrates public clinical bulk transcriptomics, stage-adjusted sensitivity analysis, independent-cohort survival context, single-cell RNA sequencing, quantitative proteomics, and functional in vitro experiments to characterize the tumor immune microenvironment (TIME) programs associated with pre-treatment plasma EBV DNA groups in nasopharyngeal carcinoma (NPC).

All findings are exploratory and were derived from secondary analyses of publicly available datasets; no new human data were collected for this study.

## Data availability

- **GSE224450** — NanoString PanCancer IO 360 profiles of 22 pre-treatment biopsies from a published phase II trial of nivolumab plus ipilimumab in recurrent/metastatic EBV-associated NPC (Lim et al., Nat Commun 2023;14:2781).
- **GSE102349** — bulk RNA-seq of treatment-naive NPC tumors, used for exploratory survival context.
- **GSE266679** — single-cell RNA-seq of the EBV-positive NPC43 model.
- Experimental (proteomics and in vitro) data: deposited in PRIDE/ProteomeXchange upon acceptance; available from the corresponding author on reasonable request until then.

## Repository structure

```
experiments/
  cytokines/           Luminex 18-plex cytokine analyses (NPC43 monoculture and co-culture)
  protemic/            Targeted proteomics analyses (11-protein panel)
  qPCR/                qPCR for BZLF1 and inflammatory genes
sequencing/
  rna/                 Bulk transcriptomics pipeline (setup, data, figures, enrichment,
                       immune deconvolution, hub genes, machine learning, WGCNA,
                       pathways, survival)
  machine-learning/    Classification and gene-immune correlation analyses
  single-cell/         NPC43 single-cell analysis
  survival analysis/   Survival analyses in GSE102349
```

## Running the pipelines

- **Bulk RNA (sequencing/rna):** modules in `R/` are run in order `00_setup.R` → `01_data.R` → `02_figures.R` → `03_enrichment.R` → `04_immune.R` → `05_hub_genes.R` → `06_machine_learning.R` → `07_wgcna.R` → `08_pathways.R` → `09_survival.R`; `run_analysis.R` is the top-level driver.
- **Machine learning (sequencing/machine-learning):** `run_ml_analysis.py` (also available as `EBV_ML.ipynb`).
- **Single-cell (sequencing/single-cell):** `run_single_cell_jupyterhub_fixed.py` (also available as `sc_umap_fixed_vscode (1).ipynb`).
- **Experiments (experiments/):** `barplot_cytokines.R` / `heatmap_cytokines.R` (cytokines), `analyze_new_cytokines.R` (co-culture panel), `protemic_graph.py` (proteomics), `final_qpcr.r` (qPCR).

## Contact

Corresponding author: Xiaoming Lyu, PhD — xiaomlyu@smu.edu.cn
