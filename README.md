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
  qPCR/                qPCR figures; raw Ct calculation code is not archived
sequencing/
  rna/                 Bulk transcriptomics pipeline (setup, data, figures, enrichment,
                       immune deconvolution, hub genes, machine learning, WGCNA,
                       pathways, survival)
  machine-learning/    Classification and gene-immune correlation analyses
  single-cell/         NPC43 analysis; downstream state-marker, enrichment, pseudotime, and LIANA code is in downstream_v3/
  figure_redraw/       Figure-only scripts and saved source data; see its audit for missing analyses
  survival analysis/   Survival analyses in GSE102349
```

## Running the pipelines

From the repository root, all scripts locate the project root themselves (via the `--file=` argument or `RNA_PROJECT_DIR`).

- **Bulk RNA (sequencing/rna):** run the full pipeline from that directory:
  ```bash
  cd sequencing/rna
  RNA_INPUT_DIR=/path/to/inputs Rscript run_full_pipeline.R
  ```
  `run_analysis.R` runs selected modules (`--modules=data,figures,enrichment`); `run_sensitivity_analysis.R` (needs `RNA_SERIES_MATRIX` or the file inside the input dir) and `run_ppi_analysis.R` run the extra analyses; `tests/validate_output.R` checks outputs. Modules in `R/` run in order `00_setup.R` → `01_data.R` → ... → `09_survival.R`, plus `10_sensitivity.R` and `11_ppi.R`.
- **Input data:** GEO raw/normalized files live outside this repository. Set `RNA_INPUT_DIR` (required; a local `data/inputs/` inside the repo is the fallback and an error is raised if missing), `RNA_OUTPUT_DIR` (default `<repo>/output`), and optionally `RNA_EXTERNAL_LIBRARY`.
- **Machine learning (sequencing/machine-learning):** `run_ml_analysis.py` (also available as `EBV_ML.ipynb`).
- **Single-cell (sequencing/single-cell):** `run_single_cell_jupyterhub_fixed.py` (also available as `sc_umap_fixed_vscode (1).ipynb`).
- **Experiments (experiments/):** `barplot_cytokines.R` / `heatmap_cytokines.R` (cytokines), `analyze_new_cytokines.R` (co-culture panel), `protemic_graph.py` (proteomics), `final_qpcr.r` (qPCR).

## Contact

Corresponding author: Xiaoming Lyu, PhD — xiaomlyu@smu.edu.cn
