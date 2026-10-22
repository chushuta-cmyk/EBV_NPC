# NPC43 downstream single-cell analysis

This folder archives the recovered version 3 downstream analysis script and its saved tables. The script was found in the local `NPC43_results_with_downstream` work folder and copied here unchanged.

## Analyses in the script

- Rank cell-state markers with Scanpy's Wilcoxon test and Benjamini–Hochberg correction.
- Run Enrichr with GO Biological Process 2023, KEGG 2021 Human, and Reactome 2022 libraries.
- Test gene–pseudotime associations with Spearman correlation and Benjamini–Hochberg correction.
- Optionally run LIANA rank-aggregate analysis among untreated states with at least 30 cells.
- Summarize cell-state proportions and map bulk/proteomic candidate genes across states.

These are exploratory cell-level summaries. The script states that individual cells are not independent biological replicates and does not treat the cell-level tests as confirmatory treatment inference.

## Inputs and running

The script reads outputs from `run_single_cell_jupyterhub_fixed.py`; it does not rerun QC, normalization, clustering, UMAP, or Palantir. The required `--run-dir` must contain:

- `adata_all_cells_hvg_fixed.h5ad`
- `adata_ut_full_genes_fixed.h5ad`
- `input_validation.json`
- `run_audit.json`

The H5AD files are not copied into this code repository. The recovered originals remain under `~/Downloads/NPC43_results_with_downstream/results/03_full_with_palantir/` on this computer. A fresh run should write to a new output folder:

```bash
python run_npc43_downstream_clean_v3.py \
  --run-dir /path/to/03_full_with_palantir \
  --out-dir /path/to/new_downstream_output \
  --with-liana
```

Enrichr runs by default and uses an online service. `--skip-enrichment` disables it. `--with-liana` enables LIANA. The saved Enrichr and LIANA results were generated in separate archived output folders (`results/05_downstream_enrichment_v3/` and `results/06_downstream_liana_v3/`).

The recorded environment is in each run's `audit/software_environment.txt`. It records Python 3.10.20, Scanpy 1.11.5, GSEApy 1.3.1, and LIANA 1.10.0, among other packages. `statsmodels` is imported for multiple-testing correction but its version was not recorded, so the environment record is not a complete lockfile.

No analysis was rerun during this archive update. The original result tables and audit scope files were copied without content changes; see `RUN_PROVENANCE.md`.
