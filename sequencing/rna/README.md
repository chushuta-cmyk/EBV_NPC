# EBV RNA analysis

This folder contains a modular version of the original RNA analysis. The original script is kept in `legacy/EBV_resequencing_legacy.R`.

## Modules

The modules run in this order:

1. `data` loads counts, selects samples, runs DESeq2, and creates rlog data.
2. `figures` creates normalization, PCA, volcano, and heatmap figures.
3. `enrichment` runs GO and KEGG enrichment.
4. `immune` loads or runs CIBERSORT.
5. `hub` compares hub genes with immune cells.
6. `ml` compares classification models.
7. `wgcna` builds a co-expression network.
8. `pathways` compares immune and EBV pathway scores with immune cells.
9. `survival` runs Cox models and Kaplan-Meier plots.

The runner adds required earlier modules automatically.

## Quick start

Run the core analysis with a local count file:

```bash
RNA_COUNTS_FILE=/path/to/GSE224450_raw_counts.csv \
Rscript run_analysis.R --modules=core
```

Run selected modules:

```bash
Rscript run_analysis.R --modules=data,enrichment,wgcna
```

Run every module:

```bash
RNA_COUNTS_FILE=/path/to/GSE224450_raw_counts.csv \
RNA_CIBERSORT_FILE=/path/to/CIBERSORT_Results.csv \
RNA_SURVIVAL_PHENOTYPE_FILE=/path/to/GSE102349_samples.json \
RNA_SURVIVAL_EXPRESSION_FILE=/path/to/GSE102349_NPC_mRNA_processed.txt.gz \
Rscript run_analysis.R --modules=all
```

Without `RNA_COUNTS_FILE`, the data module downloads GSE224450. Without `RNA_SURVIVAL_PHENOTYPE_FILE`, the survival module downloads GSE102349 phenotype data.

Use `RNA_OUTPUT_DIR` to change where results are saved. The default is `output/` in this folder.

## Checks

List modules without running them:

```bash
Rscript run_analysis.R --list-modules
```

Check a planned run:

```bash
Rscript run_analysis.R --modules=hub --dry-run
```

Run the local tests:

```bash
Rscript tests/test_pipeline.R
```

Each module reports missing R packages before it starts. Install only the packages needed for the modules you plan to run.
