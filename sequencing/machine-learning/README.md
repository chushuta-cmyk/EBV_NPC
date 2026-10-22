# EBV machine-learning analysis

`run_ml_analysis.py` is the checked command-line version of the notebook analysis. It reads the newly recalculated CIBERSORT results and raw-RCC-derived hub-gene expression by default. It keeps preprocessing inside each validation fold and calculates balanced accuracy after all leave-one-out predictions are collected.

Run it with the default project data:

```bash
python3 run_ml_analysis.py --jobs 1
```

Use other inputs when needed:

```bash
python3 run_ml_analysis.py \
  --cibersort /path/to/cibersort_results.csv \
  --genes /path/to/hub_gene_expression.csv \
  --output-dir /path/to/output
```

The original notebook used a stale hard-coded input and averaged balanced accuracy inside one-sample test folds. This command-line version fixes both issues. The data contain only 22 samples, including 5 Low_CP samples. Treat model rankings, correlations, and feature weights as exploratory results that require external validation.
