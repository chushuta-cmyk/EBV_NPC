# EBV machine-learning analysis

`run_ml_analysis.py` is the checked command-line version of the notebook analysis. It keeps preprocessing inside each validation fold and calculates balanced accuracy after all leave-one-out predictions are collected.

Run it with the default project data:

```bash
python3 run_ml_analysis.py
```

Use other inputs when needed:

```bash
python3 run_ml_analysis.py \
  --cibersort /path/to/CIBERSORT_Results_Filtered.csv \
  --genes /path/to/Hub_Genes_Expression.csv \
  --output-dir /path/to/output
```

The filtered data contain very few samples. Treat model rankings and feature weights as exploratory results.
