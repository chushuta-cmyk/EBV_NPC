# NPC43 single-cell notebook audit

## Files checked

- `code/sequencing/single-cell/sc_umap_fixed_vscode (1).ipynb`
- `single_cell/sc_umap_fixed_vscode (1).ipynb`
- `single_cell/sc_umap_fixed_vscode (1) (1).ipynb`
- `single_cell/singcell_umap.py`
- `single_cell/adata_all_cells_computed.h5ad`
- `single_cell/adata_ut_palantir_completed.h5ad`

The first two notebooks are byte-for-byte identical.

## Findings

1. The original notebook is Colab-specific although its title says VS Code ready. It imports `google.colab`, mounts Google Drive, installs packages with notebook shell commands, and uses `/content/drive/MyDrive/tan`. These cells will not run unchanged on a normal JupyterHub.

2. The original workflow creates a 2,000-highly-variable-gene object and then uses that reduced object for gene-set scoring. Its saved output confirms the problem: only CXCL8, IL1A, and TREM1 were available for the five-gene inflammatory score; no EMT genes were available; only KRT14 was available for the epithelial score; and cell-cycle scoring used only one S-phase and two G2M control genes. Those original scores should not be used.

3. The later notebook corrects the untreated subset by taking it from the full 27,590-gene object. The existing `adata_ut_palantir_completed.h5ad` confirms that this corrected path was run: it contains 4,951 untreated cells, all 27,590 filtered genes, Palantir pseudotime, and complete inflammatory, EMT, epithelial, S, and G2M scores. The pseudotime starts in a cell annotated as `UT Differentiated cells` and spans 0 to 1. These corrected untreated results can be retained without rerunning.

4. The later notebook has execution counts missing for cells 14–18 while some cells still contain saved output. This makes the notebook state unreliable as a provenance record. Cells 16–18 were not executed, and the boxplot cell imports `statannotations` but does not perform any statistical annotation.

5. The standalone Python script assigns `adata.obs = cell_anno.copy()` by row order. This should be replaced by barcode-based reindexing with a complete-match check. The raw files currently contain 17,741 cells and 27,595 genes, and all 17,741 barcodes match the annotation.

6. Differential expression in the notebook uses individual cells as replicates. The data contain two UT batches, one T24 batch, and two T48 batches. A cell-level Wilcoxon test can describe cell distributions but cannot provide treatment-level biological replication, especially for T24. It should not be presented as confirmatory treatment-level inference.

7. Code for the manuscript's GO and LIANA panels is not present in the checked notebook or `singcell_umap.py`. Existing LIANA figures are available, but their source calculations cannot be reproduced from the supplied code. They may be retained only as exploratory results if the original LIANA result table or script is available elsewhere.

## Corrected script

`run_single_cell_jupyterhub_fixed.py` was added to `code/sequencing/single-cell/`. It:

- accepts JupyterHub paths through command-line options;
- accepts compressed or uncompressed GEO files;
- matches annotation by barcode;
- uses 2,000 highly variable genes only for PCA, neighbors, Leiden, and UMAP;
- uses the full filtered gene matrix for marker plots and gene-set scoring;
- saves a compact all-cell HVG object and a full-gene untreated object;
- makes Palantir optional;
- limits threads by default to avoid the local Python/Numba neighbor-search lock;
- does not run treatment differential expression because the batch design does not support confirmatory inference.

Example JupyterHub command:

```bash
python run_single_cell_jupyterhub_fixed.py \
  --data-dir /path/to/GEO \
  --output-dir /path/to/single_cell_fixed \
  --num-threads 1
```

Use `--skip-palantir` on a lower-memory node. Use `--validate-only` to check matrix dimensions and barcode matching without running neighbors or UMAP.

## Local validation

The validation-only run completed successfully:

- raw matrix: 17,741 cells × 27,595 genes;
- barcode matches: 17,741/17,741;
- treatment counts: UT 5,121; T24 2,571; T48 10,049;
- peak memory after raw loading: approximately 2.78 GB.

The full local rerun passed filtering, normalization, HVG selection, and PCA, then stalled in the local Python 3.13/Numba neighbor-search backend after a thread-lock message. It was abandoned rather than forcing a high-memory or unstable rerun. This local environment issue does not invalidate the already completed corrected full-gene untreated object.

## Update on 2026-10-22

The later downstream script `run_npc43_downstream_clean_v3.py` was recovered from the local `Downloads/NPC43_results_with_downstream/code/` folder and archived under `sequencing/single-cell/downstream_v3/`. It implements the exploratory state-marker Wilcoxon ranking, Enrichr analysis, gene–pseudotime Spearman tests, and optional LIANA rank-aggregate analysis. Its saved result tables and run notes were archived with it. The earlier statement above that the downstream GO and LIANA calculation code was not found is superseded by this recovery.
