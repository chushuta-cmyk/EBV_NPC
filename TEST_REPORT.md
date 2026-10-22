# TEST_REPORT — EBV_NPC analysis validation (work in progress)

Status legend: `passed` / `failed` / `not run` / `blocked`.

**Overall status: Layer 1 (pre-run) and Layer 2 (data module) passed; enrichment/immune/ML/WGCNA/survival/sensitivity/PPI not run.**
This report is the single validation log for the portabilized code. It is not a
claim that the paper results have been reproduced.

## Layer 1 — Pre-run checks (completed)

| Check | Result | Notes |
|---|---|---|
| Baseline (paper-result) output identified | passed | `/Users/tan/Documents/Codex/2026-09-14/cou/outputs/nanostring_raw_rcc_full_pipeline/output/` — full run finished 2025-09-15 05:16–05:21 (tables, figures, cache, `session_info.txt`). Baseline chosen by "the run that produced the paper numbers" (730 probes / 78 DE / 11 up / 67 down per its `README_运行说明.md`), not by newest mtime. |
| Current code == baseline code + path-only edits | passed | Files byte-identical to the latest pipeline source (17/17 R, py, tests, READMEs); only path/config defaults changed. No statistical logic, thresholds or parameters modified. |
| Input directory exists & read-only in practice | passed | `.../work/public_reanalysis/inputs/` contains all required files (GSE224450_RAW.tar, GSE102349 samples/expression, GSE224450 raw counts/probes/samples, technical metadata). RCC extraction uses `tempfile()` and is deleted on exit — inputs are never written to. |
| Input checksums recorded | passed | See "Input checksums" below (SHA-256, recorded 2026-10-17). |
| `GSE224450_series_matrix.txt.gz` (sensitivity) | blocked | Not in the input dir; exists at the old Codex path `/Users/tan/Documents/Codex/2026-09-14/cou/work/GEO/...` (SHA-256 `b774a624…`). Sensitivity runs need `RNA_SERIES_MATRIX` set to it (or the file placed in the input dir). |
| Env-var routing (INPUT/OUTPUT/SERIES_MATRIX) | passed | `run_full_pipeline.R` sets `RNA_PROJECT_DIR/MODULES/INPUT_DIR/SURVIVAL_*`; `make_config()` reads `RNA_INPUT_DIR/RNA_OUTPUT_DIR/RNA_EXTERNAL_LIBRARY/RNA_CIBERSORT_PERM/RNA_SERIES_MATRIX`; `10_sensitivity.R` defaults `RNA_SERIES_MATRIX` to `<input_dir>/GSE224450_series_matrix.txt.gz`; `validate_output.R` reads `RNA_OUTPUT_DIR`. All entry points route through `config$tables_dir/figures_dir/cache_dir`. |
| Output isolation | passed | No module writes outside `config$*_dir`. New test output will go to a fresh directory (e.g. `/tmp/ebv_npc_test_out`), disjoint from inputs and from the baseline `output/`. |
| Auto-install / download / side effects on load | passed (with caveats) | No `install.packages` / `download.file` / `system()` in `sequencing/rna/` or `sequencing/machine-learning/`. **Caveats:** `legacy/EBV_resequencing_legacy.R` contains `install.packages` calls — not part of the paper path, not run, not validated. `03_enrichment.R` calls `clusterProfiler::enrichKEGG` (KEGG REST, **network**; no local BiocFileCache found) and `11_ppi.R` uses `STRINGdb` (STRING server, **network**; no `~/.STRINGdb` cache found). |
| R version | passed | R 4.5.1 (2025-06-13) — identical to baseline `session_info.txt`. |
| Package versions — default library vs baseline | **flagged** | Default library is NEWER Bioconductor than baseline: edgeR 4.8.2 vs 4.6.3, limma 3.66.0 vs 3.64.3, GSVA 2.4.9 vs 2.2.1, DOSE 4.4.0 vs 4.2.0, BiocFileCache 3.0.0 vs 2.16.2, KEGGREST 1.50.0 vs 1.48.1, preprocessCore 1.72.0 vs 1.70.0, DESeq2 1.50.2 vs 1.48.2 (ext), clusterProfiler 4.18.4 vs 4.16.0 (ext), org.Hs.eg.db 3.22.0 vs 3.21.0 (ext). |
| Package versions — external library (baseline lib) | passed | `/Users/tan/Developer/envs/R/lib/R/library` matches baseline for the packages it provides: immunedeconv 2.1.0, WGCNA 1.74, caret 7.0-1, ggplot2 4.0.3, preprocessCore 1.70.0, genefilter 1.90.0, org.Hs.eg.db 3.21.0, EPIC 1.1.7, MCPcounter 1.2.0, quantiseqr 1.16.0, xCell 1.1.0, CIBERSORT 0.1.0, xgboost 3.2.1.1, kernlab 0.9-33, DESeq2 1.48.2, clusterProfiler 4.16.0. **Missing from default lib: immunedeconv, xgboost, kernlab** → must set `RNA_EXTERNAL_LIBRARY` for the immune/ML modules. |
| Consequence of removing the default external-library path | **flagged** | Baseline run mounted the external lib by default. After portabilization it only mounts when `RNA_EXTERNAL_LIBRARY` is set. Without it, R loads newer Bioconductor defaults → results may differ from baseline. Regression runs must set `RNA_EXTERNAL_LIBRARY=/Users/tan/Developer/envs/R/lib/R/library`; residual version drift (edgeR/limma/GSVA/DOSE/BiocFileCache/KEGGREST from default lib) cannot be removed without downgrading — results are compared with tolerance, not byte-for-byte. |
| Python environment | flagged | Default `python3` = 3.14.7 (no sklearn). Candidate `/usr/local/bin/python3` has sklearn 1.7.2 / numpy 2.3.3 / pandas 2.3.3. Baseline ML environment was not recorded → ML re-run is a fresh run, not a byte-identical check. |
| `validate_output.R` scope | passed | Checks 13 hard-coded assertions by reading current `RNA_OUTPUT_DIR` (module status, 730 probes, 78 DE, 11/67 directions, integer counts, no rounding, ≥4 immune methods, GO non-empty, KEGG background ≤730, finite ML ROC, WGCNA/survival non-empty, non-empty PNGs). **Limitation:** it does not verify timestamps/provenance — if pointed at the baseline or an old output dir it can pass on stale files. It must only be run against this round's fresh output directory. |
| No credentials/sensitive content | passed | Scans of worktree and of unpushed commit `9bdaa20`: no keys/tokens/passwords. Notebook embeds 12 UMAP PNGs (normal outputs). |
| No commit / no push | passed | Worktree has uncommitted changes only; `9bdaa20` preserved; nothing pushed. |

## Input checksums (SHA-256, recorded 2026-10-17)

| File | SHA-256 |
|---|---|
| GSE224450_RAW.tar | `06902fa3d532d87d546f830a0f1a531ee49db16b4a25d5dd9912ad6ea06361fd` |
| GSE102349_samples.json | `ca0ee315b27ac13c35ca14e1510c4c0f623de4ca4d9c93fc30c2e57c1a215bf7` |
| GSE102349_NPC_mRNA_processed.txt.gz | `b85a16f7bc6d3441aec047cf5015b9a41633ec56822499c01495116716845773` |
| GSE224450_raw_counts.csv | `e176ef38adccb23e62d5e5d7ebeee6a83fd08c856dc70e8353b9936e9b83989f` |
| GSE224450_probes.csv | `26c9097fa2249ad9e9fe4797df902acee1ca60d28d1abd38973725ed1a2d3b05` |
| GSE224450_samples.json | `9cfd1afc04650861efbdab694a8a0bc6e02ad7003f723e9bcdd5886f4e20ecb8` |
| GSE224450_series_matrix.txt.gz (old Codex path) | `b774a6249349768c03adffb756ec7a4aae966656cce214558b01b97839cadf0d` |

## Layer 2 — data-module smoke/regression run (COMPLETED 2026-10-17, approved scope)

**Result: passed.** Command, logs, comparison, and outputs persisted in
`/Users/tan/Desktop/EBV_publication/validation_run_20261017/` (`L2_command.txt`,
`L2_run.log`, `L2_output/session_info.txt`, `L2_comparison.txt`, `compare_l2.py`).

Command (exit 0):
`RNA_INPUT_DIR=<inputs> RNA_OUTPUT_DIR=<.../validation_run_20261017/L2_output> RNA_EXTERNAL_LIBRARY=/Users/tan/Developer/envs/R/lib/R/library RNA_SERIES_MATRIX=<old Codex path> Rscript run_analysis.R --modules=data`

What ran (per `R/01_data.R`): raw RCC read from `GSE224450_RAW.tar` (22 samples,
no rounding, integer check) → endogenous filter `rowSums(≥10) ≥ 2` → DESeq2 with
housekeeping-only size factors, `design = ~Group`, `contrast = c("Group","High_CP","Low_CP")`
→ rlog (after DE). This is a read/preprocess/differential-analysis reproduction;
no enrichment/immune/ML/WGCNA/survival/sensitivity/PPI, no network, no commit/push.

Environment actually loaded: `.libPaths()` = [default R 4.5.1 library (priority),
external baseline library]. DESeq2 1.50.2 / edgeR 4.8.2 / SummarizedExperiment 1.40.0
(default-lib, newer than baseline — see Layer 1 flag). External lib provides
immunedeconv/xgboost/kernlab etc. for later layers. `R_library` inside the repo
does not exist, so no override.

Comparison vs baseline `output/tables/`:

| Check | New vs baseline | Result |
|---|---|---|
| sample_groups (22 IDs, encoding) | identical | passed |
| Low_CP / High_CP | 5 / 17 | passed |
| probe_audit (8 metrics) | identical | passed |
| DE rows / gene set (730) | identical | passed |
| Strict DE (padj<0.01, \|log2FC\|>2) | 78 = 78; up 11 / down 67; gene sets identical (0 diff) | passed |
| log2FoldChange diff | 382 rows differ, max abs 1.07e-14 (LILRA5) | passed (rounding-level) |
| pvalue diff | 528 rows differ, max abs 1.3e-14 (BAX) | passed (rounding-level) |
| padj diff | 548 rows differ, max abs 1.81e-14 (NFKB1) | passed (rounding-level) |
| filtered/rlog/normalized matrices | 730×23 = 730×23 | passed |

All numeric differences are floating-point rounding (<2e-14); no gene crosses any
threshold differently. Outputs confirmed freshly generated 2026-10-17 19:02 in the
new isolated dir (no baseline/cache reuse). No tolerance loosening, no parameter or
sample changes.

## Planned test layers (awaiting approval)

- **Layer 3 — regression vs baseline:** enrichment, immune, ML, WGCNA, survival, sensitivity, PPI, each compared to baseline `output/tables/...` with tolerance. Network steps: enrichment (KEGG), PPI (STRINGdb). ML via `/usr/local/bin/python3 run_ml_analysis.py`.
- **Layer 4 — full run:** `run_full_pipeline.R` once in a fresh output dir; key tables/figures cross-checked against the manuscript.

On failure: record the exact error and the completed stage; no parameter/sample/method changes, no package installs, no silent retries.
