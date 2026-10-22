# Manuscript-to-code archive audit

## Scope

This audit compares the 17 October `EBV_manuscript_JTM_revised.docx` draft with the code and outputs located in the EBV project and relevant Codex work folders. It is a static inventory. The current validation log reports that the raw-RCC data module passed on 17 October 2026, while enrichment, immune analysis, machine learning, WGCNA, survival, sensitivity analysis, and PPI were not run in that validation round (`../../TEST_REPORT.md`).

## Coverage

| Manuscript analysis | Code or saved evidence | Archive status |
|---|---|---|
| Raw RCC processing, housekeeping size factors, DESeq2, post-DE rlog | `../rna/R/01_data.R`; the Layer 2 comparison passed for 22 samples, 730 probes, and 78 strict DE transcripts | Code present; Layer 2 verified. The draft states DESeq2 1.40.2, but the validation log records DESeq2 1.50.2 in the active library. Confirm the version used for the manuscript and correct the Methods. |
| GO and KEGG enrichment | `../rna/R/03_enrichment.R`; saved bulk enrichment tables and figure inputs are available | Code present; not rerun in the current validation. KEGG requires network access. |
| STRING protein interaction network | `../rna/R/11_ppi.R`, `../rna/run_ppi_analysis.R`; saved `source_data/ppi_edges.csv` and `source_data/ppi_node_centrality.csv` | Code and saved figure inputs present; not rerun in the current validation. STRING access requires a network connection. |
| Immune deconvolution and gene–cell correlations | `../rna/R/04_immune.R` and the Python ML workflow | Code present; not rerun in the current validation. The draft's result values still need to be tied to a named output run. |
| Seven-classifier LOOCV and WGCNA | `../machine-learning/run_ml_analysis.py`; `../rna/R/07_wgcna.R` | Code present; not rerun in the current validation. The Python script defines the seven classifiers named in the draft. |
| Survival analysis | `../rna/R/09_survival.R` and `../survival analysis/` | Code present; not rerun in the current validation. |
| Stage and sample-exclusion sensitivity analyses | `../rna/R/10_sensitivity.R` and `../rna/run_sensitivity_analysis.R` | Code present; not rerun in the current validation. The validation log says the series matrix is outside the configured input directory and must be supplied through `RNA_SERIES_MATRIX`. |
| Single-cell QC, embedding, and untreated-cell scoring | `../single-cell/run_single_cell_jupyterhub_fixed.py`; original notebook issues are described in `SINGLE_CELL_NOTEBOOK_AUDIT.md` | Core processing code is present, but this script does not reproduce all downstream methods claimed in draft §2.9. |
| Single-cell state-marker Wilcoxon tests, Enrichr, pseudotime Spearman tests, and LIANA rank-aggregate analysis | `../single-cell/downstream_v3/run_npc43_downstream_clean_v3.py`; staged results are archived under `../single-cell/downstream_v3/results/` | Code and result tables found in the Downloads work folder and now archived. The script uses preprocessed H5AD files; those large inputs remain outside the code repository. The run record does not preserve the exact original command line. |
| Proteomics differential testing | `../../experiments/protemic/protemic_graph.py` reads a precomputed result-summary workbook and makes plots | Plotting code is present; the upstream DIA abundance-testing workflow is not in this code archive. |
| qPCR analysis | `../../experiments/qPCR/final_qpcr.r`; three older plotting scripts are preserved under `../../experiments/qPCR/legacy/` | The scripts plot fixed RQ values; they do not calculate ΔΔCt from raw Ct measurements. One legacy script adds random noise to demonstration ΔCt values, so that plot is synthetic and must not be treated as experimental data. No Ct-level analysis script was located in the code repository. |

## Manuscript issues to resolve

1. **Single-cell denominator conflict.** Methods §2.9 reports 17,329 cells and 27,590 genes after QC, which matches the archived downstream input audit. Results §3.9 reports the raw pre-QC condition counts (5,121 UT, 2,571 T24, and 10,049 T48; total 17,741). The saved post-QC state table gives 4,951 UT, 2,530 T24, and 9,848 T48 cells (total 17,329). The state percentages in the draft match the post-QC table. Replace the Results counts with the post-QC totals so its denominators match the Methods and plotted percentages.
2. **Single-cell analysis code recovered.** The downstream v3 script implements the state-marker Wilcoxon tests, Enrichr libraries, gene–pseudotime Spearman tests, and optional LIANA rank aggregation described in §2.9. The relevant result tables and run-scope files are now archived under `../single-cell/downstream_v3/results/`. The script consumes preprocessed H5AD objects that remain outside the code repository, and the original command lines were not preserved. The current LIANA result table and top-interaction CSV are included; a separate script for deriving the manuscript's final unique-pair shortlist was not found, so document that selection rule if it differs from the archived CSV.
3. **Interpretation boundary.** The Results describe NPC43 expression as “direct evidence” for the source of the patient Low_CP signature. The Discussion's limitations correctly state that the cell-line analysis does not validate the clinical plasma-EBV grouping. Keep the claim at the level of expression detected in the NPC43 model; it does not establish the cellular source of the patient-tumor signal.
4. **Software version.** The DESeq2 version in §2.3 does not match the active version recorded in the current validation log. Verify the version associated with the reported result tables before submission.
5. **Current rerun scope.** The bulk validation log verifies only the raw-RCC data module. Results for bulk enrichment, immune deconvolution, ML, WGCNA, survival, sensitivity, and PPI remain historical outputs until their modules are run to fresh, isolated output directories and compared with the manuscript tables and figures. The separate single-cell downstream outputs were run previously and are archived here; they were not rerun in the current bulk validation round.

## Search limits

The search covered the EBV code repository, relevant Codex work folders, the 17 October submission folder, and `Downloads/NPC43_results_with_downstream/`. The downstream v3 script and its outputs were found in Downloads and copied into the repository. The figures and displayed result tables were treated as outputs; only the recovered script's documented steps are treated as rerunnable calculations.

## Audit guidance

The evidence-tracing approach was informed by Kassis T, Agarwal V, He Y, Patel D, Brueckner AM. *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. 2026. doi: [10.48550/arXiv.2609.00065](https://doi.org/10.48550/arXiv.2609.00065). This citation describes the audit guidance; it is not evidence for any study result.
