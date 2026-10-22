# Archived downstream run provenance

The source outputs were saved in two downstream run folders generated from the same finalized NPC43 objects:

- `results/05_downstream_enrichment_v3/` contains state-marker, Enrichr, and pseudotime-association outputs. Its scope file records a runtime of 11.24 minutes.
- `results/06_downstream_liana_v3/` contains the LIANA result tables and candidate interaction table. Its scope file records a runtime of 1.38 minutes.

The original run folders recorded the following input dimensions:

| Data stage | Cells | Genes | Evidence |
|---|---:|---:|---|
| Raw matrix after barcode matching | 17,741 | 27,595 | 17,741/17,741 barcodes matched |
| All-cell object after QC | 17,329 | 2,000 HVGs | Palantir and UMAP fields present |
| Untreated full-gene object | 4,951 | 27,590 | Full filtered gene matrix and all five scores present |

The raw condition totals were UT 5,121, T24 2,571, and T48 10,049. After QC, the saved state-count table gives UT 4,951, T24 2,530, and T48 9,848, totaling 17,329. The manuscript Methods §2.9 reports the post-QC cell and gene totals, while Results §3.9 reports the pre-QC condition totals. The state-composition percentages in the draft match the post-QC output. Update the Results sample counts to the post-QC totals and retain the matching percentages.

The upstream audit records two UT batches, one T24 batch, and two T48 batches. It also states that Palantir completed and that treatment-level differential expression was not tested because T24 has one batch. The downstream scope files describe state-marker tests, enrichment, pseudotime, and LIANA as exploratory.

The exact command lines and upstream H5AD input objects are not included in this archive. The command example in `README.md` describes how to rerun the workflow from the named inputs; it is not a record of the original invocation.
