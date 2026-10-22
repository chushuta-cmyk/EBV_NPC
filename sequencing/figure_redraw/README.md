# Figure redraw archive

This folder preserves the scripts and saved inputs used to redraw manuscript Figures 3–9. These scripts replot existing results; they do not rerun differential expression, enrichment, deconvolution, pseudotime, or ligand–receptor analysis.

## Files

- `redraw_figures.py` redraws Figures 3–7 and 9 from `source_data/`.
- `redraw_singlecell.py` redraws Figure 8 from saved single-cell coordinates and result tables.
- `source_data/` contains 30 copied source and derived files, including the pseudotime gene tables, cell-state expression tables, PPI edge and node tables, and figure inputs.
- `source_data/provenance.json` records the source package and input roles without machine-specific absolute paths.
- `SINGLE_CELL_NOTEBOOK_AUDIT.md` preserves the earlier notebook and rerun audit.
- `MANUSCRIPT_CODE_AUDIT.md` maps the current manuscript methods to code and records unresolved gaps.
- The recovered state-marker, Enrichr, pseudotime, and LIANA calculations are archived in `../single-cell/downstream_v3/` with their saved result tables.

The two original source manifests contained personal absolute paths. They were not copied; their input roles are summarized in `source_data/provenance.json`.

## Replotting

Run these commands from this folder after installing the packages in `requirements.txt`:

```bash
python redraw_figures.py
python redraw_singlecell.py
```

The first command draws Figures 3–7 and 9 by default. The scripts save PDF and PNG files here. No redraw was run as part of this archive update.

## Scope

The saved tables support figure reconstruction and review. They do not replace the missing calculation scripts listed in `MANUSCRIPT_CODE_AUDIT.md`, and they do not establish that the manuscript's statistical results were reproduced by the current pipeline.
