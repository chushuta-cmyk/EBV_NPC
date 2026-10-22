#!/usr/bin/env python3
"""Clean downstream analysis for finalized NPC43 single-cell outputs.

Version 3 fixes expression_summary so that absent genes are represented as
one-row pandas DataFrames before concatenation. This script reads outputs from
run_single_cell_jupyterhub_fixed.py and does not rerun QC, normalization, HVG
selection, PCA, neighbors, Leiden clustering, UMAP, or Palantir.

Required files in --run-dir:
  adata_all_cells_hvg_fixed.h5ad
  adata_ut_full_genes_fixed.h5ad
  input_validation.json
  run_audit.json

Default analyses use real metadata only: Treatment, Batch, cellular_states,
native_cellular_states, leiden, existing scores, and existing pseudotime.
No artificial EBV-high/EBV-low groups are created.

Optional modules:
  --with-treatment-exploratory  Cell-level UT versus T48 ranking, descriptive only.
  --with-liana                  Candidate epithelial-state LIANA analysis, exploratory.
  --skip-enrichment             Skip GO/KEGG/Reactome annotation.
"""

from __future__ import annotations

import argparse
import importlib
import json
import os
import platform
import sys
import textwrap
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Clean downstream analysis for corrected NPC43 single-cell outputs."
    )
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=20260609)
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument("--with-treatment-exploratory", action="store_true")
    parser.add_argument("--with-liana", action="store_true")
    parser.add_argument("--skip-enrichment", action="store_true")
    return parser.parse_args()


# Parse arguments and configure thread variables before NumPy/Scanpy/Numba imports.
ARGS = parse_args()
RUN_DIR = ARGS.run_dir.expanduser().resolve()
OUT_DIR = (ARGS.out_dir or RUN_DIR / "downstream_clean_v3").expanduser().resolve()
FIG_DIR = OUT_DIR / "figures"
TABLE_DIR = OUT_DIR / "tables"
AUDIT_DIR = OUT_DIR / "audit"
ENRICH_DIR = OUT_DIR / "enrichment"
for _directory in (FIG_DIR, TABLE_DIR, AUDIT_DIR, ENRICH_DIR):
    _directory.mkdir(parents=True, exist_ok=True)

os.environ["NUMBA_CACHE_DIR"] = str(AUDIT_DIR / "numba_cache")
os.environ["MPLCONFIGDIR"] = str(AUDIT_DIR / "matplotlib_cache")
os.environ["NUMBA_NUM_THREADS"] = str(ARGS.num_threads)
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"
os.environ["OMP_NUM_THREADS"] = str(ARGS.num_threads)
os.environ["OPENBLAS_NUM_THREADS"] = str(ARGS.num_threads)
os.environ["MKL_NUM_THREADS"] = str(ARGS.num_threads)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(ARGS.num_threads)
os.environ["BLIS_NUM_THREADS"] = str(ARGS.num_threads)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

if ARGS.skip_enrichment:
    gp = None
else:
    try:
        import gseapy as gp
    except ImportError:
        gp = None

if ARGS.with_liana:
    try:
        import liana as li
    except ImportError:
        li = None
else:
    li = None

np.random.seed(ARGS.seed)
sc.settings.verbosity = 2
sc.settings.figdir = str(FIG_DIR)
sc.set_figure_params(dpi=250, facecolor="white", frameon=True)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
sns.set_theme(style="whitegrid", context="notebook")


def load_json(path: Path) -> dict:
    if not path.exists():
        return {"missing_file": str(path)}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception as error:
        return {"read_error": f"{type(error).__name__}: {error}"}


def version_of(import_name: str) -> str:
    try:
        module = importlib.import_module(import_name)
        return str(getattr(module, "__version__", "installed"))
    except Exception as error:
        return f"NOT AVAILABLE: {type(error).__name__}: {error}"


def write_environment_report(path: Path) -> None:
    package_names = [
        "numpy", "pandas", "scipy", "matplotlib", "seaborn", "anndata",
        "scanpy", "numba", "umap", "pynndescent", "igraph", "leidenalg",
        "gseapy", "palantir", "liana",
    ]
    lines = [
        "NPC43 clean downstream workflow version: v3",
        f"Python executable: {sys.executable}",
        f"Python version: {sys.version.replace(chr(10), ' ')}",
        f"Platform: {platform.platform()}",
        "", "Package versions:",
    ]
    lines.extend(f"{name}: {version_of(name)}" for name in package_names)
    lines.extend(["", "Thread environment:"])
    for key in [
        "NUMBA_NUM_THREADS", "NUMBA_THREADING_LAYER", "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "BLIS_NUM_THREADS",
    ]:
        lines.append(f"{key}: {os.environ.get(key, '<not set>')}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def dense_vector(adata, gene: str) -> np.ndarray:
    values = adata[:, gene].X
    return values.toarray().ravel() if sp.issparse(values) else np.asarray(values).ravel()


def ordered_categories(values: pd.Series, preferred: list[str] | None = None) -> list[str]:
    observed = [str(x) for x in values.dropna().unique().tolist()]
    if preferred is None:
        return sorted(observed)
    initial = [x for x in preferred if x in observed]
    remainder = sorted([x for x in observed if x not in initial])
    return initial + remainder


def expression_summary(adata, genes: list[str], groupby: str) -> pd.DataFrame:
    """Summarize detection and log-normalized expression by an existing metadata field.

    Every element appended to rows is a DataFrame. This explicitly handles genes
    that are not represented in the AnnData variable names and prevents a dict/
    DataFrame concatenation error.
    """
    rows: list[pd.DataFrame] = []
    group_values = adata.obs[groupby].astype(str).values
    empty_columns = [
        "gene", "present_in_matrix", "group", "n_cells", "pct_detected",
        "mean_log_normalized_expression", "median_log_normalized_expression",
    ]

    for gene in genes:
        if gene not in adata.var_names:
            absent = pd.DataFrame(
                [{
                    "gene": gene,
                    "present_in_matrix": False,
                    "group": "not_in_matrix",
                    "n_cells": np.nan,
                    "pct_detected": np.nan,
                    "mean_log_normalized_expression": np.nan,
                    "median_log_normalized_expression": np.nan,
                }],
                columns=empty_columns,
            )
            rows.append(absent)
            continue

        values = dense_vector(adata, gene)
        frame = pd.DataFrame({"group": group_values, "expression": values})
        summary = frame.groupby("group", observed=True)["expression"].agg(
            n_cells="size",
            mean_log_normalized_expression="mean",
            median_log_normalized_expression="median",
        )
        summary["pct_detected"] = (
            frame.assign(detected=frame["expression"] > 0)
            .groupby("group", observed=True)["detected"]
            .mean()
            .mul(100)
        )
        summary = summary.reset_index()
        summary.insert(0, "present_in_matrix", True)
        summary.insert(0, "gene", gene)
        rows.append(summary[empty_columns])

    if not rows:
        return pd.DataFrame(columns=empty_columns)
    return pd.concat(rows, ignore_index=True)


def rank_groups_to_dataframe(adata, key: str) -> pd.DataFrame:
    result = adata.uns[key]
    groups = result["names"].dtype.names
    frames = []
    for group in groups:
        frames.append(pd.DataFrame({
            "gene": result["names"][group],
            "scores": result["scores"][group],
            "logfoldchanges": result["logfoldchanges"][group],
            "pvals": result["pvals"][group],
            "pvals_adj": result["pvals_adj"][group],
            "group": group,
        }))
    return pd.concat(frames, ignore_index=True)


def plot_stacked_composition(composition: pd.DataFrame, title: str, outpath: Path) -> None:
    ax = composition.plot(
        kind="bar", stacked=True, figsize=(10, 5.5), colormap="tab20",
        width=0.78, edgecolor="white", linewidth=0.3,
    )
    ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("Cell proportion (%)", fontsize=10)
    ax.tick_params(axis="x", rotation=0)
    ax.legend(
        title="Cellular state", bbox_to_anchor=(1.02, 1), loc="upper left",
        frameon=False, fontsize=8, title_fontsize=9,
    )
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close("all")


def plot_marker_dotplot(adata, genes: list[str], groupby: str, title: str, outpath: Path) -> None:
    found = [gene for gene in genes if gene in adata.var_names]
    if not found:
        print(f"[skip] No requested genes found for: {title}", flush=True)
        return
    dotplot = sc.pl.dotplot(
        adata, var_names=found, groupby=groupby, standard_scale="var",
        dendrogram=False, show=False, return_fig=True,
    )
    dotplot.style(cmap="viridis")
    dotplot.make_figure()
    dotplot.fig.suptitle(title, y=1.02, fontsize=12, fontweight="bold")
    dotplot.fig.savefig(outpath, bbox_inches="tight")
    plt.close("all")


def plot_state_marker_heatmap(adata, key: str, outpath: Path) -> None:
    sc.pl.rank_genes_groups_heatmap(
        adata, key=key, n_genes=8, standard_scale="var", swap_axes=True, show=False
    )
    plt.savefig(outpath, bbox_inches="tight")
    plt.close("all")


def enrich_gene_list(gene_list: list[str], label: str, output: Path, libraries: list[str]):
    if gp is None:
        print("[skip] gseapy unavailable or enrichment was disabled.", flush=True)
        return None
    output.mkdir(parents=True, exist_ok=True)
    genes = sorted(set(gene_list))
    if len(genes) < 5:
        print(f"[skip] {label}: fewer than five genes available.", flush=True)
        return None
    try:
        enrichment = gp.enrichr(
            gene_list=genes,
            gene_sets=libraries,
            organism="human",
            outdir=str(output),
            no_plot=True,
        )
        enrichment.results.to_csv(output / f"{label}_all.csv", index=False)
        fdr = enrichment.results[enrichment.results["Adjusted P-value"] < 0.05].copy()
        fdr.to_csv(output / f"{label}_fdr05.csv", index=False)
        return fdr
    except Exception as error:
        print(f"[warning] Enrichment failed for {label}: {type(error).__name__}: {error}", flush=True)
        return None


def plot_enrichment_table(result: pd.DataFrame | None, title: str, outpath: Path, palette: str) -> None:
    if result is None or result.empty:
        return
    plot_data = result.sort_values("Adjusted P-value").head(10).copy()
    plot_data["minus_log10_fdr"] = -np.log10(plot_data["Adjusted P-value"].clip(lower=1e-300))
    plot_data["term"] = plot_data["Term"].astype(str).map(lambda value: textwrap.fill(value, width=42))
    fig, ax = plt.subplots(figsize=(8, 5.2), dpi=300)
    sns.barplot(data=plot_data, x="minus_log10_fdr", y="term", palette=palette, ax=ax)
    ax.axvline(-np.log10(0.05), color="gray", linestyle="--", linewidth=1)
    ax.set_title(title, fontsize=11, fontweight="bold", pad=10)
    ax.set_xlabel("-log10 adjusted P value")
    ax.set_ylabel("")
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close("all")


def pseudotime_correlations(adata, min_pct_detected: float = 5.0) -> pd.DataFrame:
    if "pseudotime" not in adata.obs.columns:
        return pd.DataFrame()
    pseudotime = adata.obs["pseudotime"].to_numpy(dtype=float)
    if not np.isfinite(pseudotime).all():
        raise ValueError("pseudotime contains non-finite values")
    matrix = adata.X
    if sp.issparse(matrix):
        detected = np.asarray((matrix > 0).mean(axis=0)).ravel() * 100
    else:
        detected = (np.asarray(matrix) > 0).mean(axis=0) * 100
    valid_idx = np.where(detected >= min_pct_detected)[0]
    records = []
    for idx in valid_idx:
        values = matrix[:, idx].toarray().ravel() if sp.issparse(matrix) else np.asarray(matrix[:, idx]).ravel()
        rho, pvalue = spearmanr(values, pseudotime)
        if np.isfinite(rho) and np.isfinite(pvalue):
            records.append({
                "gene": str(adata.var_names[idx]),
                "spearman_rho": float(rho),
                "pvalue": float(pvalue),
                "pct_detected": float(detected[idx]),
            })
    result = pd.DataFrame(records)
    if result.empty:
        return result
    result["qvalue"] = multipletests(result["pvalue"].values, method="fdr_bh")[1]
    return result.sort_values(["qvalue", "spearman_rho"], ascending=[True, False])


def plot_pseudotime_trends(adata, genes: list[str], outpath: Path) -> None:
    found = [gene for gene in genes if gene in adata.var_names]
    if not found or "pseudotime" not in adata.obs.columns:
        return
    ncols = 3
    nrows = int(np.ceil(len(found) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.2 * nrows), dpi=220)
    axes = np.atleast_1d(axes).ravel()
    pseudotime = adata.obs["pseudotime"].to_numpy(dtype=float)
    for axis, gene in zip(axes, found):
        frame = pd.DataFrame({"pseudotime": pseudotime, "expression": dense_vector(adata, gene)}).sort_values("pseudotime")
        sns.scatterplot(data=frame, x="pseudotime", y="expression", s=5, alpha=0.10, linewidth=0, ax=axis)
        sns.regplot(
            data=frame, x="pseudotime", y="expression", lowess=True, scatter=False,
            color="#C62828", line_kws={"linewidth": 2}, ax=axis,
        )
        axis.set_title(gene, fontweight="bold")
        axis.set_ylabel("Log-normalized expression")
    for axis in axes[len(found):]:
        axis.axis("off")
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close("all")


def plot_score_correlations(adata, columns: list[str], outpath: Path) -> None:
    columns = [column for column in columns if column in adata.obs.columns]
    if len(columns) < 2:
        return
    correlation = adata.obs[columns].corr(method="spearman")
    fig, axis = plt.subplots(figsize=(7, 5.5), dpi=300)
    sns.heatmap(correlation, annot=True, fmt=".2f", cmap="vlag", vmin=-1, vmax=1, ax=axis)
    axis.set_title("Spearman correlations among untreated-cell scores", fontweight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close("all")


def batch_gene_summary(adata, genes: list[str], groupby: list[str]) -> pd.DataFrame:
    records = []
    for gene in genes:
        if gene not in adata.var_names:
            continue
        frame = adata.obs[groupby].astype(str).copy()
        frame["expression"] = dense_vector(adata, gene)
        aggregate = frame.groupby(groupby, observed=True)["expression"].agg(
            n_cells="size", mean_expression="mean", median_expression="median"
        )
        detected = (
            frame.assign(detected=frame["expression"] > 0)
            .groupby(groupby, observed=True)["detected"]
            .mean()
            .mul(100)
            .rename("pct_detected")
        )
        summary = pd.concat([aggregate, detected], axis=1).reset_index()
        summary.insert(0, "gene", gene)
        records.append(summary)
    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()


def main() -> None:
    start_time = time.time()
    write_environment_report(AUDIT_DIR / "software_environment.txt")

    all_path = RUN_DIR / "adata_all_cells_hvg_fixed.h5ad"
    ut_path = RUN_DIR / "adata_ut_full_genes_fixed.h5ad"
    if not all_path.exists():
        raise FileNotFoundError(f"Missing upstream all-cell object: {all_path}")
    if not ut_path.exists():
        raise FileNotFoundError(f"Missing upstream UT object: {ut_path}")

    print("Loading finalized upstream AnnData objects", flush=True)
    adata_all = sc.read_h5ad(all_path)
    adata_ut = sc.read_h5ad(ut_path)

    required_all = {"Treatment", "Batch", "cellular_states", "leiden"}
    absent_all = sorted(required_all.difference(adata_all.obs.columns))
    if absent_all:
        raise ValueError(f"Missing metadata in all-cell object: {absent_all}")
    if "cellular_states" not in adata_ut.obs.columns:
        raise ValueError("UT object is missing cellular_states")

    validation = load_json(RUN_DIR / "input_validation.json")
    upstream_audit = load_json(RUN_DIR / "run_audit.json")
    score_columns = [
        column for column in ["S_score", "G2M_score", "Inflammatory_Score", "EMT_score", "Epithelial_score"]
        if column in adata_ut.obs.columns
    ]
    downstream_audit = {
        "workflow_version": "v3",
        "run_dir": str(RUN_DIR),
        "all_cell_shape": [int(adata_all.n_obs), int(adata_all.n_vars)],
        "ut_shape": [int(adata_ut.n_obs), int(adata_ut.n_vars)],
        "all_cell_has_umap": "X_umap" in adata_all.obsm,
        "ut_has_umap": "X_umap" in adata_ut.obsm,
        "has_pseudotime": "pseudotime" in adata_ut.obs.columns,
        "available_scores": score_columns,
        "input_validation": validation,
        "upstream_run_audit": upstream_audit,
        "scope_note": (
            "Real Treatment, Batch, and cellular-state metadata only. No artificial EBV grouping. "
            "Treatment contrasts, state-marker ranking, pseudotime and LIANA analyses are exploratory."
        ),
    }
    (AUDIT_DIR / "downstream_input_audit.json").write_text(json.dumps(downstream_audit, indent=2), encoding="utf-8")
    pd.DataFrame({"all_cell_obs_column": list(adata_all.obs.columns)}).to_csv(AUDIT_DIR / "all_cell_obs_columns.csv", index=False)
    pd.DataFrame({"ut_obs_column": list(adata_ut.obs.columns)}).to_csv(AUDIT_DIR / "ut_obs_columns.csv", index=False)

    print("Module 1: all-cell UMAP and composition", flush=True)
    sc.pl.umap(adata_all, color=["leiden", "Treatment", "Batch", "cellular_states"], ncols=2, wspace=0.45, show=False)
    plt.savefig(FIG_DIR / "Figure8A_all_cell_umap_metadata.pdf", bbox_inches="tight")
    plt.close("all")

    state_order = ordered_categories(adata_all.obs["cellular_states"])
    treatment_order = ordered_categories(adata_all.obs["Treatment"], ["UT", "T24", "T48"])
    counts = adata_all.obs.groupby(["Treatment", "Batch", "cellular_states"], observed=True).size().rename("n_cells").reset_index()
    counts.to_csv(TABLE_DIR / "treatment_batch_cellular_state_counts.csv", index=False)
    treatment_composition = pd.crosstab(adata_all.obs["Treatment"], adata_all.obs["cellular_states"], normalize="index").mul(100).reindex(index=treatment_order, columns=state_order, fill_value=0)
    treatment_composition.to_csv(TABLE_DIR / "treatment_by_cellular_state_percent.csv")
    plot_stacked_composition(treatment_composition, "Descriptive cellular-state composition across experimental conditions", FIG_DIR / "Figure8B_treatment_state_composition.pdf")

    batch_order = ordered_categories(adata_all.obs["Batch"])
    batch_composition = pd.crosstab(adata_all.obs["Batch"], adata_all.obs["cellular_states"], normalize="index").mul(100).reindex(index=batch_order, columns=state_order, fill_value=0)
    batch_composition.to_csv(TABLE_DIR / "batch_by_cellular_state_percent.csv")
    plot_stacked_composition(batch_composition, "Cellular-state composition by batch", FIG_DIR / "Supplementary_batch_state_composition.pdf")

    print("Module 2: state markers and hub genes", flush=True)
    state_marker_key = "state_markers_exploratory"
    sc.tl.rank_genes_groups(adata_all, groupby="cellular_states", method="wilcoxon", use_raw=False, key_added=state_marker_key, n_genes=min(300, adata_all.n_vars))
    state_markers = rank_groups_to_dataframe(adata_all, state_marker_key)
    state_markers["analysis_note"] = "Exploratory cell-state marker ranking; cells are not independent biological replicates."
    state_markers.to_csv(TABLE_DIR / "exploratory_state_markers_all_cells.csv", index=False)
    state_markers.groupby("group", observed=True).head(20).to_csv(TABLE_DIR / "exploratory_top20_markers_per_state.csv", index=False)
    plot_state_marker_heatmap(adata_all, state_marker_key, FIG_DIR / "Supplementary_state_marker_heatmap.pdf")

    hub_genes = ["CENPF", "TOP2A", "TPX2", "MKI67", "BIRC5", "CDK1", "UBE2C", "NUSAP1", "CCNB2"]
    plot_marker_dotplot(adata_ut, hub_genes, "cellular_states", "Mitotic-marker expression across untreated NPC43 cellular states", FIG_DIR / "Figure8C_mitotic_hub_genes_by_state.pdf")
    expression_summary(adata_ut, hub_genes, "cellular_states").to_csv(TABLE_DIR / "mitotic_hub_gene_expression_by_ut_state.csv", index=False)

    state_marker_panel = [
        "EPCAM", "KRT5", "KRT14", "KRT17", "KRT19", "KRT15", "KRT13",
        "VIM", "FN1", "THBS1", "ITGA3", "S100A10",
        "CXCL8", "IL1A", "CXCL1", "CXCL2", "CCL20", "TNFAIP3", "NDRG1",
        "MKI67", "TOP2A", "CENPF", "BIRC5",
    ]
    plot_marker_dotplot(adata_ut, state_marker_panel, "cellular_states", "Selected epithelial, ECM, inflammatory, and cell-cycle markers in untreated NPC43 states", FIG_DIR / "Supplementary_ut_state_marker_dotplot.pdf")

    print("Module 3: bulk/proteomic candidate mapping", flush=True)
    gene_groups = {
        "Low_CP_myeloid_inflammatory": ["CXCL8", "CXCL2", "IL1A", "IL1B", "IL6", "TREM1", "CSF3R", "S100A8", "S100A9", "S100A12", "SLC11A1", "IL17A"],
        "High_CP_lymphoid_chemokine": ["CCL19", "CCL21", "CXCL9", "CD3D", "CD247", "MS4A1", "SPIB", "FASLG", "VCAM1", "C7", "BCL2"],
        "TPA_proteomic_context": ["FN1", "S100P", "EPHA2", "ITGA3", "S100A10", "THBS1", "STMN1", "STAT3", "STAT6", "STAT1", "PFN2"],
    }
    summaries = []
    for set_name, genes in gene_groups.items():
        state_summary = expression_summary(adata_ut, genes, "cellular_states")
        state_summary.insert(0, "gene_set", set_name)
        summaries.append(state_summary)
        state_summary.to_csv(TABLE_DIR / f"{set_name}_gene_detection_by_ut_state.csv", index=False)
        treatment_summary = expression_summary(adata_all, genes, "Treatment")
        treatment_summary.insert(0, "gene_set", set_name)
        treatment_summary.to_csv(TABLE_DIR / f"{set_name}_gene_detection_by_treatment.csv", index=False)
        plot_marker_dotplot(adata_ut, genes, "cellular_states", f"{set_name.replace('_', ' ')} genes across untreated NPC43 states", FIG_DIR / f"Supplementary_{set_name}_by_ut_state.pdf")
    pd.concat(summaries, ignore_index=True).to_csv(TABLE_DIR / "clinical_bulk_and_proteomic_gene_detection_by_ut_state.csv", index=False)
    combined_genes = [gene for genes in gene_groups.values() for gene in genes]
    plot_marker_dotplot(adata_ut, combined_genes, "cellular_states", "Clinical bulk and proteomic candidate genes in untreated NPC43 states", FIG_DIR / "Figure8D_clinical_and_proteomic_gene_mapping.pdf")

    print("Module 4: UT scores and pseudotime", flush=True)
    if score_columns:
        sc.pl.umap(adata_ut, color=score_columns, ncols=3, show=False)
        plt.savefig(FIG_DIR / "Figure8F_ut_full_gene_scores_umap.pdf", bbox_inches="tight")
        plt.close("all")
        sc.pl.violin(adata_ut, keys=score_columns, groupby="cellular_states", stripplot=False, rotation=45, show=False)
        plt.savefig(FIG_DIR / "Supplementary_ut_scores_by_state.pdf", bbox_inches="tight")
        plt.close("all")
        plot_score_correlations(adata_ut, score_columns, FIG_DIR / "Supplementary_ut_score_spearman_correlation.pdf")

    if "pseudotime" in adata_ut.obs.columns:
        sc.pl.umap(adata_ut, color=["pseudotime", "cellular_states"], ncols=2, show=False)
        plt.savefig(FIG_DIR / "Figure8E_ut_pseudotime_and_states.pdf", bbox_inches="tight")
        plt.close("all")
        adata_ut.obs.groupby("cellular_states", observed=True)["pseudotime"].agg(["size", "mean", "median", "min", "max"]).reset_index().to_csv(TABLE_DIR / "ut_pseudotime_by_state.csv", index=False)

    if not ARGS.skip_enrichment:
        print("Module 5: enrichment", flush=True)
        libraries = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
        selected_states = [state for state in ["UT Cycling-G2M", "UT Stressed cells", "UT Differentiated cells", "UT ECM-related (1) cells", "UT ECM-related (2) cells"] if state in set(state_markers["group"])]
        enriched_states = []
        for state in selected_states:
            genes = state_markers.loc[(state_markers["group"] == state) & (state_markers["pvals_adj"] < 0.05) & (state_markers["logfoldchanges"] > 0.25), "gene"].head(150).tolist()
            safe_name = state.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", "")
            result = enrich_gene_list(genes, f"state_markers_{safe_name}", ENRICH_DIR / f"state_markers_{safe_name}", libraries)
            if result is not None and not result.empty:
                result = result.copy()
                result.insert(0, "cellular_state", state)
                enriched_states.append(result)
                plot_enrichment_table(result, f"Functional annotation of markers: {state}", FIG_DIR / f"Supplementary_GO_KEGG_{safe_name}.pdf", "viridis")
        if enriched_states:
            pd.concat(enriched_states, ignore_index=True).to_csv(TABLE_DIR / "state_marker_enrichment_fdr05.csv", index=False)

        if "pseudotime" in adata_ut.obs.columns:
            print("Computing gene-pseudotime correlations", flush=True)
            pt_result = pseudotime_correlations(adata_ut, min_pct_detected=5.0)
            pt_result.to_csv(TABLE_DIR / "ut_gene_pseudotime_spearman_correlations.csv", index=False)
            increasing = pt_result[(pt_result["qvalue"] < 0.05) & (pt_result["spearman_rho"] >= 0.15)]
            decreasing = pt_result[(pt_result["qvalue"] < 0.05) & (pt_result["spearman_rho"] <= -0.15)]
            increasing.to_csv(TABLE_DIR / "ut_pseudotime_increasing_genes.csv", index=False)
            decreasing.to_csv(TABLE_DIR / "ut_pseudotime_decreasing_genes.csv", index=False)
            up = enrich_gene_list(increasing.head(300)["gene"].tolist(), "pseudotime_increasing", ENRICH_DIR / "pseudotime_increasing", libraries)
            down = enrich_gene_list(decreasing.head(300)["gene"].tolist(), "pseudotime_decreasing", ENRICH_DIR / "pseudotime_decreasing", libraries)
            plot_enrichment_table(up, "Functional annotation of genes increasing along inferred pseudotime", FIG_DIR / "Figure8G_pseudotime_increasing_GO_KEGG.pdf", "Reds_r")
            plot_enrichment_table(down, "Functional annotation of genes decreasing along inferred pseudotime", FIG_DIR / "Figure8G_pseudotime_decreasing_GO_KEGG.pdf", "Blues_r")
            trend_genes = increasing.head(6)["gene"].tolist() + decreasing.head(6)["gene"].tolist()
            plot_pseudotime_trends(adata_ut, trend_genes, FIG_DIR / "Supplementary_top_pseudotime_gene_trends.pdf")

    if ARGS.with_treatment_exploratory:
        print("Module 6: exploratory UT versus T48 ranking", flush=True)
        treatment_data = adata_all[adata_all.obs["Treatment"].isin(["UT", "T48"])].copy()
        treatment_data.obs["Treatment"] = treatment_data.obs["Treatment"].astype("category")
        treatment_data.obs["Treatment"] = treatment_data.obs["Treatment"].cat.reorder_categories(["UT", "T48"])
        treatment_key = "exploratory_T48_vs_UT_cell_level"
        sc.tl.rank_genes_groups(treatment_data, groupby="Treatment", reference="UT", method="wilcoxon", use_raw=False, key_added=treatment_key, n_genes=min(300, treatment_data.n_vars))
        treatment_rank = rank_groups_to_dataframe(treatment_data, treatment_key)
        treatment_rank["analysis_note"] = "Exploratory cell-level ranking only; cells are not independent biological replicates. Do not interpret p values as confirmatory treatment-level inference."
        treatment_rank.to_csv(TABLE_DIR / "exploratory_cell_level_T48_vs_UT_rankings.csv", index=False)
        t48 = treatment_rank[treatment_rank["group"] == "T48"].copy()
        t48_up = t48[(t48["pvals_adj"] < 0.05) & (t48["logfoldchanges"] > 0.25)]
        t48_down = t48[(t48["pvals_adj"] < 0.05) & (t48["logfoldchanges"] < -0.25)]
        t48_up.to_csv(TABLE_DIR / "exploratory_T48_associated_positive_rankings.csv", index=False)
        t48_down.to_csv(TABLE_DIR / "exploratory_T48_associated_negative_rankings.csv", index=False)
        batch_gene_summary(treatment_data, ["CXCL8", "IL1A", "CXCL1", "CXCL2", "CCL20", "TNFAIP3", "IL1RN", "IL1B", "FN1", "THBS1", "ITGA3", "S100A10", "MKI67", "TOP2A", "CENPF", "BIRC5"], ["Treatment", "Batch", "cellular_states"]).to_csv(TABLE_DIR / "exploratory_key_gene_batch_treatment_state_summary.csv", index=False)

    if ARGS.with_liana:
        print("Module 7: exploratory LIANA", flush=True)
        if li is None:
            print("[warning] liana is not installed; skipping LIANA.", flush=True)
        else:
            try:
                liana_data = adata_ut.copy()
                state_column = "native_cellular_states" if "native_cellular_states" in liana_data.obs.columns else "cellular_states"
                state_counts = liana_data.obs[state_column].value_counts()
                retained = state_counts[state_counts >= 30].index.tolist()
                liana_data = liana_data[liana_data.obs[state_column].isin(retained)].copy()
                if len(retained) < 2:
                    raise ValueError("Fewer than two states contain at least 30 cells")
                li.mt.rank_aggregate(liana_data, groupby=state_column, use_raw=False, verbose=True)
                liana_result = liana_data.uns["liana_res"].copy()
                liana_result["analysis_note"] = "Exploratory ligand-receptor co-expression prediction among NPC43 epithelial states; not direct evidence of functional communication or immune-cell crosstalk."
                liana_result.to_csv(TABLE_DIR / "exploratory_liana_epithelial_state_results.csv", index=False)
                if "ligand_complex" in liana_result.columns:
                    candidates = ["CXCL8", "IL1A", "CXCL1", "CXCL2", "CCL20", "TNF", "TGFB1", "FN1", "THBS1"]
                    hits = liana_result[liana_result["ligand_complex"].isin(candidates)].copy()
                    sort_col = "magnitude_rank" if "magnitude_rank" in hits.columns else hits.columns[0]
                    hits.sort_values(sort_col).head(50).to_csv(TABLE_DIR / "exploratory_liana_candidate_ligand_hits.csv", index=False)
            except Exception as error:
                (AUDIT_DIR / "liana_error.txt").write_text(f"{type(error).__name__}: {error}\n", encoding="utf-8")
                print(f"[warning] LIANA failed: {type(error).__name__}: {error}", flush=True)

    manifest = []
    for path in sorted(OUT_DIR.rglob("*")):
        if path.is_file():
            manifest.append({"relative_path": str(path.relative_to(OUT_DIR)), "bytes": path.stat().st_size})
    pd.DataFrame(manifest).to_csv(AUDIT_DIR / "downstream_output_manifest.csv", index=False)
    elapsed_minutes = round((time.time() - start_time) / 60, 2)
    (AUDIT_DIR / "README_scope.txt").write_text(
        "\n".join([
            "NPC43 clean downstream workflow scope (v3)",
            "",
            "- Reads finalized outputs generated by run_single_cell_jupyterhub_fixed.py.",
            "- Does not rerun upstream preprocessing or trajectory inference.",
            "- Uses real Treatment, Batch, and cellular-state metadata only.",
            "- Does not create EBV-high/EBV-low labels from cellular states.",
            "- State-marker ranking, pseudotime, treatment ranking, enrichment, and LIANA are exploratory.",
            "- Cell-level rankings are not confirmatory treatment-level biological inference.",
            "- LIANA, if run, represents candidate epithelial-state ligand-receptor co-expression only.",
            "- NPC43 provides epithelial cell-intrinsic context and does not directly validate plasma EBV DNA groups in tumors.",
            "",
            f"Runtime minutes: {elapsed_minutes}",
        ]) + "\n",
        encoding="utf-8",
    )
    print(f"Completed clean downstream analysis v3: {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
