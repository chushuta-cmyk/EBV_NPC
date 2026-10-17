#!/usr/bin/env python3
"""Run the corrected NPC43 single-cell workflow on JupyterHub or locally."""

from __future__ import annotations

import argparse
import json
import os
import resource
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260609)
    parser.add_argument("--num-threads", type=int, default=1)
    parser.add_argument("--skip-palantir", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    return parser.parse_args()


def find_input(folder: Path, stem: str) -> Path:
    # Accept compressed and uncompressed GEO files.
    for name in (stem, f"{stem}.gz"):
        path = folder / name
        if path.exists():
            return path
    raise FileNotFoundError(f"Missing {stem} or {stem}.gz in {folder}")


def peak_memory_gb() -> float:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes. Linux reports KiB.
    if value > 10_000_000:
        return value / 1024**3
    return value / 1024**2


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = args.output_dir / "figures"
    figure_dir.mkdir(exist_ok=True)

    # Put caches in the output folder for JupyterHub and restricted systems.
    os.environ.setdefault("NUMBA_CACHE_DIR", str(args.output_dir / ".numba_cache"))
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".mpl_cache"))
    os.environ.setdefault("NUMBA_NUM_THREADS", str(args.num_threads))
    os.environ.setdefault("OMP_NUM_THREADS", str(args.num_threads))
    os.environ.setdefault("OPENBLAS_NUM_THREADS", str(args.num_threads))

    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc

    np.random.seed(args.seed)
    sc.settings.verbosity = 2
    sc.settings.figdir = str(figure_dir)
    sc.set_figure_params(dpi=180, facecolor="white", frameon=True)
    plt.rcParams["pdf.fonttype"] = 42

    matrix_path = find_input(args.data_dir, "npc43_gene_exp.mtx")
    gene_path = find_input(args.data_dir, "npc43_genes.tsv")
    barcode_path = find_input(args.data_dir, "npc43_barcodes.tsv")
    annotation_path = find_input(args.data_dir, "npc43_cell_annotation.csv")

    print("Loading raw sparse matrix")
    matrix = sc.read_mtx(matrix_path).X
    genes = pd.read_csv(gene_path, sep="\t", header=None, compression="infer")
    barcodes = pd.read_csv(barcode_path, sep="\t", header=None, compression="infer")
    annotation = pd.read_csv(annotation_path, index_col=0, compression="infer")

    gene_names = genes.iloc[:, 1] if genes.shape[1] > 1 else genes.iloc[:, 0]
    gene_names = gene_names.astype(str)
    barcode_names = barcodes.iloc[:, 0].astype(str)

    if matrix.shape == (len(gene_names), len(barcode_names)):
        matrix = matrix.T
    elif matrix.shape != (len(barcode_names), len(gene_names)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match "
            f"{len(gene_names)} genes and {len(barcode_names)} barcodes"
        )

    adata = ad.AnnData(X=matrix)
    adata.var_names = gene_names.to_numpy()
    adata.obs_names = barcode_names.to_numpy()
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # Match annotations by barcode, not row position.
    annotation.index = annotation.index.astype(str)
    if not annotation.index.is_unique:
        raise ValueError("Annotation barcodes are not unique")
    missing = adata.obs_names.difference(annotation.index)
    if len(missing):
        raise ValueError(f"{len(missing)} expression barcodes lack annotation")
    adata.obs = annotation.reindex(adata.obs_names).copy()

    required = {"Treatment", "cellular_states", "percent.mt"}
    absent = sorted(required.difference(adata.obs.columns))
    if absent:
        raise ValueError(f"Missing required annotation columns: {absent}")

    validation = {
        "matrix_shape": [int(adata.n_obs), int(adata.n_vars)],
        "matched_barcodes": int(adata.n_obs),
        "treatment_counts_raw": {
            str(k): int(v) for k, v in adata.obs["Treatment"].value_counts().items()
        },
        "peak_memory_gb_after_load": round(peak_memory_gb(), 3),
    }
    (args.output_dir / "input_validation.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )
    print(json.dumps(validation, indent=2))
    if args.validate_only:
        return

    print("Filtering and normalizing")
    sc.pp.filter_cells(adata, min_genes=1000)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["percent.mt"] < 20].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Use highly variable genes only for PCA, neighbors, clustering, and UMAP.
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack", random_state=args.seed)
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=15, random_state=args.seed)
    sc.tl.leiden(adata_hvg, resolution=0.5, random_state=args.seed)
    sc.tl.umap(adata_hvg, random_state=args.seed)

    # Keep coordinates and clusters with the full-gene object.
    adata.obs["leiden"] = adata_hvg.obs["leiden"].copy()
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"].copy()
    adata.uns["umap"] = adata_hvg.uns["umap"].copy()
    sc.pl.umap(
        adata,
        color=["leiden", "Treatment", "cellular_states"],
        wspace=0.4,
        show=False,
    )
    plt.savefig(figure_dir / "umap_all_cells_fixed.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    marker_genes = [
        "EPCAM", "KRT5", "KRT14", "KRT17", "KRT19",
        "VIM", "FN1", "SNAI1", "TWIST1", "ZEB1",
        "CXCL8", "IL1A", "TREM1", "CXCR2", "CSF3R",
        "MKI67", "PCNA", "CDK1",
    ]
    marker_genes = [gene for gene in marker_genes if gene in adata.var_names]
    sc.pl.dotplot(adata, marker_genes, groupby="leiden", show=False)
    plt.savefig(figure_dir / "marker_genes_all_genes.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    # The compact all-cell file contains the graph features only.
    adata_hvg.write_h5ad(
        args.output_dir / "adata_all_cells_hvg_fixed.h5ad", compression="gzip"
    )

    print("Preparing the untreated full-gene object")
    adata_ut = adata[adata.obs["Treatment"].astype(str) == "UT"].copy()
    sc.pp.highly_variable_genes(adata_ut, n_top_genes=2000, flavor="seurat")
    sc.tl.pca(
        adata_ut,
        n_comps=30,
        mask_var="highly_variable",
        random_state=args.seed,
    )
    sc.pp.neighbors(adata_ut, n_neighbors=30, n_pcs=30, random_state=args.seed)
    sc.tl.umap(adata_ut, random_state=args.seed)

    if not args.skip_palantir:
        import palantir

        print("Running Palantir")
        palantir.utils.run_diffusion_maps(adata_ut)
        multiscale = palantir.utils.determine_multiscale_space(adata_ut)
        root_mask = adata_ut.obs["cellular_states"].astype(str).eq(
            "UT Differentiated cells"
        )
        if not root_mask.any():
            raise ValueError("No UT Differentiated cells are available for the root")
        root_cell = adata_ut.obs_names[root_mask][0]
        result = palantir.core.run_palantir(
            multiscale, early_cell=root_cell, seed=args.seed
        )
        adata_ut.obs["pseudotime"] = result.pseudotime.reindex(adata_ut.obs_names)

    def score(genes_to_use: list[str], name: str) -> list[str]:
        found = [gene for gene in genes_to_use if gene in adata_ut.var_names]
        if len(found) < 2:
            print(f"Skipping {name}: only {len(found)} genes found")
            return found
        sc.tl.score_genes(
            adata_ut,
            gene_list=found,
            score_name=name,
            random_state=args.seed,
        )
        return found

    s_genes = [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
        "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2",
        "RPA2", "NASP", "RAD51AP1", "GMNN",
    ]
    g2m_genes = [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
        "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
        "TACC3", "SMC4", "CCNB2",
    ]
    s_found = [gene for gene in s_genes if gene in adata_ut.var_names]
    g2m_found = [gene for gene in g2m_genes if gene in adata_ut.var_names]
    if len(s_found) >= 2 and len(g2m_found) >= 2:
        sc.tl.score_genes_cell_cycle(
            adata_ut, s_genes=s_found, g2m_genes=g2m_found
        )

    coverage = {
        "S_score": [len(s_found), len(s_genes)],
        "G2M_score": [len(g2m_found), len(g2m_genes)],
        "Inflammatory_Score": [
            len(score(["CXCL8", "IL1A", "TREM1", "CXCR2", "CSF3R"], "Inflammatory_Score")),
            5,
        ],
        "EMT_score": [
            len(score(["VIM", "FN1", "SNAI1", "TWIST1", "ZEB1"], "EMT_score")),
            5,
        ],
        "Epithelial_score": [
            len(score(["EPCAM", "KRT5", "KRT14", "KRT17", "KRT19"], "Epithelial_score")),
            5,
        ],
    }

    plot_columns = [
        column
        for column in (
            "S_score", "G2M_score", "Inflammatory_Score", "EMT_score",
            "Epithelial_score", "pseudotime",
        )
        if column in adata_ut.obs
    ]
    sc.pl.umap(adata_ut, color=plot_columns, show=False)
    plt.savefig(figure_dir / "ut_scores_fixed.png", dpi=300, bbox_inches="tight")
    plt.close("all")

    adata_ut.write_h5ad(
        args.output_dir / "adata_ut_full_genes_fixed.h5ad", compression="gzip"
    )

    audit = {
        "cells_after_qc": int(adata.n_obs),
        "genes_after_qc": int(adata.n_vars),
        "hvg_count": int(adata_hvg.n_vars),
        "untreated_cells": int(adata_ut.n_obs),
        "gene_set_coverage": coverage,
        "palantir_completed": bool("pseudotime" in adata_ut.obs),
        "peak_memory_gb": round(peak_memory_gb(), 3),
        "note": (
            "Treatment differential expression is not tested here because T24 has one batch "
            "and cell-level tests would treat cells as independent biological replicates."
        ),
    }
    (args.output_dir / "run_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
