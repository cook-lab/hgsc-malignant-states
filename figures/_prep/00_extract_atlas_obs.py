#!/usr/bin/env python3
"""
00_extract_atlas_obs — cache atlas obs metadata to parquet for fast figure rendering
====================================================================================
PURPOSE
    Read the atlas h5ad files, extract obs columns (computing QC metrics where
    needed), and write lightweight parquet caches that the QC supplemental
    figure scripts read instead of the multi-GB h5ad. Run once.

INPUTS
    - obj("atlas_scanvi")                       post-filter final atlas (obs)
    - <data_root>/2026_final_atlas/pre-processing and integration/processed/
        atlas_concat_counts_only_X.h5ad         pre-filter raw concat (X for QC)

OUTPUTS (under fig_data_dir/_prep_caches/)
    - atlas_obs_prefilter.parquet   (2.73M cells, pre-QC, QC metrics computed)
    - atlas_obs_postfilter.parquet  (final filtered atlas obs)

MANUSCRIPT PANEL(S): SF1A-C (QC violins by study) — consumed by
    supplementary/atlas_qc_metrics_by_study.py.

RUNTIME TIER: heavy (loads the full pre-filter matrix to compute QC metrics).
"""
from __future__ import annotations

import gc
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# --- Central config (this file is 2 levels under the repo root) -------------
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path  # noqa: E402

# ============================================================================
# PATHS
# ============================================================================

RAW_H5AD = obj("atlas_scanvi").replace(
    "hgsc_atlas_scanvi.h5ad",
    os.path.join("pre-processing and integration", "processed",
                 "atlas_concat_counts_only_X.h5ad"),
)
POST_H5AD = obj("atlas_scanvi")

CACHE_PRE = path("output_root", "_prep_caches", "atlas_obs_prefilter.parquet")
CACHE_POST = path("output_root", "_prep_caches", "atlas_obs_postfilter.parquet")

# ============================================================================
# COLUMNS
# ============================================================================

POST_OBS_COLS = [
    "study", "sample_id", "patient_id",
    "total_counts", "n_genes_by_counts", "doublet_score_scrublet",
]
NUMERIC_COLS = ["total_counts", "n_genes_by_counts", "doublet_score_scrublet"]


def extract_postfilter(h5ad_path, cache_path):
    """Load final atlas backed, extract obs, cache as parquet."""
    if os.path.exists(cache_path):
        print(f"  SKIP post-filter — cache exists: {cache_path}")
        return

    import scanpy as sc

    print(f"  Loading post-filter (backed): {h5ad_path}", flush=True)
    adata = sc.read_h5ad(h5ad_path, backed="r")
    print(f"    Shape: {adata.shape[0]:,} x {adata.shape[1]:,}")

    cols_present = [c for c in POST_OBS_COLS if c in adata.obs.columns]
    df = adata.obs[cols_present].copy()
    for col in NUMERIC_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    adata.file.close()
    del adata
    gc.collect()

    df.to_parquet(cache_path)
    print(f"    Saved: {cache_path} ({len(df):,} cells)")


def extract_prefilter(h5ad_path, cache_path):
    """Load raw concatenated atlas, compute QC metrics, cache as parquet."""
    if os.path.exists(cache_path):
        print(f"  SKIP pre-filter — cache exists: {cache_path}")
        return

    import scanpy as sc
    import scipy.sparse as sp

    print(f"  Loading pre-filter: {h5ad_path}", flush=True)
    adata = sc.read_h5ad(h5ad_path)
    print(f"    Shape: {adata.shape[0]:,} x {adata.shape[1]:,}")

    print("    Computing total_counts and n_genes_by_counts...", flush=True)
    X = adata.X
    if sp.issparse(X):
        total_counts = np.array(X.sum(axis=1)).flatten()
        n_genes = np.array((X > 0).sum(axis=1)).flatten()
    else:
        total_counts = np.array(X.sum(axis=1)).flatten()
        n_genes = np.array((X > 0).sum(axis=1)).flatten()

    df = adata.obs[["study", "sample_id", "patient_id"]].copy()
    df["total_counts"] = total_counts
    df["n_genes_by_counts"] = n_genes

    del adata, X, total_counts, n_genes
    gc.collect()

    df.to_parquet(cache_path)
    print(f"    Saved: {cache_path} ({len(df):,} cells)")


if __name__ == "__main__":
    print("=" * 60)
    print("00 — Extract atlas obs to parquet")
    print("=" * 60)
    extract_prefilter(RAW_H5AD, CACHE_PRE)
    extract_postfilter(POST_H5AD, CACHE_POST)
    print("\nDone. Figure scripts can now read from the _prep_caches parquets.")
