#!/usr/bin/env python3
"""
00 — Extract atlas obs metadata to parquet (supplementary-figure cache)
======================================================================

Purpose
    Read the atlas h5ads and cache lightweight obs parquets so the
    supplementary QC figure (SF1) can render without loading the full atlas.
    Run once; SF1 reads from these cached parquets.

INPUTS  (via config)
    - atlas concatenated pre-QC counts:
        data_root/2026_final_atlas/pre-processing and integration/processed/
        atlas_concat_counts_only_X.h5ad   (2.73M cells, pre-QC; QC metrics computed here)
    - integration/trust-boundary object: obj("atlas_scanvi")  (final filtered, post doublet removal)

OUTPUTS
    - output_root/figures/data/atlas_obs_prefilter.parquet
    - output_root/figures/data/atlas_obs_postfilter.parquet

MANUSCRIPT PANEL(S)
    Upstream cache for SF1A-C (QC metrics by study).

RUNTIME TIER
    moderate (loads the pre-QC concat matrix to compute total_counts / n_genes).
"""

import gc
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

RAW_H5AD = path(
    "data_root", "2026_final_atlas", "pre-processing and integration",
    "processed", "atlas_concat_counts_only_X.h5ad",
)
POST_H5AD = obj("atlas_scanvi")

CACHE_PRE = path("output_root", "figures", "data", "atlas_obs_prefilter.parquet")
CACHE_POST = path("output_root", "figures", "data", "atlas_obs_postfilter.parquet")

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
    total_counts = np.array(X.sum(axis=1)).flatten()
    if sp.issparse(X):
        n_genes = np.array((X > 0).sum(axis=1)).flatten()
    else:
        n_genes = np.array((X > 0).sum(axis=1)).flatten()

    df = adata.obs[["study", "sample_id", "patient_id"]].copy()
    df["total_counts"] = total_counts
    df["n_genes_by_counts"] = n_genes

    del adata, X, total_counts, n_genes
    gc.collect()

    df.to_parquet(cache_path)
    print(f"    Saved: {cache_path} ({len(df):,} cells)")


print("=" * 60)
print("00 — Extract atlas obs to parquet")
print("=" * 60)

extract_prefilter(RAW_H5AD, CACHE_PRE)
extract_postfilter(POST_H5AD, CACHE_POST)

print("\nDone. SF1 can now read from output_root/figures/data/*.parquet")
