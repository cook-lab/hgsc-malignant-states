#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/03b_refilter_umap.py, the script that ran.
# Kept: the Scrublet < 0.25 re-filter, neighbours and UMAP that write hgsc_atlas_scanvi.h5ad (read by level-1 annotation).
# Removed: the QC UMAP plots for approval, their style settings and imports, and the closing summary print.
# The untrimmed script is in the lab project folder.
"""
Step 03b — Re-filter doublets (threshold 0.25) and recompute UMAP

Replicates scripts/03_filter_umap.ipynb with a stricter doublet score
threshold (0.25 vs 0.275).

Pipeline:
  1. Load integrated_scanvi.h5ad (2,326,532 cells)
  2. Remove doublet_scrublet == True cells
  3. Remove doublet_score_scrublet >= 0.25 (was 0.275)
  4. Recompute neighbors (n=10, X_scanvi)
  5. Recompute UMAP (min_dist=0.2)
  6. Save as hgsc_atlas_scanvi.h5ad (overwrites previous)
  7. Quick UMAP plot for visual check

Style: Cook Lab v1.2
"""

import os
import warnings

import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

PROJECT_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas"
INPUT_H5AD  = os.path.join(
    PROJECT_DIR,
    "pre-processing and integration",
    "20260213_integration",
    "output", "anndata", "integrated_scanvi.h5ad",
)
OUTPUT_H5AD = os.path.join(PROJECT_DIR, "hgsc_atlas_scanvi.h5ad")

SCORE_THRESHOLD = 0.25  # stricter than previous 0.275

# ============================================================================
# LOAD
# ============================================================================

print("=" * 60)
print("Step 03b — Re-filter doublets (threshold 0.25) + UMAP")
print("=" * 60)

print(f"\n1. Loading: {INPUT_H5AD}", flush=True)
adata = sc.read_h5ad(INPUT_H5AD)
print(f"   Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# ============================================================================
# FILTER 1 — Remove doublet_scrublet == True
# ============================================================================

print("\n2. Removing doublet_scrublet == True...", flush=True)
n_before = adata.n_obs

doublet_col = adata.obs["doublet_scrublet"]
if doublet_col.dtype == object:
    mask_remove = doublet_col.astype(str).str.lower() == "true"
else:
    mask_remove = doublet_col == True

n_doublets = int(mask_remove.sum())
mask_keep = ~mask_remove
adata._inplace_subset_obs(mask_keep)

print(f"   Before:  {n_before:,}")
print(f"   Removed: {n_doublets:,} (doublet_scrublet == True)")
print(f"   After:   {adata.n_obs:,}")

# ============================================================================
# FILTER 2 — Score threshold 0.25
# ============================================================================

print(f"\n3. Removing doublet_score_scrublet >= {SCORE_THRESHOLD}...", flush=True)
n_before2 = adata.n_obs

mask_keep2 = adata.obs["doublet_score_scrublet"].fillna(0) < SCORE_THRESHOLD
n_removed2 = int((~mask_keep2).sum())

adata._inplace_subset_obs(mask_keep2)

print(f"   Before:  {n_before2:,}")
print(f"   Removed: {n_removed2:,} (score >= {SCORE_THRESHOLD})")
print(f"   After:   {adata.n_obs:,}")

total_removed = n_doublets + n_removed2
print(f"\n   Total removed: {total_removed:,} ({total_removed / n_before * 100:.2f}%)")
print(f"   Final shape:   {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# ============================================================================
# NEIGHBORS + UMAP
# ============================================================================

print("\n4. Computing neighbors (n=10, X_scanvi)...", flush=True)
sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=10)

print("\n5. Computing UMAP (min_dist=0.2)...", flush=True)
sc.tl.umap(adata, min_dist=0.2)

# ============================================================================
# SAVE
# ============================================================================

print(f"\n6. Saving: {OUTPUT_H5AD}", flush=True)
adata.write_h5ad(OUTPUT_H5AD)
fsize = os.path.getsize(OUTPUT_H5AD) / (1024**3)
print(f"   File size: {fsize:.1f} GB")
