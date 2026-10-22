#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/09d_publication_figures.py, the script that ran.
# Kept: step 1, the final filter that drops the level-2 "Excluded" clusters and writes hgsc_atlas_final.h5ad.
# Removed: the level-1, study and level-2 UMAP figures with their style, palettes and plotting helpers.
# The untrimmed script is in the lab project folder.
"""
Step 09d — Publication-ready UMAPs from filtered atlas
========================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Creates a QC-filtered h5ad (Excluded cells removed) and renders
three publication-ready UMAP figures:

  A — celltype_level1 (Cook Lab palette v1.2)
  B — study           (Kelly 22 palette)
  C — celltype_level2 (shades of level1 hues)

Each figure is saved three ways:
  *.pdf  — rasterized dots at 600 dpi (all cells), vector text/legend
  *.svg  — same as PDF but in SVG format
  *_vector.svg — pure vector (150k subsampled), every dot editable
                 in Adobe Illustrator

Sizing: individual panels fit on half a US letter page (8.5 × 5.5 in).

Usage:
    python scripts/09d_publication_figures.py

Output:
    hgsc_atlas_celltype_level2_filtered.h5ad  (project root)
    output/09d_publication_figures/*.pdf, *.svg
"""

import os
import time
import gc
import numpy as np
import pandas as pd
import anndata as ad

# ── paths ────────────────────────────────────────────────────
BASE = "/Volumes/CookLab/Sarah/2026_final_atlas"
INPUT_H5AD    = os.path.join(BASE, "hgsc_atlas_celltype_level2.h5ad")
FILTERED_H5AD = os.path.join(BASE, "hgsc_atlas_final.h5ad")

# ============================================================================
# MAIN
# ============================================================================

print("=" * 70)
print("  Step 09d — Publication-ready figures (filtered atlas)")
print(f"  {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 70)

# ──────────────────────────────────────────────────────────────
# STEP 1: Create filtered h5ad (remove Excluded cells)
# ──────────────────────────────────────────────────────────────

print("\n1. Creating filtered h5ad...")

if os.path.exists(FILTERED_H5AD):
    print(f"   Already exists: {FILTERED_H5AD}")
    print(f"   Loading filtered atlas (backed)...")
    t0 = time.time()
    adata = ad.read_h5ad(FILTERED_H5AD, backed="r")
    umap   = np.array(adata.obsm["X_umap"])
    level1 = np.array(adata.obs["celltype_level1"])
    level2 = np.array(adata.obs["celltype_level2"])
    study  = np.array(adata.obs["study"])
    n_cells = adata.n_obs
    print(f"   {n_cells:,} cells loaded in {time.time()-t0:.0f}s")
    del adata
    gc.collect()
else:
    print(f"   Loading full atlas: {INPUT_H5AD}")
    t0 = time.time()
    adata = ad.read_h5ad(INPUT_H5AD)
    n_before = adata.n_obs
    print(f"   Full atlas: {n_before:,} cells (loaded in {time.time()-t0:.0f}s)")

    # Filter — excluded labels start with "Excluded" (e.g. "Excluded B cell")
    mask = ~adata.obs["celltype_level2"].astype(str).str.startswith("Excluded")
    adata = adata[mask].copy()
    n_after = adata.n_obs
    n_removed = n_before - n_after
    print(f"   Removed {n_removed:,} Excluded cells → {n_after:,} retained")

    # Save
    print(f"   Writing: {FILTERED_H5AD}")
    t1 = time.time()
    adata.write_h5ad(FILTERED_H5AD)
    print(f"   Saved in {time.time()-t1:.0f}s")

    # Extract arrays
    umap   = np.array(adata.obsm["X_umap"])
    level1 = np.array(adata.obs["celltype_level1"])
    level2 = np.array(adata.obs["celltype_level2"])
    study  = np.array(adata.obs["study"])
    n_cells = n_after

    del adata
    gc.collect()

print(f"   n_cells:  {n_cells:,}")
print(f"   level1:   {len(np.unique(level1))} categories")
print(f"   level2:   {len(np.unique(level2))} categories")
print(f"   studies:  {len(np.unique(study))} studies")
