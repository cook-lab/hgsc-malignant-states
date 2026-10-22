#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/9a_atlas_celltypelevel1_v2.py, the script that ran.
# Kept: the split of hgsc_atlas_final.h5ad into celltype_h5ad/hgsc_atlas_final_<lineage>.h5ad (13 lineage subsets,
# read by the epithelial, CNV, deconvolution and Fig 6 steps). Removed: the level-1 figure panels A-E, their style and palette.
# The untrimmed script is in the lab project folder.
"""
Step 9a v2 — Atlas celltype_level1 Figures + Per-Compartment H5AD Split
========================================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Speed-optimized vs v1:
  - Pre-computed RGBA color array (skip 2M hex→RGBA in matplotlib scatter)
  - Panel E background downsampled to 200k points (10× less render work)
  - Sorted index slicing for efficient HDF5 reads during split

New functionality:
  - Splits hgsc_atlas_final.h5ad into 13 per-celltype_level1 h5ad files
    via backed-mode row slicing (safe for 96 GB RAM + 96 GB file)

Panels (unchanged from v1):
  A. UMAP direct-labeled (centroid shadowtext, no legend)
  B. UMAP with legend (clean, no on-plot labels)
  C. Cell count bar chart (horizontal, ordered)
  D. Composition per study (stacked bar, 13 studies)
  E. Per-compartment highlight UMAPs (4×4 grid)

Source:  hgsc_atlas_final.h5ad (backed mode)
Output:  output/09_final_figures/  (figures)
         celltype_h5ad/            (per-compartment h5ad)

Usage:
    python scripts/9a_atlas_celltypelevel1_v2.py
    python scripts/9a_atlas_celltypelevel1_v2.py --figures-only
    python scripts/9a_atlas_celltypelevel1_v2.py --split-only
"""

import argparse
import gc
import os
import time

import anndata as ad
import numpy as np
import pandas as pd

# ============================================================================
# PATHS
# ============================================================================

BASE      = "/Volumes/CookLab/Sarah/202601_final_atlas"
ATLAS     = os.path.join(BASE, "hgsc_atlas_final.h5ad")
SPLIT_DIR = os.path.join(BASE, "celltype_h5ad")

# ============================================================================
# COOK LAB CELL TYPE PALETTE v1.2  (13 celltype_level1 compartments)
# ============================================================================

LEVEL1_ORDER = [
    "Epithelial",
    "Mesothelial",
    "Fibroblast",
    "Smooth muscle",
    "Pericyte",
    "Endothelial",
    "T/NK cell",
    "B cell",
    "Plasma cell",
    "Macrophage",
    "DC",
    "Neutrophil",
    "Mast cell",
]

# ============================================================================
# FIGURE HELPERS
# ============================================================================

def ct_to_filename(ct):
    """Epithelial → epithelial, T/NK cell → t_nk_cell, etc."""
    return ct.lower().replace("/", "_").replace(" ", "_")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="9a v2: figures + split")
    parser.add_argument("--figures-only", action="store_true",
                        help="Generate figures only (skip h5ad split)")
    parser.add_argument("--split-only", action="store_true",
                        help="Split h5ad only (skip figure generation)")
    args = parser.parse_args()

    do_figures = not args.split_only
    do_split   = not args.figures_only

    t0 = time.time()
    print("=" * 65)
    print("  9a v2 — Atlas celltype_level1 (speed-optimized)")
    print("=" * 65)

    # ------------------------------------------------------------------
    # Phase 1: Load atlas backed, extract metadata
    # ------------------------------------------------------------------
    print("\n  [LOAD] Opening atlas (backed)...", flush=True)
    adata = ad.read_h5ad(ATLAS, backed="r")
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    print(f"    {n_cells:,} cells × {n_genes:,} genes")
    print(f"    X format: {type(adata.X).__name__}")
    print(f"    Opened in {time.time()-t0:.0f}s")

    # Extract arrays we need for figures
    print("    Extracting UMAP + metadata...", flush=True)
    umap   = np.array(adata.obsm["X_umap"])
    level1 = adata.obs["celltype_level1"].values.astype(str)
    study  = adata.obs["study"].values.astype(str)

    # Per-celltype indices (shared by figures + split)
    ct_indices = {}
    for ct in LEVEL1_ORDER:
        ct_indices[ct] = np.where(level1 == ct)[0]

    counts = pd.Series(level1).value_counts()
    rng = np.random.RandomState(42)
    order = rng.permutation(n_cells)

    print(f"    Metadata extracted in {time.time()-t0:.0f}s")

    # ==================================================================
    # Phase 3: Split H5AD by celltype_level1
    # ==================================================================
    if do_split:
        os.makedirs(SPLIT_DIR, exist_ok=True)
        t_split = time.time()
        print(f"\n  [SPLIT] Writing per-celltype h5ad → {SPLIT_DIR}")
        print(f"    Source X format: {type(adata.X).__name__}")

        for ct in LEVEL1_ORDER:
            ts = time.time()
            idx = np.sort(ct_indices[ct])   # sorted for sequential HDF5 reads
            n_ct = len(idx)
            fname = f"hgsc_atlas_final_{ct_to_filename(ct)}.h5ad"
            out_path = os.path.join(SPLIT_DIR, fname)

            print(f"\n    {ct}: {n_ct:,} cells", flush=True)
            print(f"      Slicing backed adata...", flush=True)
            sub = adata[idx, :].to_memory()

            print(f"      Writing {fname}...", flush=True)
            sub.write_h5ad(out_path)

            size_mb = os.path.getsize(out_path) / (1024 ** 2)
            print(f"      Done: {size_mb:,.0f} MB  ({time.time()-ts:.0f}s)")

            del sub
            gc.collect()

        elapsed_split = time.time() - t_split
        print(f"\n    Split done in {elapsed_split:.0f}s")

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------
    adata.file.close()
    del adata
    gc.collect()

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    print(f"\n{'=' * 65}")
    print(f"  9a v2 COMPLETE ({elapsed:.0f}s)")
    if do_split:
        print(f"  Split h5ad: {SPLIT_DIR}")
        for ct in LEVEL1_ORDER:
            fname = f"hgsc_atlas_final_{ct_to_filename(ct)}.h5ad"
            print(f"    {fname}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
