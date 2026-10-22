#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/9e_celltypelevel2_umaps.py, the script that ran.
# Kept: the per-lineage UMAP on X_scanvi (n = 15, min_dist 0.3) saved as obsm X_umap_local; the epithelial one is the
# epithelial UMAP (meta.parquet UMAP1/UMAP2). Removed: the level-2 UMAP plots with their palettes, style and plot helpers.
# The untrimmed script is in the lab project folder.
"""
Step 9e — Per-Compartment celltype_level2 UMAPs
=================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

For each of the 13 celltype_level1 compartments:
  1. Load per-compartment h5ad from celltype_h5ad/
  2. Recompute local UMAP from X_scanvi (neighbors n=15, min_dist=0.3)
  3. Plot UMAP colored by celltype_level2 subtypes (level1 divergent palette)
  4. Save recomputed UMAP as X_umap_local in the h5ad

Source:  celltype_h5ad/hgsc_atlas_final_{compartment}.h5ad
Output:  output/09_final_figures/9e_{compartment}_level2_umap.svg/pdf

Usage:
    python scripts/9e_celltypelevel2_umaps.py
    python scripts/9e_celltypelevel2_umaps.py --compartment "B cell"
    python scripts/9e_celltypelevel2_umaps.py --skip-recompute
    python scripts/9e_celltypelevel2_umaps.py --figures-only
"""

import argparse
import gc
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# ============================================================================
# PATHS
# ============================================================================

BASE      = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD_DIR  = os.path.join(BASE, "celltype_h5ad")

# ============================================================================
# PALETTES & CONSTANTS
# ============================================================================

# Processing order: smallest compartment first
COMPARTMENT_ORDER = [
    "Neutrophil",
    "Mast cell",
    "DC",
    "Pericyte",
    "Plasma cell",
    "Smooth muscle",
    "B cell",
    "Endothelial",
    "Mesothelial",
    "Fibroblast",
    "Macrophage",
    "T/NK cell",
    "Epithelial",
]

# Files above this size (bytes) use backed-mode loading
BACKED_THRESHOLD = 3 * 1024**3  # 3 GB


# ============================================================================
# HELPERS
# ============================================================================

def ct_to_filename(ct):
    """Epithelial → epithelial, T/NK cell → t_nk_cell, etc."""
    return ct.lower().replace("/", "_").replace(" ", "_")


# ============================================================================
# CORE: plot_compartment_umap
# ============================================================================

def plot_compartment_umap(compartment, skip_recompute=False, figures_only=False):
    """Load compartment h5ad, recompute local UMAP, plot, and save."""
    fname = f"hgsc_atlas_final_{ct_to_filename(compartment)}.h5ad"
    h5ad_path = os.path.join(H5AD_DIR, fname)

    if not os.path.exists(h5ad_path):
        print(f"    WARNING: {fname} not found, skipping")
        return

    file_size = os.path.getsize(h5ad_path)
    use_backed = file_size > BACKED_THRESHOLD
    t0 = time.time()

    # ------------------------------------------------------------------
    # Load
    # ------------------------------------------------------------------
    if use_backed:
        print(f"    Loading (backed, {file_size / 1024**3:.1f} GB)...", flush=True)
        bdata = ad.read_h5ad(h5ad_path, backed="r")
        n_cells = bdata.n_obs

        # Extract only what we need
        scanvi = np.array(bdata.obsm["X_scanvi"])
        level2 = bdata.obs["celltype_level2"].values.astype(str)
        barcodes = bdata.obs.index.values.copy()

        # Check for existing local UMAP
        existing_umap = None
        if skip_recompute and "X_umap_local" in bdata.obsm:
            existing_umap = np.array(bdata.obsm["X_umap_local"])
        elif skip_recompute and "X_umap" in bdata.obsm:
            existing_umap = np.array(bdata.obsm["X_umap"])

        bdata.file.close()
        del bdata
        gc.collect()

        # Build minimal AnnData for UMAP computation
        adata = ad.AnnData(
            obs=pd.DataFrame({"celltype_level2": level2}, index=barcodes),
        )
        adata.obsm["X_scanvi"] = scanvi

        if existing_umap is not None:
            adata.obsm["X_umap"] = existing_umap

    else:
        print(f"    Loading ({file_size / 1024**3:.1f} GB)...", flush=True)
        adata = ad.read_h5ad(h5ad_path)
        n_cells = adata.n_obs

    print(f"      {n_cells:,} cells, loaded in {time.time()-t0:.0f}s")

    # ------------------------------------------------------------------
    # Recompute local UMAP
    # ------------------------------------------------------------------
    if not skip_recompute:
        print("      Computing neighbors + UMAP (X_scanvi)...", flush=True)
        t_umap = time.time()
        sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=15, random_state=0)
        sc.tl.umap(adata, min_dist=0.3, random_state=42)
        print(f"      UMAP done in {time.time()-t_umap:.0f}s")
    else:
        if "X_umap" not in adata.obsm:
            print("      WARNING: No UMAP found and --skip-recompute set, skipping")
            del adata
            gc.collect()
            return
        print("      Using existing UMAP coordinates")

    umap_coords = np.array(adata.obsm["X_umap"])

    # ------------------------------------------------------------------
    # Save X_umap_local back to h5ad
    # ------------------------------------------------------------------
    if not skip_recompute and not figures_only:
        print("      Saving X_umap_local to h5ad...", flush=True)
        t_save = time.time()
        if use_backed:
            # Reopen file to write the new obsm
            import h5py
            with h5py.File(h5ad_path, "a") as f:
                obsm = f["obsm"] if "obsm" in f else f.create_group("obsm")
                key = "X_umap_local"
                if key in obsm:
                    del obsm[key]
                obsm.create_dataset(key, data=umap_coords)
            print(f"      Saved X_umap_local ({time.time()-t_save:.0f}s)")
        else:
            adata.obsm["X_umap_local"] = umap_coords
            adata.write_h5ad(h5ad_path)
            print(f"      Saved X_umap_local ({time.time()-t_save:.0f}s)")

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------
    del adata, umap_coords
    gc.collect()

    elapsed = time.time() - t0
    print(f"      {compartment} complete ({elapsed:.0f}s)")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="9e: Per-compartment celltype_level2 UMAPs")
    parser.add_argument("--compartment", type=str, default=None,
                        help="Run a single compartment (e.g. 'B cell')")
    parser.add_argument("--skip-recompute", action="store_true",
                        help="Use existing X_umap instead of recomputing")
    parser.add_argument("--figures-only", action="store_true",
                        help="Skip saving X_umap_local back to h5ad")
    args = parser.parse_args()

    t0 = time.time()
    print("=" * 65)
    print("  9e — Per-Compartment celltype_level2 UMAPs")
    print("=" * 65)

    if args.compartment:
        compartments = [args.compartment]
    else:
        compartments = COMPARTMENT_ORDER

    n_total = len(compartments)
    for i, comp in enumerate(compartments, 1):
        print(f"\n  [{i}/{n_total}] {comp}")
        print(f"  {'-' * 50}")
        plot_compartment_umap(comp,
                              skip_recompute=args.skip_recompute,
                              figures_only=args.figures_only)

    elapsed = time.time() - t0
    print(f"\n{'=' * 65}")
    print(f"  9e COMPLETE ({elapsed:.0f}s)")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
