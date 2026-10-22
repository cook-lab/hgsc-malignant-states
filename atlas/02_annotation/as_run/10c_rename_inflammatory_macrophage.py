#!/usr/bin/env python3
"""
10c_rename_inflammatory_macrophage.py — Rename Inflammatory → Monocyte-derived
================================================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Based on 10c macrophage validation results, renames "Inflammatory macrophage"
to "Monocyte-derived macrophage" in:

  1. hgsc_atlas_final.h5ad          (obs["celltype_level2"])
  2. hgsc_atlas_celltype_level2.h5ad (obs["celltype_level2"])
  3. Barcode map CSV                 (content replacement)
  4. Marker CSV                      (filename rename)

Usage:
    python scripts/10c_rename_inflammatory_macrophage.py
"""

import os
import shutil
from datetime import datetime

import anndata as ad
import pandas as pd

# ============================================================================
# CONFIG
# ============================================================================

PROJECT_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas"

OLD_NAME = "Inflammatory macrophage"
NEW_NAME = "Monocyte-derived macrophage"

ATLAS_FILES = [
    os.path.join(PROJECT_DIR, "hgsc_atlas_final.h5ad"),
    os.path.join(PROJECT_DIR, "hgsc_atlas_celltype_level2.h5ad"),
]

BARCODE_DIR = os.path.join(PROJECT_DIR, "output", "08_celltype_level2", "barcode_maps")
MARKER_DIR  = os.path.join(PROJECT_DIR, "output", "08_celltype_level2", "markers")


def _safe_name(label):
    return label.replace(" ", "_").replace("/", "-")


# ============================================================================
# 1. RENAME IN H5AD FILES
# ============================================================================

def rename_h5ad(path):
    """Rename label in a single h5ad file."""
    print(f"\n  Loading: {os.path.basename(path)}")
    adata = ad.read_h5ad(path)

    col = adata.obs["celltype_level2"].astype(str)
    n_match = (col == OLD_NAME).sum()

    if n_match == 0:
        print(f"    '{OLD_NAME}' not found — skipping")
        del adata
        return

    print(f"    Found {n_match:,} cells with '{OLD_NAME}'")
    col = col.replace(OLD_NAME, NEW_NAME)
    adata.obs["celltype_level2"] = pd.Categorical(col)

    # Verify
    assert (adata.obs["celltype_level2"] == OLD_NAME).sum() == 0
    assert (adata.obs["celltype_level2"] == NEW_NAME).sum() == n_match

    print(f"    Writing back to: {os.path.basename(path)}")
    adata.write_h5ad(path)
    del adata
    print(f"    Done — {n_match:,} cells renamed")


# ============================================================================
# 2. RENAME IN BARCODE MAP
# ============================================================================

def rename_barcode_map():
    """Rename label in macrophage barcode map CSV."""
    csv_path = os.path.join(BARCODE_DIR, "macrophage_barcode_to_level2.csv")
    if not os.path.exists(csv_path):
        print(f"\n  Barcode map not found: {csv_path}")
        return

    print(f"\n  Updating barcode map: {os.path.basename(csv_path)}")
    df = pd.read_csv(csv_path)
    n_match = (df["celltype_level2"] == OLD_NAME).sum()
    df["celltype_level2"] = df["celltype_level2"].replace(OLD_NAME, NEW_NAME)
    df.to_csv(csv_path, index=False)
    print(f"    Renamed {n_match:,} entries")


# ============================================================================
# 3. RENAME MARKER CSV FILE
# ============================================================================

def rename_marker_file():
    """Rename the marker CSV filename."""
    old_file = os.path.join(MARKER_DIR, f"macrophage_{_safe_name(OLD_NAME)}.csv")
    new_file = os.path.join(MARKER_DIR, f"macrophage_{_safe_name(NEW_NAME)}.csv")

    if os.path.exists(old_file):
        print(f"\n  Renaming marker file:")
        print(f"    {os.path.basename(old_file)} → {os.path.basename(new_file)}")
        shutil.move(old_file, new_file)
    elif os.path.exists(new_file):
        print(f"\n  Marker file already renamed: {os.path.basename(new_file)}")
    else:
        print(f"\n  WARNING: marker file not found: {os.path.basename(old_file)}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("  10c — Rename 'Inflammatory macrophage' → 'Monocyte-derived macrophage'")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # 1. H5AD files
    for path in ATLAS_FILES:
        if os.path.exists(path):
            rename_h5ad(path)
        else:
            print(f"\n  SKIP (not found): {path}")

    # 2. Barcode map
    rename_barcode_map()

    # 3. Marker file
    rename_marker_file()

    print(f"\n{'=' * 70}")
    print("  DONE — all data files updated")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
