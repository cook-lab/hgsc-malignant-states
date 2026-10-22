#!/usr/bin/env python3
"""
08c_rename_level2_labels_v2.py — Post-10b review renames
=========================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

After 10b naming review (agent-assisted marker + GSEA analysis), 19
celltype_level2 labels are renamed to better reflect their marker/enrichment
profiles.  This script applies the renames atomically to:

  1. The atlas h5ad  (obs["celltype_level2"])
  2. Barcode map CSVs (output/08_celltype_level2/barcode_maps/)
  3. Marker CSV files (output/08_celltype_level2/markers/) — filename only

Usage:
    python scripts/08c_rename_level2_labels_v2.py

No arguments required.  Prints a summary of all changes.
"""

import os
import shutil
import sys
from datetime import datetime

import anndata as ad
import pandas as pd

# ============================================================================
# PATHS
# ============================================================================

PROJECT_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas"
ATLAS_H5AD  = os.path.join(PROJECT_DIR, "hgsc_atlas_celltype_level2.h5ad")
BARCODE_DIR = os.path.join(PROJECT_DIR, "output", "08_celltype_level2", "barcode_maps")
MARKER_DIR  = os.path.join(PROJECT_DIR, "output", "08_celltype_level2", "markers")

# ============================================================================
# RENAME MAP — old label → new label  (19 renames from 10b review)
# ============================================================================

RENAME_MAP = {
    # Dendritic cell
    "Conventional dendritic cell":           "Conventional dendritic cell type 1",
    # Smooth muscle
    "ECM-producing smooth muscle cell":      "Inflammatory fibroblast-like smooth muscle cell",
    "Homeostatic smooth muscle":             "Stress-response smooth muscle cell",
    # B cell
    "Germinal centre B cell":               "Cycling B cell",
    "Naive B cell":                          "IFN-activated B cell",
    # Endothelial
    "Activated capillary endothelial cell":  "Venous endothelial cell",
    # Fibroblast
    "Complement secreting fibroblast":       "PI16+ universal fibroblast",
    "Fibroblast":                            "Ovarian steroidogenic cell",
    "Inflammatory cancer-associated fibroblast": "Hypoxic inflammatory cancer-associated fibroblast",
    "Neuronal cell":                         "Schwann cell",
    # Macrophage
    "Myeloid-derived dendritic cell":        "Conventional dendritic cell type 2",
    "Perivascular macrophage":               "Cycling C1Q+ tissue-resident macrophage",
    # T/NK cell
    "CD8 memory T cell":                     "CD8 effector/exhausted T cell",
    "Innate lymphoid cell":                  "Hematopoietic stem cell",
    "T/NK cell":                             "Metallothionein-high stress-response T cell",
    # Epithelial
    "Epithelial cell_1":                     "Secretory epithelial cell",
    "Epithelial cell_2":                     "Stress-response secretory epithelial cell",
    "Transitioning secretory epithelial cell": "Cycling secretory epithelial cell",
    "Proliferative epithelial cell":         "Transitioning epithelial cell",
}

# Map old labels to their compartment key (for marker CSV filenames)
LABEL_TO_COMPARTMENT = {
    "Conventional dendritic cell":           "dc",
    "ECM-producing smooth muscle cell":      "smoothmuscle",
    "Homeostatic smooth muscle":             "smoothmuscle",
    "Germinal centre B cell":               "bcell",
    "Naive B cell":                          "bcell",
    "Activated capillary endothelial cell":  "endothelial",
    "Complement secreting fibroblast":       "fibroblast",
    "Fibroblast":                            "fibroblast",
    "Inflammatory cancer-associated fibroblast": "fibroblast",
    "Neuronal cell":                         "fibroblast",
    "Myeloid-derived dendritic cell":        "macrophage",
    "Perivascular macrophage":               "macrophage",
    "CD8 memory T cell":                     "tnkcell",
    "Innate lymphoid cell":                  "tnkcell",
    "T/NK cell":                             "tnkcell",
    "Epithelial cell_1":                     "epithelial",
    "Epithelial cell_2":                     "epithelial",
    "Transitioning secretory epithelial cell": "epithelial",
    "Proliferative epithelial cell":         "epithelial",
}


def _safe_name(label):
    """Convert label to safe filename component (matches step 08 convention)."""
    return label.replace(" ", "_").replace("/", "-")


# ============================================================================
# 1. RENAME IN H5AD
# ============================================================================

def rename_h5ad():
    print("\n1. Renaming labels in atlas h5ad...")
    print(f"   Loading: {ATLAS_H5AD}")
    adata = ad.read_h5ad(ATLAS_H5AD)

    col_original = adata.obs["celltype_level2"].astype(str).copy()
    old_labels = set(col_original.unique())

    # Check if already renamed (idempotent)
    # Since no old key is also a new value, all old keys should be present
    # before rename and absent after.
    rename_targets = set(RENAME_MAP.values())
    old_only_keys = [k for k in RENAME_MAP if k not in rename_targets]
    old_only_missing = [k for k in old_only_keys if k not in old_labels]
    if old_only_missing:
        # If *all* old-only keys are gone AND all new values exist → already done
        if len(old_only_missing) == len(old_only_keys) and all(
            v in old_labels for v in rename_targets
        ):
            print("   SKIP: Labels already renamed (h5ad is up-to-date)")
            del adata
            return True
        else:
            missing = [k for k in RENAME_MAP if k not in old_labels]
            print(f"   ERROR: Labels not found in h5ad: {missing}")
            sys.exit(1)

    # Apply renames (single pass — no collision because no old key == new value)
    new_col = col_original.map(lambda x: RENAME_MAP.get(x, x))
    adata.obs["celltype_level2"] = pd.Categorical(new_col)

    # Validate — all old keys should be gone (none is also a new value)
    new_labels = set(adata.obs["celltype_level2"].unique())
    leftover = [k for k in RENAME_MAP if k in new_labels and k not in rename_targets]
    if leftover:
        print(f"   ERROR: Old labels still present after rename: {leftover}")
        sys.exit(1)

    expected_new = set(RENAME_MAP.values())
    missing_new = expected_new - new_labels
    if missing_new:
        print(f"   WARNING: Expected new labels not found: {missing_new}")

    # Save
    print(f"   Saving: {ATLAS_H5AD}")
    adata.write_h5ad(ATLAS_H5AD)

    # Count renamed cells
    n_renamed = (col_original != adata.obs["celltype_level2"].astype(str)).sum()
    print(f"   Renamed {n_renamed:,} cell labels")
    print(f"   New label count: {adata.obs['celltype_level2'].nunique()}")

    # Print summary
    for old, new in RENAME_MAP.items():
        n = (col_original == old).sum()
        print(f"   {old:50s} → {new:55s} ({n:>7,} cells)")

    del adata
    return True


# ============================================================================
# 2. RENAME IN BARCODE MAP CSVs
# ============================================================================

def rename_barcode_maps():
    print("\n2. Renaming labels in barcode map CSVs...")

    affected_comps = sorted(set(LABEL_TO_COMPARTMENT.values()))
    for comp_key in affected_comps:
        csv_path = os.path.join(BARCODE_DIR, f"{comp_key}_barcode_to_level2.csv")
        if not os.path.exists(csv_path):
            print(f"   WARNING: {csv_path} not found — skipping")
            continue

        df = pd.read_csv(csv_path)
        n_before = df["celltype_level2"].nunique()
        df["celltype_level2"] = df["celltype_level2"].replace(RENAME_MAP)
        n_after = df["celltype_level2"].nunique()
        df.to_csv(csv_path, index=False)
        print(f"   {comp_key}: {n_before} → {n_after} unique labels")

    return True


# ============================================================================
# 3. RENAME MARKER CSV FILES (two-phase move to avoid collisions)
# ============================================================================

def rename_marker_csvs():
    print("\n3. Renaming marker CSV files...")

    # Build file rename map
    file_renames = {}  # old_path → new_path
    for old_label, new_label in RENAME_MAP.items():
        comp_key = LABEL_TO_COMPARTMENT[old_label]
        old_safe = _safe_name(old_label)
        new_safe = _safe_name(new_label)
        old_path = os.path.join(MARKER_DIR, f"{comp_key}_{old_safe}.csv")
        new_path = os.path.join(MARKER_DIR, f"{comp_key}_{new_safe}.csv")
        if os.path.exists(old_path):
            file_renames[old_path] = new_path
        else:
            # Check if already renamed
            if os.path.exists(new_path):
                print(f"   SKIP (already renamed): {os.path.basename(new_path)}")
            else:
                print(f"   WARNING: {old_path} not found — skipping")

    if not file_renames:
        print("   No files to rename (all already renamed or missing)")
        return True

    # Phase 1: move to staging
    staging = {}
    for old_path, new_path in file_renames.items():
        staging_path = old_path + ".staging"
        shutil.move(old_path, staging_path)
        staging[staging_path] = new_path

    # Phase 2: move to final
    for staging_path, new_path in staging.items():
        shutil.move(staging_path, new_path)
        old_name = os.path.basename(staging_path.replace(".staging", ""))
        new_name = os.path.basename(new_path)
        print(f"   {old_name:60s} → {new_name}")

    print(f"   Renamed {len(file_renames)} files")
    return True


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("  08c_rename_level2_labels_v2.py — Post-10b review renames")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    print(f"\n  {len(RENAME_MAP)} label renames to apply:")
    for old, new in RENAME_MAP.items():
        print(f"    {old:50s} → {new}")

    rename_h5ad()
    rename_barcode_maps()
    rename_marker_csvs()

    print(f"\n{'='*70}")
    print("  DONE — All renames applied successfully")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
