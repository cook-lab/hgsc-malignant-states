#!/usr/bin/env python3
"""
Atlas 02 — Step 08c: post-review celltype_level2 label renames

PURPOSE
    After the 10b naming review (marker + GSEA analysis), apply 19 celltype_level2
    label renames atomically to (1) the level2 atlas obs, (2) the barcode→level2
    CSV maps, and (3) the per-label marker CSV filenames. Idempotent.

INPUTS
    obj("atlas_celltype_l2")  = hgsc_atlas_celltype_level2.h5ad
    output_root/02_annotation/08_celltype_level2/barcode_maps/*.csv
    output_root/02_annotation/08_celltype_level2/markers/*.csv

OUTPUTS
    rewrites the above in place with the reviewed labels

MANUSCRIPT PANEL(S)
    Annotation backend; fixes the level2 vocabulary used by Fig 1B-E and SF4.

RUNTIME TIER
    heavy (loads + rewrites the level2 atlas).

NOTE ON LABEL NAMING
    The RENAME_MAP values are the VALIDATED level2 cell-type annotation strings
    (they key the barcode-maps and marker filenames). The manuscript "Transitioning"
    -> "Intermediate" standardisation applies to the epithelial POLARIZATION display
    labels (SecA / Intermediate / SecB; see 12d and the NMF schema), NOT to the
    level2 cell-type name "Transitioning epithelial cell". The level2 strings are
    preserved verbatim here to keep the validated annotation schema intact.
"""

import os
import sys
import shutil
from datetime import datetime
from pathlib import Path

import anndata as ad
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

import numpy as np  # noqa: E402
np.random.seed(SEED)

# ============================================================================
# PATHS (resolved via central config)
# ============================================================================

ATLAS_H5AD  = obj("atlas_celltype_l2")
BARCODE_DIR = path("output_root", "02_annotation", "08_celltype_level2", "barcode_maps")
MARKER_DIR  = path("output_root", "02_annotation", "08_celltype_level2", "markers")

# ============================================================================
# RENAME MAP — old label → new label  (19 renames from 10b review)
# ============================================================================

RENAME_MAP = {
    "Conventional dendritic cell":           "Conventional dendritic cell type 1",
    "ECM-producing smooth muscle cell":      "Inflammatory fibroblast-like smooth muscle cell",
    "Homeostatic smooth muscle":             "Stress-response smooth muscle cell",
    "Germinal centre B cell":               "Cycling B cell",
    "Naive B cell":                          "IFN-activated B cell",
    "Activated capillary endothelial cell":  "Venous endothelial cell",
    "Complement secreting fibroblast":       "PI16+ universal fibroblast",
    "Fibroblast":                            "Ovarian steroidogenic cell",
    "Inflammatory cancer-associated fibroblast": "Hypoxic inflammatory cancer-associated fibroblast",
    "Neuronal cell":                         "Schwann cell",
    "Myeloid-derived dendritic cell":        "Conventional dendritic cell type 2",
    "Perivascular macrophage":               "Cycling C1Q+ tissue-resident macrophage",
    "CD8 memory T cell":                     "CD8 effector/exhausted T cell",
    "Innate lymphoid cell":                  "Hematopoietic stem cell",
    "T/NK cell":                             "Metallothionein-high stress-response T cell",
    "Epithelial cell_1":                     "Secretory epithelial cell",
    "Epithelial cell_2":                     "Stress-response secretory epithelial cell",
    "Transitioning secretory epithelial cell": "Cycling secretory epithelial cell",
    "Proliferative epithelial cell":         "Transitioning epithelial cell",
}

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
    return label.replace(" ", "_").replace("/", "-")


# ============================================================================
# 1. RENAME IN H5AD
# ============================================================================

def rename_h5ad():
    print("\n1. Renaming labels in atlas h5ad...")
    adata = ad.read_h5ad(ATLAS_H5AD)
    col_original = adata.obs["celltype_level2"].astype(str).copy()
    old_labels = set(col_original.unique())

    rename_targets = set(RENAME_MAP.values())
    old_only_keys = [k for k in RENAME_MAP if k not in rename_targets]
    old_only_missing = [k for k in old_only_keys if k not in old_labels]
    if old_only_missing:
        if len(old_only_missing) == len(old_only_keys) and all(
            v in old_labels for v in rename_targets
        ):
            print("   SKIP: Labels already renamed (h5ad is up-to-date)")
            return True
        missing = [k for k in RENAME_MAP if k not in old_labels]
        print(f"   ERROR: Labels not found in h5ad: {missing}")
        sys.exit(1)

    new_col = col_original.map(lambda x: RENAME_MAP.get(x, x))
    adata.obs["celltype_level2"] = pd.Categorical(new_col)

    new_labels = set(adata.obs["celltype_level2"].unique())
    leftover = [k for k in RENAME_MAP if k in new_labels and k not in rename_targets]
    if leftover:
        print(f"   ERROR: Old labels still present after rename: {leftover}")
        sys.exit(1)
    missing_new = set(RENAME_MAP.values()) - new_labels
    if missing_new:
        print(f"   WARNING: Expected new labels not found: {missing_new}")

    adata.write_h5ad(ATLAS_H5AD)
    n_renamed = (col_original != adata.obs["celltype_level2"].astype(str)).sum()
    print(f"   Renamed {n_renamed:,} cell labels; new label count: "
          f"{adata.obs['celltype_level2'].nunique()}")
    for old, new in RENAME_MAP.items():
        n = (col_original == old).sum()
        print(f"   {old:50s} → {new:55s} ({n:>7,} cells)")
    return True


# ============================================================================
# 2. RENAME IN BARCODE MAP CSVs
# ============================================================================

def rename_barcode_maps():
    print("\n2. Renaming labels in barcode map CSVs...")
    for comp_key in sorted(set(LABEL_TO_COMPARTMENT.values())):
        csv_path = os.path.join(BARCODE_DIR, f"{comp_key}_barcode_to_level2.csv")
        if not os.path.exists(csv_path):
            print(f"   WARNING: {csv_path} not found — skipping")
            continue
        df = pd.read_csv(csv_path)
        n_before = df["celltype_level2"].nunique()
        df["celltype_level2"] = df["celltype_level2"].replace(RENAME_MAP)
        df.to_csv(csv_path, index=False)
        print(f"   {comp_key}: {n_before} → {df['celltype_level2'].nunique()} unique labels")
    return True


# ============================================================================
# 3. RENAME MARKER CSV FILES (two-phase move to avoid collisions)
# ============================================================================

def rename_marker_csvs():
    print("\n3. Renaming marker CSV files...")
    file_renames = {}
    for old_label, new_label in RENAME_MAP.items():
        comp_key = LABEL_TO_COMPARTMENT[old_label]
        old_path = os.path.join(MARKER_DIR, f"{comp_key}_{_safe_name(old_label)}.csv")
        new_path = os.path.join(MARKER_DIR, f"{comp_key}_{_safe_name(new_label)}.csv")
        if os.path.exists(old_path):
            file_renames[old_path] = new_path
        elif os.path.exists(new_path):
            print(f"   SKIP (already renamed): {os.path.basename(new_path)}")
        else:
            print(f"   WARNING: {old_path} not found — skipping")

    if not file_renames:
        print("   No files to rename (all already renamed or missing)")
        return True

    staging = {}
    for old_path, new_path in file_renames.items():
        staging_path = old_path + ".staging"
        shutil.move(old_path, staging_path)
        staging[staging_path] = new_path
    for staging_path, new_path in staging.items():
        shutil.move(staging_path, new_path)
        print(f"   {os.path.basename(staging_path.replace('.staging', '')):60s} → "
              f"{os.path.basename(new_path)}")
    print(f"   Renamed {len(file_renames)} files")
    return True


def main():
    print("=" * 70)
    print("  08c_rename_level2_labels — Post-10b review renames")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)
    print(f"\n  {len(RENAME_MAP)} label renames to apply")
    rename_h5ad()
    rename_barcode_maps()
    rename_marker_csvs()
    print(f"\n{'='*70}\n  DONE — All renames applied successfully\n{'='*70}")


if __name__ == "__main__":
    main()
