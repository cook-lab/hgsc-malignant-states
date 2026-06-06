#!/usr/bin/env python3
"""
Atlas 02 — Step 08 (merge): write celltype_level2 onto the atlas

PURPOSE
    Read the celltype_level1 atlas and add celltype_level2 from the per-compartment
    barcode→level2 CSV maps produced by 08_finalize, then flag low-complexity
    (stressed/low-quality) clusters. Produces the level2 atlas object.

INPUTS
    DATA_ROOT/2026_final_atlas/hgsc_atlas_celltype_level1.h5ad
    output_root/02_annotation/08_celltype_level2/barcode_maps/*.csv

OUTPUTS
    obj("atlas_celltype_l2")  = hgsc_atlas_celltype_level2.h5ad
        (retains all obs cols incl. celltype_level1; adds celltype_level2 + is_low_complexity)

MANUSCRIPT PANEL(S)
    Annotation backend; underpins Fig 1B-E composition and SF4 panels.

RUNTIME TIER
    heavy (loads + rewrites the ~2.3M-cell atlas).
"""

import os
import sys
import time
from pathlib import Path

import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

import numpy as np  # noqa: E402
np.random.seed(SEED)

# ── paths ────────────────────────────────────────────────────
INPUT_H5AD  = path("data_root", "2026_final_atlas", "hgsc_atlas_celltype_level1.h5ad")
MAP_DIR     = path("output_root", "02_annotation", "08_celltype_level2", "barcode_maps")
OUTPUT_H5AD = obj("atlas_celltype_l2")

print("=" * 60)
print("Step 08 — Merge celltype_level2 into atlas")
print("=" * 60)

# ── 1. Concatenate all barcode-to-level2 maps ───────────────
print("\n1. Loading barcode-to-level2 maps...")
map_files = sorted(f for f in os.listdir(MAP_DIR)
                   if f.endswith(".csv") and not f.startswith("."))
print(f"   Found {len(map_files)} map files")

dfs = []
for f in map_files:
    df = pd.read_csv(os.path.join(MAP_DIR, f))
    compartment = f.replace("_barcode_to_level2.csv", "")
    print(f"   {compartment:20s}: {len(df):>8,} barcodes, "
          f"{df.celltype_level2.nunique()} level2 categories")
    dfs.append(df)

merged = pd.concat(dfs, ignore_index=True)
print(f"\n   Total merged barcodes: {len(merged):,}")
assert merged["barcode"].is_unique, (
    f"Duplicate barcodes found! {merged['barcode'].duplicated().sum()} duplicates"
)
barcode_to_level2 = merged.set_index("barcode")["celltype_level2"]

# ── 2. Load the celltype_level1 atlas ────────────────────────
print("\n2. Loading celltype_level1 atlas...")
t0 = time.time()
adata = ad.read_h5ad(INPUT_H5AD)
print(f"   Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes  ({time.time()-t0:.0f}s)")

# ── 3. Map celltype_level2 onto atlas ────────────────────────
print("\n3. Mapping celltype_level2 onto atlas...")
adata.obs["celltype_level2"] = barcode_to_level2.reindex(adata.obs.index).values
n_mapped = adata.obs["celltype_level2"].notna().sum()
n_missing = adata.obs["celltype_level2"].isna().sum()
print(f"   Mapped:  {n_mapped:,} ({100 * n_mapped / adata.n_obs:.1f}%)")
print(f"   Missing: {n_missing:,} ({100 * n_missing / adata.n_obs:.1f}%)")
if n_missing > 0:
    print("\n   WARNING: Some cells have no level2 mapping!")
    missing_mask = adata.obs["celltype_level2"].isna()
    print(adata.obs.loc[missing_mask, "celltype_level1"].value_counts().to_string())
    adata.obs["celltype_level2"] = adata.obs["celltype_level2"].fillna("Unmapped")
adata.obs["celltype_level2"] = pd.Categorical(adata.obs["celltype_level2"])

# ── 4. Summary ───────────────────────────────────────────────
print("\n4. celltype_level2 distribution:")
vc = adata.obs["celltype_level2"].value_counts()
for label, count in vc.items():
    print(f"   {label:45s} {count:>8,}  ({100 * count / adata.n_obs:5.1f}%)")
print(f"\n   Total unique level2 labels: {vc.shape[0]}")
n_excluded = (adata.obs["celltype_level2"] == "Excluded").sum()
n_annotated = adata.n_obs - n_excluded
print(f"   Annotated cells: {n_annotated:,} | Excluded: {n_excluded:,}")

# ── 5. Add is_low_complexity flag ─────────────────────────────
print("\n5. Adding is_low_complexity flag...")
LOW_COMPLEXITY_LABELS = [
    "Stressed Epithelial (ribosomal)",
    "Stressed Epithelial (small)",
    "Stressed Macrophage",
    "Stressed Fibroblast",
    "Low-complexity Plasma Cell",
    "Heat-shock T/NK Cell",
    "Stressed Endothelial",
    "Low-complexity DC",
    "Stressed",
    "Stressed Mesothelial",
    "Ribosomal-high Fibroblasts",
    "Stressed B cell",
    "Stress-activated B cell",
    "Ribosome-high plasma cell",
    "Heat shock plasma cell",
    "Heat shock smooth muscle cell",
    "Ribosome-high smooth muscle cell",
    "Low-quality T Cells",
    "Low-complexity pericyte",
]
adata.obs["is_low_complexity"] = adata.obs["celltype_level2"].isin(LOW_COMPLEXITY_LABELS)
n_lc = adata.obs["is_low_complexity"].sum()
print(f"   Low-complexity cells: {n_lc:,} ({100 * n_lc / adata.n_obs:.1f}%)")

# ── 6. Save ──────────────────────────────────────────────────
print(f"\n6. Saving: {OUTPUT_H5AD}")
t0 = time.time()
adata.write_h5ad(OUTPUT_H5AD)
print(f"   File size: {os.path.getsize(OUTPUT_H5AD) / 1e9:.1f} GB  ({time.time()-t0:.0f}s)")

print(f"\n{'=' * 60}")
print("DONE — Step 08: celltype_level2 atlas created")
print(f"{'=' * 60}")
