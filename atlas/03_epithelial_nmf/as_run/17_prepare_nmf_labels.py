#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/17_prepare_nmf_labels.py, the script that ran.
# Kept: the Factor 2 percentile labels (SecA, Transitioning = Intermediate, SecB, ciliated) written to celltype_nmf_mapping.csv,
# read by the UCell export and the TCGA reference. Removed: the label-summary CSV and the in-place celltype_nmf column in
# hgsc_atlas_final.h5ad, which nothing reads. The untrimmed script is in the lab project folder.
"""
Step 17 — Prepare NMF-based epithelial labels for LIANA re-run
================================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Creates a new `celltype_nmf` column in hgsc_atlas_final.h5ad that replaces
the 6 Leiden-based epithelial celltype_level2 labels with NMF-derived
SecA / Transitioning / SecB / Ciliated labels (matching manuscript figures).

Non-epithelial cells retain their existing celltype_level2 labels unchanged.

NMF Schema (from 00_prepare_data.py):
  - Factor_2 from 11d NMF, thresholds computed on non-ciliated epithelial cells
  - SecA:          Factor_2 < p50
  - Transitioning: p50 ≤ Factor_2 < p75
  - SecB:          Factor_2 ≥ p75
  - Ciliated:      celltype_level2 == "Ciliated epithelial cell"

Usage:
    python scripts/17_prepare_nmf_labels.py

Output:
    - Updated hgsc_atlas_final.h5ad with `celltype_nmf` obs column
    - output/17_nmf_labels/celltype_nmf_mapping.csv
    - output/17_nmf_labels/celltype_nmf_summary.csv
"""

import os
import warnings

import numpy as np
import pandas as pd
import anndata as ad

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────
PROJECT_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas"
ATLAS_H5AD  = os.path.join(PROJECT_DIR, "hgsc_atlas_final.h5ad")
NMF_USAGE   = os.path.join(PROJECT_DIR, "output", "11d_epithelial_nmf",
                            "11d_nmf_usage.csv")
OUT_DIR     = os.path.join(PROJECT_DIR, "output", "17_nmf_labels")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================
# 1. LOAD METADATA
# ============================================================================

print("=" * 60)
print("  Step 17 — Prepare NMF-based epithelial labels")
print("=" * 60)

print("\n[1] Loading atlas metadata (backed mode)...")
adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
obs = adata.obs[["celltype_level1", "celltype_level2"]].copy()
n_total = len(obs)
print(f"    Atlas: {n_total:,} cells")

# ============================================================================
# 2. LOAD NMF FACTOR_2
# ============================================================================

print("\n[2] Loading NMF Factor_2 usage...")
nmf_usage = pd.read_csv(NMF_USAGE, index_col=0)
f2 = nmf_usage["Factor_2"]
print(f"    NMF data: {len(f2):,} epithelial cells")

# ============================================================================
# 3. COMPUTE THRESHOLDS & ASSIGN LABELS
# ============================================================================

print("\n[3] Computing NMF thresholds and assigning labels...")

# Identify epithelial cells in the atlas
is_epithelial = obs["celltype_level1"] == "Epithelial"
epi_barcodes = obs.index[is_epithelial]
print(f"    Epithelial cells in atlas: {is_epithelial.sum():,}")

# Match epithelial barcodes with NMF data
shared = epi_barcodes.intersection(f2.index)
print(f"    Matched with NMF data: {len(shared):,}")

# Get Factor_2 for matched epithelial cells
f2_matched = f2.loc[shared]

# Identify ciliated among matched cells
is_ciliated = obs.loc[shared, "celltype_level2"] == "Ciliated epithelial cell"
non_cil_f2 = f2_matched[~is_ciliated]

# Compute percentiles on non-ciliated epithelial cells
p50 = np.percentile(non_cil_f2, 50)
p75 = np.percentile(non_cil_f2, 75)
print(f"    NMF thresholds: p50={p50:.4f}, p75={p75:.4f}")

# Assign NMF labels for epithelial cells
epi_labels = pd.Series("SecA epithelium", index=shared)
epi_labels[f2_matched >= p50] = "Transitioning epithelium"
epi_labels[f2_matched >= p75] = "SecB epithelium"
epi_labels[is_ciliated] = "Ciliated epithelial cell"

# Print epithelial distribution
for g in ["SecA epithelium", "Transitioning epithelium", "SecB epithelium",
           "Ciliated epithelial cell"]:
    n = (epi_labels == g).sum()
    print(f"    {g}: {n:,} ({n/len(epi_labels)*100:.1f}%)")

# ============================================================================
# 4. BUILD FULL celltype_nmf COLUMN
# ============================================================================

print("\n[4] Building full celltype_nmf column...")

# Start with celltype_level2 for all cells (convert to string to allow new labels)
celltype_nmf = obs["celltype_level2"].astype(str).copy()

# Override epithelial cells with NMF labels
celltype_nmf.loc[shared] = epi_labels

# Handle any epithelial cells NOT in NMF data (shouldn't happen, but be safe)
epi_not_in_nmf = epi_barcodes.difference(shared)
if len(epi_not_in_nmf) > 0:
    print(f"    WARNING: {len(epi_not_in_nmf)} epithelial cells not in NMF data")
    # Keep their celltype_level2 label
    print(f"    These will retain their celltype_level2 labels")

print(f"    Total unique labels: {celltype_nmf.nunique()}")
print(f"    Non-epithelial labels: {celltype_nmf[~is_epithelial].nunique()}")

# ============================================================================
# 5. SAVE MAPPING CSV
# ============================================================================

print("\n[5] Saving mapping CSV...")

mapping = pd.DataFrame({
    "celltype_level2": obs["celltype_level2"],
    "celltype_nmf": celltype_nmf,
    "celltype_level1": obs["celltype_level1"],
})
mapping.to_csv(os.path.join(OUT_DIR, "celltype_nmf_mapping.csv"))
print(f"    Saved: celltype_nmf_mapping.csv ({len(mapping):,} rows)")

print(f"\n{'='*60}")
print("  Step 17 complete!")
print(f"  Output: {OUT_DIR}")
print(f"{'='*60}")
