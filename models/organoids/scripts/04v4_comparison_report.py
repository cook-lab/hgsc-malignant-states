#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202603_organoids/scripts/04v4_comparison_report.py, the script that ran. Kept:
# organoid_ucell_scores_v4_with_atlas_class.csv, whose atlas_class gives the "fewer than 1%" of the Results.
# Removed: the figures, the per-model summary table and the HTML report. The untrimmed script, which the docstring
# below describes, is in the lab project folder.

"""
04v4 — Organoid vs Atlas UCell Comparison Report
=================================================
Cook Lab · 2026

Contextualizes organoid UCell scores (v4, aligned gene space) against
atlas secretory epithelial distributions to guide cutoff decisions.

Input:
  - organoids: output/02_Secretory_Polarization_v4_UCell_aligned/organoid_ucell_scores_v4_aligned.csv
  - atlas:     /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_ucell_scores.csv
  - atlas:     /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_secretory_metadata.csv

Output:
  - output/02_Secretory_Polarization_v4_UCell_aligned/04v4_atlas_organoid_comparison.html
"""

import os
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
ORG_DIR   = "/Volumes/CookLab/Sarah/organoids"
ATLAS_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas"
OUT_DIR   = os.path.join(ORG_DIR, "output", "02_Secretory_Polarization_v4_UCell_aligned")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

print("=" * 60)
print("  Organoid vs Atlas UCell Comparison Report")
print("=" * 60)

print("\n[1] Loading data...")
org = pd.read_csv(os.path.join(OUT_DIR, "organoid_ucell_scores_v4_aligned.csv"))
atlas_scores = pd.read_csv(os.path.join(ATLAS_DIR, "atlas_ucell_scores.csv"))
atlas_meta = pd.read_csv(os.path.join(ATLAS_DIR, "atlas_secretory_metadata.csv"), index_col=0)
atlas = atlas_meta.join(atlas_scores.set_index("barcode"), how="inner")

print(f"    Organoid cells: {len(org):,}")
print(f"    Atlas secretory cells: {len(atlas):,}")

models = sorted(org["model"].unique())
print(f"    Organoid models: {', '.join(models)}")

# Compute atlas percentile for each organoid cell
atlas_secA_vals = atlas["SecA_UCell"].values
atlas_secB_vals = atlas["SecB_UCell"].values
atlas_pol_vals  = atlas["sec_polarization"].values

org["atlas_pctile_SecA"] = org["SecA_UCell"].apply(
    lambda x: np.searchsorted(np.sort(atlas_secA_vals), x) / len(atlas_secA_vals) * 100)
org["atlas_pctile_SecB"] = org["SecB_UCell"].apply(
    lambda x: np.searchsorted(np.sort(atlas_secB_vals), x) / len(atlas_secB_vals) * 100)
org["atlas_pctile_pol"] = org["sec_polarization"].apply(
    lambda x: np.searchsorted(np.sort(atlas_pol_vals), x) / len(atlas_pol_vals) * 100)

# Classify organoid cells based on atlas NMF percentile boundaries
# Atlas: SecA = bottom 50% polarization, Trans = 50-75%, SecB = top 25%
atlas_pol_sorted = np.sort(atlas_pol_vals)
atlas_p50 = atlas_pol_sorted[int(len(atlas_pol_sorted) * 0.50)]
atlas_p75 = atlas_pol_sorted[int(len(atlas_pol_sorted) * 0.75)]

org["atlas_class"] = "Transitioning"
org.loc[org["sec_polarization"] < atlas_p50, "atlas_class"] = "SecA"
org.loc[org["sec_polarization"] >= atlas_p75, "atlas_class"] = "SecB"

# Also save organoid scores with atlas classification
org.to_csv(os.path.join(OUT_DIR, "organoid_ucell_scores_v4_with_atlas_class.csv"), index=False)

print(f"\n{'='*60}")
print(f"  Report complete!")
print(f"{'='*60}")
