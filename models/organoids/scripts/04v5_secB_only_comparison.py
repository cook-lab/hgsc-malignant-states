#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202603_organoids/scripts/04v5_secB_only_comparison.py, the script that ran. Kept:
# organoid_secB_classified_v5.csv, the per-cell SecB UCell by model that figures/figure3/02_... plots (Fig 3b).
# Removed: the threshold check on atlas cells, the per-model summary, the figures and the HTML report. The untrimmed
# script, which the docstring below describes, is in the lab project folder.

"""
04v5 — SecB-Only UCell Comparison: Organoids vs Atlas
=====================================================
Cook Lab · 2026

Focuses exclusively on SecB UCell scores (dropping the weak SecA signal)
to derive a clean, atlas-grounded SecB threshold for organoid classification.

Rationale:
  - SecA UCell is noisy (TF-based genes rank poorly in UCell)
  - SecB UCell alone is the dominant discriminator in the atlas
  - A single SecB threshold is cleaner than the polarization metric

Uses existing v4 aligned scores (UCell computes signatures independently,
so SecB_UCell is identical whether scored alone or alongside SecA).

Input:
  - organoids: output/02_Secretory_Polarization_v4_UCell_aligned/organoid_ucell_scores_v4_aligned.csv
  - atlas:     /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_ucell_scores.csv
  - atlas:     /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_secretory_metadata.csv

Output:
  - output/02_Secretory_Polarization_v5_SecB_only/04v5_secB_only_comparison.html
  - output/02_Secretory_Polarization_v5_SecB_only/model_secB_summary_v5.csv
  - output/02_Secretory_Polarization_v5_SecB_only/organoid_secB_classified_v5.csv
  - figs/02_Secretory_Polarization_v5_SecB_only/*.svg + *.pdf
"""

import os
import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, auc

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
ORG_DIR   = "/Volumes/CookLab/Sarah/organoids"
ATLAS_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas"
OUT_DIR   = os.path.join(ORG_DIR, "output", "02_Secretory_Polarization_v5_SecB_only")
os.makedirs(OUT_DIR, exist_ok=True)

NMF_ORDER = ["SecA epithelium", "Transitioning epithelium", "SecB epithelium"]

SECB_GENES = ["KRT17", "KRT19", "KRT7", "LCN2", "PRSS22", "SLPI", "TACSTD2"]

# ============================================================================
# 1. LOAD DATA
# ============================================================================

print("=" * 60)
print("  04v5 — SecB-Only UCell: Organoids vs Atlas")
print("=" * 60)

print("\n[1] Loading data...")
org = pd.read_csv(os.path.join(
    ORG_DIR, "output", "02_Secretory_Polarization_v4_UCell_aligned",
    "organoid_ucell_scores_v4_aligned.csv"))
atlas_scores = pd.read_csv(os.path.join(ATLAS_DIR, "atlas_ucell_scores.csv"))
atlas_meta = pd.read_csv(os.path.join(ATLAS_DIR, "atlas_secretory_metadata.csv"), index_col=0)
atlas = atlas_meta.join(atlas_scores.set_index("barcode"), how="inner")

print(f"    Organoid cells: {len(org):,}")
print(f"    Atlas secretory cells: {len(atlas):,}")

models = sorted(org["model"].unique())
print(f"    Organoid models: {', '.join(models)}")

# ============================================================================
# 2. ATLAS SecB REFERENCE DISTRIBUTIONS
# ============================================================================

print("\n[2] Computing atlas SecB reference statistics...")

atlas_stats = {}
for label in NMF_ORDER:
    vals = atlas.loc[atlas["celltype_nmf"] == label, "SecB_UCell"]
    atlas_stats[label] = {
        "n": len(vals),
        "mean": vals.mean(),
        "median": vals.median(),
        "sd": vals.std(),
        "p5": vals.quantile(0.05),
        "p10": vals.quantile(0.10),
        "p25": vals.quantile(0.25),
        "p75": vals.quantile(0.75),
        "p90": vals.quantile(0.90),
        "p95": vals.quantile(0.95),
    }
    short = label.replace(" epithelium", "")
    print(f"    {short:15s}  median={vals.median():.4f}  "
          f"IQR=[{vals.quantile(0.25):.4f}, {vals.quantile(0.75):.4f}]  "
          f"p5={vals.quantile(0.05):.4f}  p95={vals.quantile(0.95):.4f}")

# ============================================================================
# 3. DERIVE ATLAS-GROUNDED SecB THRESHOLDS
# ============================================================================

print("\n[3] Deriving SecB thresholds from atlas...")

# Approach 1: ROC-optimal threshold (SecA vs SecB binary, ignoring Transitioning)
secA_vals = atlas.loc[atlas["celltype_nmf"] == "SecA epithelium", "SecB_UCell"].values
secB_vals = atlas.loc[atlas["celltype_nmf"] == "SecB epithelium", "SecB_UCell"].values
y_true = np.concatenate([np.zeros(len(secA_vals)), np.ones(len(secB_vals))])
y_score = np.concatenate([secA_vals, secB_vals])
fpr, tpr, thresholds = roc_curve(y_true, y_score)
roc_auc = auc(fpr, tpr)
youden_idx = np.argmax(tpr - fpr)
roc_threshold = thresholds[youden_idx]

print(f"    ROC AUC (SecA vs SecB): {roc_auc:.4f}")
print(f"    Youden optimal threshold: {roc_threshold:.4f}")
print(f"    Sensitivity at threshold: {tpr[youden_idx]:.4f}")
print(f"    Specificity at threshold: {1 - fpr[youden_idx]:.4f}")

# Approach 2: Percentile-based thresholds from the full distribution
all_secB = atlas["SecB_UCell"].values
all_sorted = np.sort(all_secB)
p50_thresh = np.median(all_secB)
p75_thresh = np.percentile(all_secB, 75)

# Approach 3: NMF-label-informed boundaries
# SecB p25 = bottom of the NMF-defined SecB cells
secB_p25 = atlas_stats["SecB epithelium"]["p25"]
# Transitioning median = midpoint of the gradient
trans_median = atlas_stats["Transitioning epithelium"]["median"]
# SecA p75 = top of the NMF-defined SecA cells
secA_p75 = atlas_stats["SecA epithelium"]["p75"]

print(f"\n    Candidate thresholds:")
print(f"      Youden (ROC optimal):     {roc_threshold:.4f}")
print(f"      Atlas overall median:     {p50_thresh:.4f}")
print(f"      Atlas overall p75:        {p75_thresh:.4f}")
print(f"      SecA p75 (top of SecA):   {secA_p75:.4f}")
print(f"      Trans median:             {trans_median:.4f}")
print(f"      SecB p25 (bottom of SecB):{secB_p25:.4f}")

# Use Youden threshold as the primary classification boundary
PRIMARY_THRESH = roc_threshold
# Use SecB p25 as the "high confidence SecB" boundary
HIGH_CONF_THRESH = secB_p25

print(f"\n    → Primary threshold (Youden): {PRIMARY_THRESH:.4f}")
print(f"    → High-confidence SecB:       {HIGH_CONF_THRESH:.4f}")

# ============================================================================
# 5. CLASSIFY ORGANOID CELLS
# ============================================================================

print("\n[5] Classifying organoid cells on SecB_UCell...")

org["secB_class"] = "SecB-low"
org.loc[org["SecB_UCell"] >= PRIMARY_THRESH, "secB_class"] = "SecB-positive"
org.loc[org["SecB_UCell"] >= HIGH_CONF_THRESH, "secB_class"] = "SecB-high"

# Atlas percentile for each organoid cell
atlas_secB_sorted = np.sort(atlas["SecB_UCell"].values)
org["atlas_secB_pctile"] = org["SecB_UCell"].apply(
    lambda x: np.searchsorted(atlas_secB_sorted, x) / len(atlas_secB_sorted) * 100)

# Per-NMF-label percentiles (where does the organoid cell fall within each NMF group?)
for label in NMF_ORDER:
    short = label.replace(" epithelium", "").replace(" ", "")
    vals_sorted = np.sort(atlas.loc[atlas["celltype_nmf"] == label, "SecB_UCell"].values)
    org[f"pctile_in_{short}"] = org["SecB_UCell"].apply(
        lambda x, vs=vals_sorted: np.searchsorted(vs, x) / len(vs) * 100)

# Save classified organoid data
org.to_csv(os.path.join(OUT_DIR, "organoid_secB_classified_v5.csv"), index=False)

print(f"\n{'='*60}")
print(f"  04v5 complete!")
print(f"  Data:     {os.path.join(OUT_DIR, 'organoid_secB_classified_v5.csv')}")
print(f"{'='*60}")
