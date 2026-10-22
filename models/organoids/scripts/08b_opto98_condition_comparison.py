#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202603_organoids/scripts/08b_opto98_condition_comparison.py, the script that ran.
# Kept: opto98_statistical_tests.csv and the printed Kruskal-Wallis test, which back the Fig 3c-e text. Removed: the
# atlas-threshold classification, the condition summaries, the figures and the HTML report. The untrimmed script,
# which the docstring below describes, is in the lab project folder.

"""
08b — OPTO-098 SecB UCell by Condition: Atlas Comparison
========================================================
Cook Lab · 2026

Applies the v5 SecB-only classification (atlas-grounded Youden + SecB-p25
thresholds) to the OPTO-098 dataset across Growth, Media, and Treatment
conditions. Asks: which conditions shift the SecB signal relative to
untreated baseline and atlas reference?

Input:
  - output/08_OPTO98_SecB_Conditions/opto98_ucell_scores_aligned.csv
  - output/02_Secretory_Polarization_v5_SecB_only/organoid_secB_classified_v5.csv  (baseline OPTO98)
  - /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_ucell_scores.csv
  - /Volumes/CookLab/Sarah/2026_final_atlas/output/18_ucell_atlas/atlas_secretory_metadata.csv

Output:
  - output/08_OPTO98_SecB_Conditions/08b_opto98_condition_report.html
  - output/08_OPTO98_SecB_Conditions/opto98_condition_summary.csv
  - output/08_OPTO98_SecB_Conditions/opto98_classified.csv
  - figs/08_OPTO98_SecB_Conditions/*.svg + *.pdf
"""

import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
PROJECT   = "/Volumes/CookLab/Sarah/2026_organoids"
OUT_DIR   = os.path.join(PROJECT, "output", "08_OPTO98_SecB_Conditions")
os.makedirs(OUT_DIR, exist_ok=True)

# Ordered sample display
GROWTH_ORDER = ["Day02", "Day06", "Day10", "Day12"]

# ============================================================================
# 1. LOAD DATA
# ============================================================================

print("=" * 60)
print("  08b — OPTO-098 SecB by Condition")
print("=" * 60)

print("\n[1] Loading data...")
opto = pd.read_csv(os.path.join(OUT_DIR, "opto98_ucell_scores_aligned.csv"))
print(f"    OPTO-098 cells (filtered): {len(opto):,}")

print(f"\n    Experiments: {sorted(opto['experiment'].unique())}")
print(f"    Samples:     {sorted(opto['sample_id'].unique())}")

# ============================================================================
# 5. STATISTICAL TESTS
# ============================================================================

print("\n[5] Statistical tests...")

stat_results = []

# Within Treatment arm: each treatment vs PBS
pbs_vals = opto.loc[opto["sample_id"] == "PBS", "SecB_UCell"].values
for treat in ["IFNg", "TGFb", "TNFa", "WNT7a"]:
    t_vals = opto.loc[opto["sample_id"] == treat, "SecB_UCell"].values
    stat, pval = mannwhitneyu(t_vals, pbs_vals, alternative="two-sided")
    direction = "up" if np.median(t_vals) > np.median(pbs_vals) else "down"
    stat_results.append({
        "comparison": f"{treat} vs PBS",
        "arm": "Treatment",
        "n1": len(t_vals), "n2": len(pbs_vals),
        "median_1": np.median(t_vals),
        "median_2": np.median(pbs_vals),
        "direction": direction,
        "U_stat": stat, "p_value": pval,
    })
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    print(f"    {treat} vs PBS: median {np.median(t_vals):.4f} vs {np.median(pbs_vals):.4f} "
          f"({direction}) p={pval:.2e} {sig}")

# Within Growth arm: Kruskal-Wallis across timepoints
growth_groups = [opto.loc[opto["sample_id"] == d, "SecB_UCell"].values for d in GROWTH_ORDER]
kw_stat, kw_p = kruskal(*growth_groups)
print(f"\n    Growth time course (Kruskal-Wallis): H={kw_stat:.2f}, p={kw_p:.2e}")

# Pairwise: each timepoint vs Day02 (earliest)
day02_vals = opto.loc[opto["sample_id"] == "Day02", "SecB_UCell"].values
for day in ["Day06", "Day10", "Day12"]:
    d_vals = opto.loc[opto["sample_id"] == day, "SecB_UCell"].values
    stat, pval = mannwhitneyu(d_vals, day02_vals, alternative="two-sided")
    direction = "up" if np.median(d_vals) > np.median(day02_vals) else "down"
    stat_results.append({
        "comparison": f"{day} vs Day02",
        "arm": "Growth",
        "n1": len(d_vals), "n2": len(day02_vals),
        "median_1": np.median(d_vals),
        "median_2": np.median(day02_vals),
        "direction": direction,
        "U_stat": stat, "p_value": pval,
    })
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    print(f"    {day} vs Day02: median {np.median(d_vals):.4f} vs {np.median(day02_vals):.4f} "
          f"({direction}) p={pval:.2e} {sig}")

# Within Media arm: each vs BM (baseline media)
bm_vals = opto.loc[opto["sample_id"] == "BM", "SecB_UCell"].values
for med in ["BM1", "BM02", "BM03", "BM2", "CGM", "CGM_Ascites", "HPLM"]:
    m_vals = opto.loc[opto["sample_id"] == med, "SecB_UCell"].values
    if len(m_vals) == 0:
        continue
    stat, pval = mannwhitneyu(m_vals, bm_vals, alternative="two-sided")
    direction = "up" if np.median(m_vals) > np.median(bm_vals) else "down"
    stat_results.append({
        "comparison": f"{med} vs BM",
        "arm": "Media",
        "n1": len(m_vals), "n2": len(bm_vals),
        "median_1": np.median(m_vals),
        "median_2": np.median(bm_vals),
        "direction": direction,
        "U_stat": stat, "p_value": pval,
    })
    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    print(f"    {med} vs BM: median {np.median(m_vals):.4f} vs {np.median(bm_vals):.4f} "
          f"({direction}) p={pval:.2e} {sig}")

stats_df = pd.DataFrame(stat_results)
stats_df.to_csv(os.path.join(OUT_DIR, "opto98_statistical_tests.csv"), index=False)

print(f"\n{'='*60}")
print(f"  08b complete!")
print(f"  Stats:    {os.path.join(OUT_DIR, 'opto98_statistical_tests.csv')}")
print(f"{'='*60}")
