#!/usr/bin/env python3
"""
Fig 5b — per-state medians of secretory nuclear area (the medians quoted in the Results).

PURPOSE
    From the per-unit medians of 03_fig5b_nuclear_area.R, prints the median nuclear area of SecA,
    Intermediate and SecB cells over the 8 whole tissues and over the TMA patients (Results:
    39.5 vs 47.7 um^2 in the whole tissues), with the n and paired Wilcoxon SecA vs SecB p of each
    arm. The panel and its tests are drawn by 11_fig5_compose_panels.py; this script draws nothing.

INPUTS
    output_root/split_recompute/fig5b_morphometrics_presplit.csv (03_fig5b_nuclear_area.R)

OUTPUTS
    stdout only

MANUSCRIPT PANEL(S): Fig 5b (medians in the Results)
RUNTIME TIER: fast
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
from scipy.stats import wilcoxon, mannwhitneyu
R = os.environ.get("OUTPUT_ROOT", "output")
d = pd.read_csv(f"{R}/split_recompute/fig5b_morphometrics_presplit.csv")
ORD = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]

def wide(src, measure):
    return (d[d.source == src].pivot_table(index="unit", columns="cell_label_split",
            values=measure).reindex(columns=ORD).dropna())

# Nuclear area: the Fig 5b measure.
MEAS = [("nucleus_area", "Nuclear area (µm²)", "epithelial_nuclear_area")]
for meas, ylab, stem in MEAS:
    # --- WT: SecA vs SecB over the 8 whole tissues ---
    w = wide("WT", meas)
    p = wilcoxon(w[ORD[0]], w[ORD[2]]).pvalue
    # --- TMA patients: SecA vs SecB ---
    t = wide("TMA", meas); pair = [ORD[0], ORD[2]]
    p2 = wilcoxon(t[pair[0]], t[pair[1]]).pvalue
    print(f"{meas:14s} WT n={len(w)} p={p:.4g} | TMA n={len(t)} p={p2:.4g}")
    print("   WT medians:", {c.split()[0]: round(float(w[c].median()),3) for c in ORD})
    print("   TMA medians:", {c.split()[0]: round(float(t[c].median()),3) for c in ORD})
