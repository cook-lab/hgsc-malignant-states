#!/usr/bin/env python3
"""
Figure 2B — Cell-cycle phase stacked bar, SecA vs SecB
======================================================
PURPOSE
    Stacked bar of G1/S/G2M phase proportions for SecA and SecB epithelium
    (Intermediate and Ciliated excluded), with in-bar phase labels.

INPUTS
    - cell_cycle_phase_proportions.csv :
        output_root/04_functional/nmf_characterization/cell_cycle_phase_proportions.csv
        (built by the Fig-2B NMF-characterization prep helper; per-epitype G1/S/G2M %).

OUTPUTS
    - figures_dir/figure2/atlas_g2m_barplot_secAB.{svg,png}

MANUSCRIPT PANEL(S): Fig 2B.

RUNTIME TIER: fast.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
DATA_DIR = path("data_root", "202605_epitype_manuscript", "20260429_figures", "data", "nmf_characterization")
OUT_SVG = path("figures_dir", "figure2", "atlas_g2m_barplot_secAB.svg")
OUT_PNG = path("figures_dir", "figure2", "atlas_g2m_barplot_secAB.png")

# ---------- Style ----------
FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FA,
    "axes.labelsize":    FA,
    "xtick.labelsize":   FK,
    "ytick.labelsize":   FK,
    "legend.fontsize":   FN,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

CLUSTER_ORDER = ["SecA", "SecB"]
CLUSTER_PALETTE = {"SecA": "#E6A141", "SecB": "#9A7D55"}
BEIGE = "#DDD5CA"

# ---------- Load ----------
df = pd.read_csv(os.path.join(DATA_DIR, "cell_cycle_phase_proportions.csv"), index_col=0)
df = df.loc[CLUSTER_ORDER, ["G1", "S", "G2M"]]

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(88 * 0.4 / 25.4, 55 / 25.4))
x = np.arange(len(CLUSTER_ORDER))
width = 0.55
g1_vals = df["G1"].values
s_vals = df["S"].values
g2m_vals = df["G2M"].values
colors = [CLUSTER_PALETTE[c] for c in CLUSTER_ORDER]

ax.bar(x, g1_vals, width, color=BEIGE, edgecolor="white", linewidth=0.3, label="G1")
ax.bar(x, s_vals, width, bottom=g1_vals, color=colors, edgecolor="white",
       linewidth=0.3, alpha=0.55, label="S")
ax.bar(x, g2m_vals, width, bottom=g1_vals + s_vals, color=colors, edgecolor="white",
       linewidth=0.3, alpha=1.0, label="G2M")

ax.set_xticks(x)
ax.set_xticklabels(CLUSTER_ORDER, fontsize=FK)
ax.set_ylabel("Proportion", fontsize=FA)
ax.set_ylim(0, 1.0)
ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1.0, decimals=0))
ax.yaxis.set_major_locator(mticker.MultipleLocator(0.25))
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)
ax.tick_params(width=0.5, length=2)

for i in range(len(CLUSTER_ORDER)):
    ax.text(x[i], g1_vals[i] / 2, "G1", ha="center", va="center",
            fontsize=FN, color="black", fontweight="bold")
    ax.text(x[i], g1_vals[i] + s_vals[i] / 2, "S", ha="center", va="center",
            fontsize=FN, color="black", fontweight="bold")
    ax.text(x[i], g1_vals[i] + s_vals[i] + g2m_vals[i] / 2, "G2M", ha="center",
            va="center", fontsize=FN, color="white", fontweight="bold")

fig.tight_layout()
fig.savefig(OUT_SVG, format="svg")
fig.savefig(OUT_PNG, format="png")
plt.close(fig)
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}")
