#!/usr/bin/env python3
"""
Figure 3F — Organoid cell-cycle stacked bar, SecA vs SecB
=========================================================
PURPOSE
    Two-bar stacked cell-cycle (G1/S/G2M) figure matching the atlas G2M panel:
    SecA = SecB-low (below p66.7 of SecB_UCell), SecB = SecB-high (>= p90); the
    mid band is dropped. Phases assigned by scanpy convention from S/G2M scores.

INPUTS
    - <organoids_root>/output/09_organoid_secB_characterization/
        09b_extended_per_cell.parquet (S_score, G2M_score) + metadata.csv (SecB_UCell)
        EXTERNAL DEPENDENCY (override ORGANOIDS_ROOT).

OUTPUTS
    - figures_dir/organoids_g2m_barplot_secAB.{svg,png}

MANUSCRIPT PANEL(S): Fig 3F.

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
from config.config import obj, path  # noqa: E402

# ---------- Paths ----------
DATA_DIR = os.path.join(obj("organoids_root"),
                        "output/09_organoid_secB_characterization")
OUT_SVG = path("figures_dir", "organoids_g2m_barplot_secAB.svg")
OUT_PNG = path("figures_dir", "organoids_g2m_barplot_secAB.png")

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

# ---------- Load + compute phases ----------
ext = pd.read_parquet(os.path.join(DATA_DIR, "09b_extended_per_cell.parquet"))
meta = pd.read_csv(os.path.join(DATA_DIR, "metadata.csv"))
meta.index = meta["barcode"]

secb_scores = meta["SecB_UCell"].values
p33 = np.nanpercentile(secb_scores, 66.7)
p67 = np.nanpercentile(secb_scores, 90)

secb_group = pd.Series(np.nan, index=meta.index, name="secb_group")
secb_group[secb_scores < p33] = "SecA"
secb_group[secb_scores >= p67] = "SecB"

keep = secb_group.dropna().index
secb_group = secb_group.loc[keep]
ext_keep = ext.loc[ext.index.intersection(keep)]

phase = pd.Series("G1", index=ext_keep.index)
phase[(ext_keep["S_score"] > 0) & (ext_keep["S_score"] > ext_keep["G2M_score"])] = "S"
phase[(ext_keep["G2M_score"] > 0) & (ext_keep["G2M_score"] > ext_keep["S_score"])] = "G2M"

df = pd.DataFrame({"group": secb_group.loc[ext_keep.index].values, "phase": phase.values})
ct = pd.crosstab(df["group"], df["phase"], normalize="index")
ct = ct.reindex(columns=["G1", "S", "G2M"]).reindex(CLUSTER_ORDER)

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(88 * 0.4 / 25.4, 55 / 25.4))
x = np.arange(len(CLUSTER_ORDER))
width = 0.55
g1_vals = ct["G1"].values
s_vals = ct["S"].values
g2m_vals = ct["G2M"].values
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
