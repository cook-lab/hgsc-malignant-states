#!/usr/bin/env python3
"""
Figure 3H — Epitype x TCGA subtype dot matrix
=============================================
PURPOSE
    Dot matrix: TCGA molecular subtypes (x) by epitype (y). Dot size encodes
    mean composition %, dot colour encodes intensity; outer ring = Q75, dotted
    ring = median (across 96 ConsensusOV-classified pseudobulk samples).

INPUTS
    - 20 ConsensusOV tables:
        output_root/07_deconvolution_survival/20_consensusov/tables/
          {20c_subtype_composition_summary.csv, 20c_per_sample_joined.csv}

OUTPUTS
    - figures_dir/atlas_tcga_subtype_epitype.{svg,png}

MANUSCRIPT PANEL(S): Fig 3H.

RUNTIME TIER: fast.

NOTE: epitype label standardized "Transitioning" -> "Intermediate". Upstream
ConsensusOV cache column names (mean_pct_trans / pct_trans) are read as-is.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
TBL_DIR = path("data_root", "2026_final_atlas", "output", "20_consensusov", "tables")
SUMMARY_CSV = os.path.join(TBL_DIR, "20c_subtype_composition_summary.csv")
SAMPLE_CSV = os.path.join(TBL_DIR, "20c_per_sample_joined.csv")
OUT_SVG = path("figures_dir", "atlas_tcga_subtype_epitype.svg")
OUT_PNG = path("figures_dir", "atlas_tcga_subtype_epitype.png")

assert os.path.exists(SUMMARY_CSV), f"Missing: {SUMMARY_CSV}"
assert os.path.exists(SAMPLE_CSV), f"Missing: {SAMPLE_CSV}"

# ---------- Style ----------
FA, FK, FN = 7, 6.5, 6
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       FA,
    "axes.labelsize":  FA,
    "axes.titlesize":  0,
    "xtick.labelsize": FK,
    "ytick.labelsize": FK,
    "legend.fontsize": FN,
    "svg.fonttype":    "none",
    "pdf.fonttype":    42,
    "figure.dpi":      450,
    "savefig.dpi":     450,
})

EPI_PALETTE = {"SecA": "#E6A141", "Intermediate": "#CF8C2E",
               "SecB": "#9A7D55", "Ciliated": "#E07850"}
EPI_ORDER = ["Ciliated", "SecB", "Intermediate", "SecA"]
EPI_COL = {"SecA": "mean_pct_secA", "Intermediate": "mean_pct_trans",
           "SecB": "mean_pct_secB", "Ciliated": "mean_pct_ciliated"}
SUBTYPE_ORDER = ["DIF", "IMR", "MES", "PRO"]


def intensity_color(base_hex, frac):
    r, g, b = to_rgb(base_hex)
    t = max(0.15, frac)
    return (1 - t + t * r, 1 - t + t * g, 1 - t + t * b)


# ---------- Load ----------
df = pd.read_csv(SUMMARY_CSV).set_index("subtype").reindex(SUBTYPE_ORDER)

samp = pd.read_csv(SAMPLE_CSV).rename(columns={
    "pct_secA": "SecA", "pct_trans": "Intermediate",
    "pct_secB": "SecB", "pct_ciliated": "Ciliated",
})

iqr_data = {}
for st in SUBTYPE_ORDER:
    mask = samp["subtype_bulk"] == st
    for ep in EPI_ORDER:
        vals = samp.loc[mask, ep].values
        iqr_data[(st, ep)] = tuple(np.percentile(vals, [25, 50, 75]))

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(2.0, 2.0))
MAX_SIZE, MIN_SIZE = 350, 15
val_max = max(df.loc[st, EPI_COL[ep]] for ep in EPI_ORDER for st in SUBTYPE_ORDER)

for yi, ep in enumerate(EPI_ORDER):
    base_color = EPI_PALETTE[ep]
    for xi, st in enumerate(SUBTYPE_ORDER):
        val = df.loc[st, EPI_COL[ep]]
        frac = val / val_max
        q25, med, q75 = iqr_data[(st, ep)]
        size = MIN_SIZE + frac * (MAX_SIZE - MIN_SIZE)
        color = intensity_color(base_color, frac)
        size_q75 = MIN_SIZE + (q75 / val_max) * (MAX_SIZE - MIN_SIZE)
        ax.scatter(xi, yi, s=size_q75, c="none", edgecolors=base_color,
                   linewidths=0.4, alpha=0.35, zorder=2)
        ax.scatter(xi, yi, s=size, c=[color], edgecolors=base_color, linewidths=0.3, zorder=3)
        size_med = MIN_SIZE + (med / val_max) * (MAX_SIZE - MIN_SIZE)
        inv = ax.transData.inverted()
        p0 = ax.transData.transform((xi, yi))
        p1 = (p0[0] + np.sqrt(size_med / np.pi), p0[1])
        r_data_x = inv.transform(p1)[0] - xi
        ax.add_patch(plt.Circle((xi, yi), r_data_x, fill=False, edgecolor=base_color,
                                linewidth=0.5, linestyle=(0, (1.5, 2)), alpha=0.8, zorder=4))

ax.set_xticks(range(len(SUBTYPE_ORDER)))
ax.set_xticklabels([f"{st}\n(n={int(df.loc[st, 'n_samples'])})" for st in SUBTYPE_ORDER],
                   fontsize=FN, linespacing=0.85)
ax.set_yticks(range(len(EPI_ORDER)))
ax.set_yticklabels(EPI_ORDER, fontsize=FN, fontweight="bold")
for yi, ep in enumerate(EPI_ORDER):
    ax.get_yticklabels()[yi].set_color(EPI_PALETTE[ep])

for yi in range(len(EPI_ORDER)):
    ax.axhline(yi, color="#EEEEEE", linewidth=0.3, zorder=0)
for xi in range(len(SUBTYPE_ORDER)):
    ax.axvline(xi, color="#EEEEEE", linewidth=0.3, zorder=0)

ax.set_xlim(-0.55, len(SUBTYPE_ORDER) - 0.45)
ax.set_ylim(-0.55, len(EPI_ORDER) - 0.45)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(axis="both", length=0, pad=3)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.04)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}")
