#!/usr/bin/env python3
"""
SF4C — CopyKAT aneuploid fraction by cell type
==============================================

Purpose
    Horizontal bar chart of the fraction of cells called aneuploid by CopyKAT,
    grouped by cell type. Epithelial (malignant) cells are highly aneuploid;
    non-malignant types are near-diploid.

INPUTS
    figures/supplementary/_data_cnv_aneuploid_by_celltype.csv
        (small committed aggregate of CopyKAT predictions across 251 samples;
         derived from output_root/19_cnv/per_sample/*; columns:
         celltype_pred, n_cells, n_aneuploid, frac_aneuploid)

OUTPUTS
    output_root/figures/supplementary/SF4C_atlas_cnv_aneuploid_by_celltype.{svg,png}

MANUSCRIPT PANEL(S)
    SF4C.

RUNTIME TIER
    fast.
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

DATA_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "_data_cnv_aneuploid_by_celltype.csv")
OUT_SVG = path("output_root", "figures", "supplementary", "SF4C_atlas_cnv_aneuploid_by_celltype.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF4C_atlas_cnv_aneuploid_by_celltype.png")

assert os.path.exists(DATA_CSV), f"Missing: {DATA_CSV}"

# ============================================================================
# STYLE
# ============================================================================

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
})

CELLTYPE_PALETTE = {
    "Epithelial": "#E6A141", "Mesothelial": "#D4A574", "Fibroblast": "#DDD5CA",
    "Smooth_Muscle": "#D14E6C", "Pericyte": "#B87A7A", "Endothelial": "#7D4E4E",
    "T_cell": "#87CEFA", "NK_cell": "#56AFC4", "B_cell": "#5665B6",
    "Plasma_cell": "#8A5DAF", "Macrophage": "#8FBC8F", "DC": "#2E8B57",
    "Neutrophil": "#6B8E23", "Mast": "#8B9B6B", "Erythrocyte": "#A0A0A0",
}
CELLTYPE_DISPLAY = {
    "Epithelial": "Epithelial", "Mesothelial": "Mesothelial", "Fibroblast": "Fibroblast",
    "Smooth_Muscle": "Smooth muscle", "Pericyte": "Pericyte", "Endothelial": "Endothelial",
    "T_cell": "T cell", "NK_cell": "NK cell", "B_cell": "B cell",
    "Plasma_cell": "Plasma cell", "Macrophage": "Macrophage", "DC": "Dendritic cell",
    "Neutrophil": "Neutrophil", "Mast": "Mast cell", "Erythrocyte": "Erythrocyte",
}
CELLTYPE_ORDER = [
    "Epithelial", "Mesothelial", "Fibroblast", "Smooth_Muscle", "Pericyte",
    "Endothelial", "T_cell", "NK_cell", "B_cell", "Plasma_cell",
    "Macrophage", "DC", "Neutrophil", "Mast", "Erythrocyte",
]

# ============================================================================
# DATA
# ============================================================================

print("Loading summary...", flush=True)
df = pd.read_csv(DATA_CSV).set_index("celltype_pred")
order = [ct for ct in CELLTYPE_ORDER if ct in df.index]
df = df.loc[order]
print(f"  {len(df)} cell types")

# ============================================================================
# PLOT
# ============================================================================

print("Plotting...", flush=True)

n = len(order)
fig, ax = plt.subplots(figsize=(7.5, 2.2))

y_pos = np.arange(n)
ax.barh(y_pos, df["frac_aneuploid"].values * 100, height=0.65,
        color=[CELLTYPE_PALETTE.get(ct, "#999999") for ct in order],
        edgecolor="white", linewidth=0.3)

for i, ct in enumerate(order):
    frac = df.loc[ct, "frac_aneuploid"] * 100
    n_cells = int(df.loc[ct, "n_cells"])
    if frac > 5:
        ax.text(frac - 1, i, f"n={n_cells:,}", ha="right", va="center",
                fontsize=FN - 0.5, color="white")
    else:
        ax.text(frac + 1, i, f"n={n_cells:,}", ha="left", va="center",
                fontsize=FN - 0.5, color="#555555")

ax.set_yticks(y_pos)
ax.set_yticklabels([CELLTYPE_DISPLAY.get(ct, ct) for ct in order], fontsize=FA)
ax.invert_yaxis()
for i, ct in enumerate(order):
    ax.get_yticklabels()[i].set_color(CELLTYPE_PALETTE.get(ct, "black"))
    ax.get_yticklabels()[i].set_fontweight("bold")

ax.set_xlabel("Aneuploid cells (%)", fontsize=FA)
ax.set_xlim(0, 100)
ax.set_xticks([0, 25, 50, 75, 100])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)
ax.tick_params(axis="both", which="both", length=2, width=0.4)
ax.axhline(0.5, color="#AAAAAA", linewidth=0.5, linestyle="--", zorder=0)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.04)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
