#!/usr/bin/env python3
"""
SF4A — Cell annotation validation (scANVI prediction probability)
=================================================================

Purpose
    Violin plot of the scANVI posterior probability per Level-1 cell type, with
    an overall "% confident (prob > 0.5)" annotation.

INPUTS
    obj("atlas_final")  (hgsc_atlas_final.h5ad; obs: celltype_pred,
        celltype_probability, celltype_confident)

OUTPUTS
    output_root/figures/supplementary/SF4A_atlas_annotation_validation.{svg,png}

MANUSCRIPT PANEL(S)
    SF4A.

RUNTIME TIER
    fast (subsamples to 20k cells per cell type for violins).
"""

import os
import sys

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

ATLAS_H5AD = obj("atlas_final")
OUT_SVG = path("output_root", "figures", "supplementary", "SF4A_atlas_annotation_validation.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF4A_atlas_annotation_validation.png")

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
CELLTYPE_ORDER = [
    "Epithelial", "Mesothelial", "Fibroblast", "Smooth_Muscle", "Pericyte",
    "Endothelial", "T_cell", "NK_cell", "B_cell", "Plasma_cell",
    "Macrophage", "DC", "Neutrophil", "Mast",
]
CELLTYPE_DISPLAY = {
    "Epithelial": "Epithelial", "Mesothelial": "Mesothelial", "Fibroblast": "Fibroblast",
    "Smooth_Muscle": "Smooth muscle", "Pericyte": "Pericyte", "Endothelial": "Endothelial",
    "T_cell": "T cell", "NK_cell": "NK cell", "B_cell": "B cell",
    "Plasma_cell": "Plasma cell", "Macrophage": "Macrophage", "DC": "Dendritic cell",
    "Neutrophil": "Neutrophil", "Mast": "Mast cell",
}

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading atlas h5ad (backed)...", flush=True)
adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
print(f"  {adata.n_obs:,} cells")

obs = adata.obs[["celltype_pred", "celltype_probability", "celltype_confident"]].copy()
del adata

MAX_PER_GROUP = 20_000
rng = np.random.default_rng(SEED)
frames = []
for ct in CELLTYPE_ORDER:
    sub = obs[obs["celltype_pred"] == ct]
    if len(sub) > MAX_PER_GROUP:
        idx = rng.choice(len(sub), size=MAX_PER_GROUP, replace=False)
        sub = sub.iloc[idx]
    frames.append(sub)
obs_sub = pd.concat(frames, axis=0)
print(f"  Subsampled to {len(obs_sub):,} cells for plotting")

# ============================================================================
# PLOT
# ============================================================================

print("Plotting...", flush=True)

order_present = [ct for ct in CELLTYPE_ORDER if ct in obs_sub["celltype_pred"].values]
n_groups = len(order_present)

fig, ax = plt.subplots(figsize=(180 / 25.4, 60 / 25.4))

data_list = [obs_sub.loc[obs_sub["celltype_pred"] == ct, "celltype_probability"].values
             for ct in order_present]

parts = ax.violinplot(data_list, positions=range(n_groups), showmeans=False,
                      showmedians=True, showextrema=False, widths=0.7)
for i, body in enumerate(parts["bodies"]):
    body.set_facecolor(CELLTYPE_PALETTE[order_present[i]])
    body.set_edgecolor("black")
    body.set_linewidth(0.3)
    body.set_alpha(0.85)
parts["cmedians"].set_colors("black")
parts["cmedians"].set_linewidths(0.6)

for i, ct in enumerate(order_present):
    n_total = int((obs["celltype_pred"] == ct).sum())
    ax.text(i, 1.04, f"n={n_total:,}", ha="center", va="bottom",
            fontsize=FN - 0.5, color="#555555", rotation=45)

n_confident = int(obs["celltype_confident"].sum())
n_total = len(obs)
frac = n_confident / n_total * 100
ax.text(0.99, 0.02, f"{frac:.1f}% confident (prob > 0.5)", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=FN, color="#555555")

ax.set_xticks(range(n_groups))
ax.set_xticklabels([CELLTYPE_DISPLAY.get(ct, ct) for ct in order_present],
                   rotation=45, ha="right", fontsize=FK)
ax.set_ylabel("scANVI prediction probability", fontsize=FA)
ax.set_ylim(0, 1.15)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)
ax.tick_params(axis="both", which="both", length=2, width=0.4)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.04)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
