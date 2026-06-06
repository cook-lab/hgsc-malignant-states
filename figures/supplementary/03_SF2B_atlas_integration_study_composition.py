#!/usr/bin/env python3
"""
SF2B — Study composition per cell type (stacked barplot)
========================================================

Purpose
    Stacked barplot of the proportional contribution of each study to each
    cell type. Good integration = each cell type draws from multiple studies
    rather than being dominated by one.

INPUTS
    output_root/figures/data/atlas_final_umap.parquet
        (columns used: study, celltype_pred; from 00b_extract_integration_umaps.py)

OUTPUTS
    output_root/figures/supplementary/SF2B_atlas_integration_study_composition.{svg,png}

MANUSCRIPT PANEL(S)
    SF2B.

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
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

DATA_PQ = path("output_root", "figures", "data", "atlas_final_umap.parquet")
OUT_SVG = path("output_root", "figures", "supplementary", "SF2B_atlas_integration_study_composition.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF2B_atlas_integration_study_composition.png")

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

STUDY_PALETTE = {
    "denisenko_2024":      "#875692",
    "geistlinger_2020":    "#F38400",
    "hornburg_2021":       "#A1CAF1",
    "loret_2022":          "#BE0032",
    "luo_2024":            "#C2B280",
    "nath_2021":           "#848482",
    "olalekan_2021":       "#008856",
    "olbrecht_2021":       "#E68FAC",
    "regner_2021":         "#0067A5",
    "vazquez_garcia_2022": "#F99379",
    "xu_2022":             "#604E97",
    "zhang_2022":          "#F6A600",
    "zheng_2023":          "#B3446C",
}
STUDY_DISPLAY = {
    "denisenko_2024": "Denisenko 2024", "geistlinger_2020": "Geistlinger 2020",
    "hornburg_2021": "Hornburg 2021", "loret_2022": "Loret 2022", "luo_2024": "Luo 2024",
    "nath_2021": "Nath 2021", "olalekan_2021": "Olalekan 2021",
    "olbrecht_2021": "Olbrecht 2021", "regner_2021": "Regner 2021",
    "vazquez_garcia_2022": "Vazquez-Garcia 2022", "xu_2022": "Xu 2022",
    "zhang_2022": "Zhang 2022", "zheng_2023": "Zheng 2023",
}
STUDY_ORDER = list(STUDY_PALETTE.keys())

CELLTYPE_ORDER = [
    "Epithelial", "Macrophage", "T_cell", "Fibroblast", "Endothelial",
    "B_cell", "Plasma_cell", "NK_cell", "Smooth_Muscle", "Mesothelial",
    "Pericyte", "DC", "Mast", "Neutrophil", "Erythrocyte",
]
CELLTYPE_DISPLAY = {
    "Epithelial": "Epithelial", "Macrophage": "Macrophage", "T_cell": "T cell",
    "Fibroblast": "Fibroblast", "Endothelial": "Endothelial", "B_cell": "B cell",
    "Plasma_cell": "Plasma cell", "NK_cell": "NK cell", "Smooth_Muscle": "Smooth muscle",
    "Mesothelial": "Mesothelial", "Pericyte": "Pericyte", "DC": "Dendritic cell",
    "Mast": "Mast cell", "Neutrophil": "Neutrophil", "Erythrocyte": "Erythrocyte",
}

# ============================================================================
# DATA
# ============================================================================

print("Loading data...", flush=True)
df = pd.read_parquet(DATA_PQ, columns=["study", "celltype_pred"])
print(f"  {len(df):,} cells")

ct = pd.crosstab(df["celltype_pred"], df["study"], normalize="index")
ct_types = [c for c in CELLTYPE_ORDER if c in ct.index]
ct = ct.loc[ct_types]
ct_studies = [s for s in STUDY_ORDER if s in ct.columns]
ct = ct[ct_studies]

# ============================================================================
# PLOT
# ============================================================================

print("Plotting...", flush=True)

fig, ax = plt.subplots(figsize=(88 / 25.4, 80 / 25.4))

y_pos = np.arange(len(ct_types))
left = np.zeros(len(ct_types))
for study in ct_studies:
    widths = ct[study].values
    ax.barh(y_pos, widths, left=left, height=0.7, color=STUDY_PALETTE[study], linewidth=0)
    left += widths

ax.set_yticks(y_pos)
ax.set_yticklabels([CELLTYPE_DISPLAY.get(c, c) for c in ct_types], fontsize=FK)
ax.set_xlabel("Proportion of cells", fontsize=FA)
ax.set_xlim(0, 1)
ax.invert_yaxis()
ax.spines[["top", "right"]].set_visible(False)

handles = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor=STUDY_PALETTE[s],
           markersize=4, linewidth=0, label=STUDY_DISPLAY[s])
    for s in ct_studies
]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12),
          ncol=4, fontsize=FN - 0.5, frameon=False, handletextpad=0.3, columnspacing=0.8)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight")
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_SVG}")
print(f"Saved: {OUT_PNG}")
