#!/usr/bin/env python3
"""
SF4B — Canonical markers x cell types dotplot
=============================================

Purpose
    Dot plot of 26 canonical marker genes across 13 Level-1 cell types.
    Dot size = fraction expressing; dot colour = normalised mean expression.

INPUTS
    output_root/09f_dotplot_canonical_markers/dotplot_stats.csv
        (atlas step 09f; columns: celltype, gene, frac_expressing, mean_expression,
         n_cells, norm_expression)
    output_root/09f_dotplot_canonical_markers/lilra4_stats.csv
        (LILRA4 / pDC marker stats appended to the dotplot; produced alongside 09f)

OUTPUTS
    output_root/figures/supplementary/SF4B_atlas_canonical_markers_dotplot.{svg,png}

MANUSCRIPT PANEL(S)
    SF4B.

RUNTIME TIER
    fast (renders from pre-computed stats CSVs).
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

DOTPLOT_CSV = path("output_root", "09f_dotplot_canonical_markers", "dotplot_stats.csv")
LILRA4_CSV = path("output_root", "09f_dotplot_canonical_markers", "lilra4_stats.csv")
OUT_SVG = path("output_root", "figures", "supplementary", "SF4B_atlas_canonical_markers_dotplot.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF4B_atlas_canonical_markers_dotplot.png")

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

# Spatial expression gradient: light grey -> burgundy/maroon
LAPAZ_CMAP = LinearSegmentedColormap.from_list(
    "lapaz_r",
    [(0.00, "#F5F5F5"), (0.20, "#F0D0B0"), (0.40, "#D4916B"),
     (0.60, "#B05A5A"), (0.80, "#7A2040"), (1.00, "#4A0C2A")],
)

# ============================================================================
# DATA
# ============================================================================

print("Loading dotplot stats...", flush=True)
df = pd.read_csv(DOTPLOT_CSV)
lilra4 = pd.read_csv(LILRA4_CSV)
df = pd.concat([df, lilra4], ignore_index=True)
print(f"  {len(df)} rows, {df['celltype'].nunique()} celltypes, {df['gene'].nunique()} genes")

CELLTYPE_ORDER = [
    "Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte",
    "Endothelial", "T/NK cell", "B cell", "Plasma cell",
    "Macrophage", "DC", "Neutrophil", "Mast cell",
]
CELLTYPE_DISPLAY = {
    "Epithelial": "Epithelial", "Mesothelial": "Mesothelial", "Fibroblast": "Fibroblast",
    "Smooth muscle": "Smooth muscle", "Pericyte": "Pericyte", "Endothelial": "Endothelial",
    "T/NK cell": "T/NK cell", "B cell": "B cell", "Plasma cell": "Plasma cell",
    "Macrophage": "Macrophage", "DC": "Dendritic cell", "Neutrophil": "Neutrophil",
    "Mast cell": "Mast cell",
}
GENE_ORDER = [
    "EPCAM", "PAX8", "MSLN", "UPK3B", "COL1A1", "DCN", "MYH11", "ACTA2",
    "PDGFRB", "RGS5", "CDH5", "VWF", "CD3E", "NKG7", "MS4A1", "CD79A",
    "MZB1", "XBP1", "CD14", "C1QA", "CD1C", "LILRA4", "CSF3R", "S100A8",
    "KIT", "TPSAB1",
]

ct_present = [ct for ct in CELLTYPE_ORDER if ct in df["celltype"].values]
gene_present = [g for g in GENE_ORDER if g in df["gene"].values]

CELLTYPE_PALETTE = {
    "Epithelial": "#E6A141", "Mesothelial": "#D4A574", "Fibroblast": "#DDD5CA",
    "Smooth muscle": "#D14E6C", "Pericyte": "#B87A7A", "Endothelial": "#7D4E4E",
    "T/NK cell": "#87CEFA", "B cell": "#5665B6", "Plasma cell": "#8A5DAF",
    "Macrophage": "#8FBC8F", "DC": "#2E8B57", "Neutrophil": "#6B8E23",
    "Mast cell": "#8B9B6B",
}

GENE_GROUP_BOUNDARIES = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26]

# ============================================================================
# PLOT
# ============================================================================

print("Plotting...", flush=True)

n_ct = len(ct_present)
n_gene = len(gene_present)

fig, ax = plt.subplots(figsize=(7.5, 2.2))
DOT_SCALE = 120

for i, ct in enumerate(ct_present):
    for j, gene in enumerate(gene_present):
        row = df[(df["celltype"] == ct) & (df["gene"] == gene)]
        if row.empty:
            continue
        frac = row["frac_expressing"].values[0]
        expr = row["norm_expression"].values[0]
        ax.scatter(j, n_ct - 1 - i, s=frac * DOT_SCALE, c=[expr], cmap=LAPAZ_CMAP,
                   vmin=0, vmax=1, edgecolors="black", linewidths=0.3, zorder=3)

for boundary in GENE_GROUP_BOUNDARIES[:-1]:
    ax.axvline(boundary - 0.5, color="#DDDDDD", linewidth=0.4, zorder=1)

# Epithelial | Stromal | Immune lymphoid | Immune myeloid
CT_GROUP_BOUNDARIES = [2, 6, 9]
for boundary in CT_GROUP_BOUNDARIES:
    ax.axhline(n_ct - boundary - 0.5, color="#DDDDDD", linewidth=0.4, zorder=1)

ax.set_xticks(range(n_gene))
ax.set_xticklabels(gene_present, rotation=45, ha="right", fontsize=FA, style="italic")
ax.set_yticks(range(n_ct))
ax.set_yticklabels([CELLTYPE_DISPLAY.get(ct, ct) for ct in reversed(ct_present)], fontsize=FA)

for i, ct in enumerate(reversed(ct_present)):
    color = CELLTYPE_PALETTE.get(ct, "black")
    ax.get_yticklabels()[i].set_color(color)
    ax.get_yticklabels()[i].set_fontweight("bold")

for i, ct in enumerate(reversed(ct_present)):
    color = CELLTYPE_PALETTE.get(ct, "#999999")
    ax.add_patch(plt.Rectangle((-1.1, i - 0.3), 0.5, 0.6, facecolor=color,
                               edgecolor="none", clip_on=False, zorder=5))

ax.set_xlim(-0.5, n_gene - 0.5)
ax.set_ylim(-0.5, n_ct - 0.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)
ax.tick_params(axis="both", which="both", length=2, width=0.4)

size_vals = [0.2, 0.4, 0.6, 0.8]
size_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#888888",
           markeredgecolor="black", markeredgewidth=0.3,
           markersize=np.sqrt(v * DOT_SCALE) * 0.55, linewidth=0, label=f"{int(v*100)}%")
    for v in size_vals
]
leg_size = ax.legend(handles=size_handles, title="% expressing", loc="upper left",
                     bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=FN,
                     title_fontsize=FA, handletextpad=0.4, labelspacing=0.8, borderpad=0.2)
ax.add_artist(leg_size)

sm = plt.cm.ScalarMappable(cmap=LAPAZ_CMAP, norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.3, pad=0.12, aspect=10, anchor=(0.0, 0.0))
cbar.ax.tick_params(labelsize=FN)
cbar.set_label("Norm. expression", fontsize=FN)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.06)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.06)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
