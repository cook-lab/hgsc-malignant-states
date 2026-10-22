#!/usr/bin/env python3
"""
Figure 1a — Atlas integrated UMAP coloured by celltype_level1 (13 lineages)
===========================================================================
PURPOSE
    Manuscript-ready whole-atlas UMAP (no in-figure title / panel letter / tick
    text): scatter + UMAP1/UMAP2 axis stubs + categorical legend. Reads the
    final-atlas obs table via the shared _prep loader.

INPUTS
    - figures/_prep/atlas_final_tables.py :: load_obs()  (UMAP1, UMAP2,
      celltype_level1 for all 1,980,703 cells, read from obj("atlas_final");
      cache output_root/_prep_caches/atlas_final_obs.parquet).

OUTPUTS
    - figures_dir/figure1/atlas_level1_umap.{svg,png}

MANUSCRIPT PANEL(S): Fig 1a.

RUNTIME TIER: moderate (renders ~2M rasterized points).
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
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402
from figures._prep.atlas_final_tables import load_obs  # noqa: E402

# ---------- Paths ----------
OUT_SVG = path("figures_dir", "figure1", "atlas_level1_umap.svg")
OUT_PNG = path("figures_dir", "figure1", "atlas_level1_umap.png")

# ---------- Style ----------
LEVEL1_ORDER = [
    "Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte",
    "Endothelial", "T/NK cell", "B cell", "Plasma cell", "Macrophage",
    "DC", "Neutrophil", "Mast cell",
]
LEVEL1_PALETTE = {
    "Epithelial":     "#E6A141",
    "Mesothelial":    "#D4A574",
    "Fibroblast":     "#DDD5CA",
    "Smooth muscle":  "#D14E6C",
    "Pericyte":       "#B87A7A",
    "Endothelial":    "#7D4E4E",
    "T/NK cell":      "#87CEFA",
    "B cell":         "#5665B6",
    "Plasma cell":    "#8A5DAF",
    "Macrophage":     "#8FBC8F",
    "DC":             "#2E8B57",
    "Neutrophil":     "#6B8E23",
    "Mast cell":      "#8B9B6B",
}

FA, FK, FN = 8, 6, 8

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        FA,
    "axes.labelsize":   FA,
    "axes.titlesize":   0,
    "xtick.labelsize":  FK,
    "ytick.labelsize":  FK,
    "legend.fontsize":  FN,
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
})

# ---------- Load final-atlas obs ----------
print("[load] load_obs() (final atlas)")
meta = load_obs(columns=["UMAP1", "UMAP2", "celltype_level1"])
print(f"  rows: {len(meta):,}")
assert len(meta) == 1_980_703, len(meta)

meta = meta.dropna(subset=["celltype_level1"]).copy()

counts = meta["celltype_level1"].value_counts().reindex(LEVEL1_ORDER, fill_value=0)
print("  per-celltype cell counts:")
for ct in LEVEL1_ORDER:
    print(f"    {ct:18s} {counts[ct]:>10,}")

# ---------- Plot order: shuffle so all categories are visible ----------
rng = np.random.default_rng(seed=SEED)
shuffle = rng.permutation(len(meta))
meta = meta.iloc[shuffle].reset_index(drop=True)

x = meta["UMAP1"].values
y = meta["UMAP2"].values
ct = meta["celltype_level1"].astype(str).values

# ---------- Figure ----------
FIG_W_IN = 3.15
FIG_H_IN = 2.40

fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN), constrained_layout=False)
fig.subplots_adjust(left=0.08, right=0.66, top=0.97, bottom=0.10)

SCATTER_SIZE = 0.06
for cat in LEVEL1_ORDER:
    mask = ct == cat
    if not mask.any():
        continue
    ax.scatter(
        x[mask], y[mask], c=LEVEL1_PALETTE[cat], s=SCATTER_SIZE, marker="o",
        linewidths=0, rasterized=True, alpha=0.9,
    )

ax.set_xlabel("UMAP1", fontsize=FA, labelpad=2)
ax.set_ylabel("UMAP2", fontsize=FA, labelpad=2)
ax.set_xticks([]); ax.set_yticks([])
for spine_name in ("top", "right"):
    ax.spines[spine_name].set_visible(False)
xr = x.max() - x.min()
yr = y.max() - y.min()
pad_x = 0.04 * xr
pad_y = 0.04 * yr
ax.set_xlim(x.min() - pad_x, x.max() + pad_x)
ax.set_ylim(y.min() - pad_y, y.max() + pad_y)
ax.spines["bottom"].set_bounds(x.min() - pad_x, x.min() - pad_x + 0.18 * xr)
ax.spines["left"].set_bounds(y.min() - pad_y, y.min() - pad_y + 0.18 * yr)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["left"].set_linewidth(0.5)
ax.set_aspect("equal")

# ---------- Legend ----------
handles = [
    Line2D([0], [0], marker="o", linestyle="",
           markerfacecolor=LEVEL1_PALETTE[c], markeredgecolor="none",
           markersize=4, label=c)
    for c in LEVEL1_ORDER
]
leg = ax.legend(
    handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
    frameon=False, fontsize=FN, handletextpad=0.4, labelspacing=0.45,
    borderaxespad=0.0,
)
for handle in leg.legend_handles:
    handle.set_markersize(4)

# ---------- Save ----------
fig.savefig(OUT_SVG, format="svg", dpi=600, bbox_inches="tight", pad_inches=0.02)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print(f"[save] {OUT_SVG}\n[save] {OUT_PNG}\nDone.")
