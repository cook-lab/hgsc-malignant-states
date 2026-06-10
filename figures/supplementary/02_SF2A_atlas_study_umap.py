#!/usr/bin/env python3
"""
SF2A — Integrated atlas UMAP coloured by study (13 source studies)
==================================================================

Purpose
    Atlas integrated UMAP coloured by `study`. Sibling to the Fig 1A level-1
    UMAP: same data, same scatter geometry, different colour column.
    Manuscript-ready panel: no title, no panel letter, no axis-tick text.

INPUTS
    output_root/figures/data_fig1/meta.parquet
        (full-atlas UMAP cache from hgsc_atlas_final.h5ad; columns: UMAP1, UMAP2, study;
         produced by the Fig 1 data-extraction step)

OUTPUTS
    output_root/figures/supplementary/SF2A_atlas_study_umap.{svg,png}

MANUSCRIPT PANEL(S)
    SF2A.

RUNTIME TIER
    moderate (1.98M-point rasterized scatter).
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
from config.config import path, SEED  # noqa: E402

np.random.seed(SEED)

# ---------- Paths (central config) ----------
META_PQ = path("data_root", "202605_epitype_manuscript", "final_publication_figures", "data_fig1", "meta.parquet")
OUT_PNG = path("output_root", "figures", "supplementary", "SF2A_atlas_study_umap.png")
OUT_SVG = path("output_root", "figures", "supplementary", "SF2A_atlas_study_umap.svg")

assert os.path.exists(META_PQ), f"Cached meta parquet missing: {META_PQ}"

# ---------- Style ----------
# STUDY_DISPLAY ordered descending by total cell count for cross-figure consistency.
STUDY_DISPLAY = {
    "vazquez_garcia_2022": "Vazquez-Garcia 2022",
    "luo_2024":            "Luo 2024",
    "loret_2022":          "Loret 2022",
    "zheng_2023":          "Zheng 2023",
    "nath_2021":           "Nath 2021",
    "zhang_2022":          "Zhang 2022",
    "hornburg_2021":       "Hornburg 2021",
    "geistlinger_2020":    "Geistlinger 2020",
    "xu_2022":             "Xu 2022",
    "olbrecht_2021":       "Olbrecht 2021",
    "denisenko_2024":      "Denisenko 2024",
    "regner_2021":         "Regner 2021",
    "olalekan_2021":       "Olalekan 2021",
}

# Canonical study palette (KELLY_22 cyclically assigned to alphabetically-sorted
# study names; source: atlas 05_umap_suite auto_palette). Hard-coded so the figure
# renders without re-importing the upstream pipeline.
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
STUDY_ORDER = list(STUDY_DISPLAY.keys())

FA, FK, FN = 6, 5.5, 5

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

# ---------- Load cached meta ----------
print(f"[load] {META_PQ}")
meta = pd.read_parquet(META_PQ, columns=["UMAP1", "UMAP2", "study"])
print(f"  rows: {len(meta):,}")

meta = meta.dropna(subset=["study"]).copy()
present = set(meta["study"].astype(str).unique())
missing_in_palette = [c for c in present if c not in STUDY_PALETTE]
assert not missing_in_palette, f"Every study must have a palette entry: {missing_in_palette}"

# ---------- Plot order: shuffle so all studies stay visible ----------
rng = np.random.default_rng(seed=SEED)
meta = meta.iloc[rng.permutation(len(meta))].reset_index(drop=True)

x = meta["UMAP1"].values
y = meta["UMAP2"].values
s = meta["study"].astype(str).values

# ---------- Figure ----------
FIG_W_IN = 88 / 25.4
FIG_H_IN = 62 / 25.4

fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN), constrained_layout=False)
fig.subplots_adjust(left=0.08, right=0.66, top=0.97, bottom=0.10)

SCATTER_SIZE = 0.06
for study in STUDY_ORDER:
    mask = s == study
    if not mask.any():
        continue
    ax.scatter(x[mask], y[mask], c=STUDY_PALETTE[study], s=SCATTER_SIZE,
               marker="o", linewidths=0, rasterized=True, alpha=0.9)

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
           markerfacecolor=STUDY_PALETTE[study], markeredgecolor="none",
           markersize=4, label=STUDY_DISPLAY[study])
    for study in STUDY_ORDER
]
leg = ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
                frameon=False, fontsize=FN, handletextpad=0.4,
                labelspacing=0.45, borderaxespad=0.0)
for handle in leg.legend_handles:
    handle.set_markersize(4)

print(f"[save] {OUT_SVG}")
fig.savefig(OUT_SVG, format="svg", dpi=600, bbox_inches="tight", pad_inches=0.02)
print(f"[save] {OUT_PNG}")
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print("Done.")
