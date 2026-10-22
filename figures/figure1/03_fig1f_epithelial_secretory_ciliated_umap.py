#!/usr/bin/env python3
"""
Figure 1f — Epithelial-only UMAP, Secretory vs Ciliated 2-class collapse
========================================================================
PURPOSE
    Epithelial-only atlas UMAP coloured by a 2-class collapse of the epithelial state:
      Secretory  <- {SecA, Intermediate, SecB}   (#E6A141)
      Ciliated   <- {Ciliated}                   (#E07850)
    Manuscript-ready (no title / panel letter / tick text). Ciliated drawn on
    top since it is the minority population.

INPUTS
    - figures/_prep/atlas_final_tables.py :: load_epithelial()  (575,366
      epithelial cells: epithelial-embedding UMAP1/UMAP2 + `state`
      SecA/Intermediate/SecB/Ciliated, metadata from obj("atlas_final");
      cache output_root/_prep_caches/atlas_final_epithelial.parquet).

OUTPUTS
    - figures_dir/figure1/atlas_epithelial_secretory_ciliated_umap.{svg,png}

MANUSCRIPT PANEL(S): Fig 1f.

RUNTIME TIER: fast.
"""
from __future__ import annotations

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
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

# ---------- Paths ----------
OUT_SVG = path("figures_dir", "figure1", "atlas_epithelial_secretory_ciliated_umap.svg")
OUT_PNG = path("figures_dir", "figure1", "atlas_epithelial_secretory_ciliated_umap.png")

# ---------- Style ----------
SECRETORY_MEMBERS = {"SecA", "Intermediate", "SecB"}
CILIATED_MEMBERS = {"Ciliated"}
CLASS_PALETTE = {"Secretory": "#E6A141", "Ciliated": "#E07850"}
CLASS_ORDER = ["Secretory", "Ciliated"]

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

# ---------- Load final-atlas epithelial table ----------
meta = load_epithelial(columns=["UMAP1", "UMAP2", "state"])
assert len(meta) == 575_366, len(meta)
meta = meta.dropna(subset=["UMAP1", "UMAP2", "state"]).copy()


def _collapse(label: str) -> str:
    if label in CILIATED_MEMBERS:
        return "Ciliated"
    if label in SECRETORY_MEMBERS:
        return "Secretory"
    return "Other"


meta["class2"] = meta["state"].astype(str).map(_collapse)
meta = meta[meta["class2"] != "Other"].copy()

counts = meta["class2"].value_counts().reindex(CLASS_ORDER, fill_value=0)
print("  per-class cell counts:")
for c in CLASS_ORDER:
    print(f"    {c:12s} {int(counts[c]):>10,}")

# ---------- Plot order ----------
rng = np.random.default_rng(seed=SEED)
sec_idx = np.where(meta["class2"].values == "Secretory")[0]
cil_idx = np.where(meta["class2"].values == "Ciliated")[0]
rng.shuffle(sec_idx)
order = np.concatenate([sec_idx, cil_idx])
meta = meta.iloc[order].reset_index(drop=True)

x = meta["UMAP1"].values
y = meta["UMAP2"].values
c = meta["class2"].astype(str).values

# ---------- Figure ----------
FIG_W_IN = 3.0
FIG_H_IN = 2.30
fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN), constrained_layout=False)
fig.subplots_adjust(left=0.08, right=0.66, top=0.97, bottom=0.10)

SCATTER_SIZE = 0.06
for cls in CLASS_ORDER:
    mask = c == cls
    if not mask.any():
        continue
    ax.scatter(x[mask], y[mask], c=CLASS_PALETTE[cls], s=SCATTER_SIZE, marker="o",
               linewidths=0, rasterized=True, alpha=0.9)

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

handles = [
    Line2D([0], [0], marker="o", linestyle="", markerfacecolor=CLASS_PALETTE[cls],
           markeredgecolor="none", markersize=4, label=cls)
    for cls in CLASS_ORDER
]
leg = ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
                frameon=False, fontsize=FN, handletextpad=0.4, labelspacing=0.45,
                borderaxespad=0.0)
for handle in leg.legend_handles:
    handle.set_markersize(4)

fig.savefig(OUT_SVG, format="svg", dpi=600, bbox_inches="tight", pad_inches=0.02)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print(f"[save] {OUT_SVG}\n[save] {OUT_PNG}\nDone.")
