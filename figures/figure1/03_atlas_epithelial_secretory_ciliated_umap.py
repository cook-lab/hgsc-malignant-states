#!/usr/bin/env python3
"""
Figure 1F — Epithelial-only UMAP, Secretory vs Ciliated 2-class collapse
========================================================================
PURPOSE
    Epithelial-only atlas UMAP coloured by a 2-class collapse of schema_nmf:
      Secretory  <- {SecA, Intermediate, SecB}   (#E6A141)
      Ciliated   <- {Ciliated}                   (#E07850)
    Manuscript-ready (no title / panel letter / tick text). Ciliated drawn on
    top since it is the minority population.

INPUTS
    - fig_data_dir/meta.parquet  (UMAP1, UMAP2, schema_nmf; epithelial-only
      embedding; produced by figures/_prep/fig_secretory_polarization_00_prepare_data.py)

OUTPUTS
    - figures_dir/atlas_epithelial_secretory_ciliated_umap.{svg,png}

MANUSCRIPT PANEL(S): Fig 1F.

RUNTIME TIER: fast.

NOTE: secretory members updated "Transitioning" -> "Intermediate".
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

# ---------- Paths ----------
META_PQ = path("data_root", "2026_final_atlas", "output", "fig_secretory_polarization", "data", "meta.parquet")
OUT_SVG = path("figures_dir", "atlas_epithelial_secretory_ciliated_umap.svg")
OUT_PNG = path("figures_dir", "atlas_epithelial_secretory_ciliated_umap.png")

assert os.path.exists(META_PQ), f"Epithelial-meta parquet missing: {META_PQ}"

# ---------- Style ----------
SECRETORY_MEMBERS = {"SecA", "Intermediate", "SecB"}
CILIATED_MEMBERS = {"Ciliated"}
CLASS_PALETTE = {"Secretory": "#E6A141", "Ciliated": "#E07850"}
CLASS_ORDER = ["Secretory", "Ciliated"]

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
meta = pd.read_parquet(META_PQ, columns=["UMAP1", "UMAP2", "schema_nmf"])
meta = meta.dropna(subset=["UMAP1", "UMAP2", "schema_nmf"]).copy()


def _collapse(label: str) -> str:
    if label in CILIATED_MEMBERS:
        return "Ciliated"
    if label in SECRETORY_MEMBERS:
        return "Secretory"
    return "Other"


meta["class2"] = meta["schema_nmf"].astype(str).map(_collapse)
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
FIG_W_IN = 88 / 25.4
FIG_H_IN = 62 / 25.4
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
