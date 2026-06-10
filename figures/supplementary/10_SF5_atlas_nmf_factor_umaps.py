#!/usr/bin/env python3
"""
SF5 — NMF factor usage on the epithelial UMAP (10 panels)
=========================================================

Purpose
    10 epithelial-UMAP panels, one per NMF factor, showing per-cell factor
    usage (weight). Portrait letter page, 2 rows x 5 columns.

INPUTS
    output_root/fig_secretory_polarization/data/meta.parquet
        (epithelial UMAP coords; n=575,366; from atlas fig_secretory_polarization module)
    output_root/11d_epithelial_nmf/11d_nmf_usage.csv
        (atlas step 11d; 10 factors x 575,366 cells)

OUTPUTS
    output_root/figures/supplementary/SF5_atlas_nmf_factor_umaps.{svg,png}

MANUSCRIPT PANEL(S)
    SF5.

RUNTIME TIER
    moderate (10 rasterized scatters of 575k points).
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

META_PQ = path("data_root", "2026_final_atlas", "output", "fig_secretory_polarization", "data", "meta.parquet")
USAGE_CSV = path("data_root", "2026_final_atlas", "output", "11d_epithelial_nmf", "11d_nmf_usage.csv")
OUT_SVG = path("output_root", "figures", "supplementary", "SF5_atlas_nmf_factor_umaps.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF5_atlas_nmf_factor_umaps.png")

assert os.path.exists(META_PQ), f"Missing: {META_PQ}"
assert os.path.exists(USAGE_CSV), f"Missing: {USAGE_CSV}"

# ============================================================================
# FACTOR NAMES (display annotations of the 11d NMF top-loaded genes)
# ============================================================================

FACTOR_NAMES = {
    "Factor_1":  "Housekeeping\n(GAPDH, FTH1, ALDOA)",
    "Factor_2":  "SecB — adaptive\n(KRT7, LCN2, TACSTD2)",
    "Factor_3":  "Progenitor A — SecA\n(WT1, RCN2, PBX1)",
    "Factor_4":  "Interferon\n(IFI27, ISG15, IFI6)",
    "Factor_5":  "Cell cycle\n(CENPF, STMN1, UBE2C)",
    "Factor_6":  "AP-1 / stress\n(FOS, JUN, ATF3)",
    "Factor_7":  "Progenitor B\n(IGFBP2, RPL41, WFDC2)",
    "Factor_8":  "Secretory maturation\n(CLU, FTH1, WFDC2)",
    "Factor_9":  "EMT\n(VIM, LGALS1, ANXA2)",
    "Factor_10": "Metallothionein\n(MT1E, MT2A, MT1G)",
}
FACTOR_ORDER = [f"Factor_{i}" for i in range(1, 11)]

# ============================================================================
# STYLE
# ============================================================================

FA, FK, FN, FG = 6, 5.5, 5, 5.5

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
    "figure.dpi":       450,
    "savefig.dpi":      450,
})

CMAP = LinearSegmentedColormap.from_list(
    "beige_burgundy",
    [(0.00, "#DDD5CA"), (0.48, "#DDD5CA"), (0.65, "#C8878A"),
     (0.85, "#A04858"), (1.00, "#6B2838")],
    N=256,
)

# ============================================================================
# LOAD DATA
# ============================================================================

print("\nLoading epithelial metadata...", flush=True)
meta = pd.read_parquet(META_PQ, columns=["UMAP1", "UMAP2"])
print(f"  {len(meta):,} cells")

print("Loading NMF usage...", flush=True)
usage = pd.read_csv(USAGE_CSV, index_col=0)
print(f"  {usage.shape[0]:,} cells × {usage.shape[1]} factors")

common = meta.index.intersection(usage.index)
assert len(common) == len(meta), f"Index mismatch: {len(common)} vs {len(meta)}"
usage = usage.loc[meta.index]

x = meta["UMAP1"].to_numpy()
y = meta["UMAP2"].to_numpy()
xr = x.max() - x.min()
yr = y.max() - y.min()
pad_x = 0.04 * xr
pad_y = 0.04 * yr

SHUFFLE_ORDER = np.random.default_rng(seed=SEED).permutation(len(meta))
x_shuf = x[SHUFFLE_ORDER]
y_shuf = y[SHUFFLE_ORDER]

# ============================================================================
# PLOT — 2 rows × 5 cols on portrait letter page
# ============================================================================

print("Plotting...", flush=True)

NROWS, NCOLS = 2, 5
PAGE_W, PAGE_H = 7.5, 4.5

fig, axes = plt.subplots(NROWS, NCOLS, figsize=(PAGE_W, PAGE_H),
                         gridspec_kw={"wspace": 0.08, "hspace": 0.18})

for idx, factor in enumerate(FACTOR_ORDER):
    ax = axes[idx // NCOLS, idx % NCOLS]
    vals = usage[factor].to_numpy()
    vo = vals[SHUFFLE_ORDER]
    vmax = float(np.quantile(vals, 0.99))
    if vmax <= 0:
        vmax = 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)

    ax.scatter(x_shuf, y_shuf, c=vo, cmap=CMAP, norm=norm, s=0.02, marker="o",
               linewidths=0, rasterized=True, alpha=0.95)

    ax.set_xlim(x.min() - pad_x, x.max() + pad_x)
    ax.set_ylim(y.min() - pad_y, y.max() + pad_y)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_bounds(x.min() - pad_x, x.min() - pad_x + 0.18 * xr)
    ax.spines["left"].set_bounds(y.min() - pad_y, y.min() - pad_y + 0.18 * yr)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.spines["left"].set_linewidth(0.4)

    label_parts = FACTOR_NAMES[factor].split("\n")
    ax.text(0.03, 0.97, label_parts[0], transform=ax.transAxes, fontsize=FG,
            fontweight="bold", ha="left", va="top")
    if len(label_parts) > 1:
        ax.text(0.03, 0.88, label_parts[1], transform=ax.transAxes, fontsize=FG - 1,
                fontstyle="italic", ha="left", va="top", color="#444444")

    cax = ax.inset_axes([0.88, 0.05, 0.04, 0.3])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="vertical")
    cb.outline.set_linewidth(0.3)
    cb.ax.tick_params(length=1.5, width=0.3, pad=1, labelsize=FN - 1)
    cb.set_ticks([0.0, vmax])
    cb.set_ticklabels([f"{0.0:.2f}", f"{vmax:.2f}"])

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.06)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.06)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
