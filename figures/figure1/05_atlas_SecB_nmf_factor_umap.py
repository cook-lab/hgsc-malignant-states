#!/usr/bin/env python3
"""
Figure 1H — SecB NMF Factor-2 usage on the epithelial UMAP
==========================================================
PURPOSE
    Single-panel epithelial UMAP coloured by continuous Factor_2 (SecB /
    adaptive) NMF usage, on the SecB brown colormap. Manuscript-ready.

INPUTS
    - fig_data_dir/meta.parquet   (UMAP1, UMAP2; epithelial cells)
    - 11d NMF usage : output_root/03_epithelial_nmf/11d_nmf_usage.csv (Factor_2)

OUTPUTS
    - figures_dir/atlas_SecB_nmf_factor_umap.{svg,png}

MANUSCRIPT PANEL(S): Fig 1H.

RUNTIME TIER: fast.
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
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import CFG, path, SEED  # noqa: E402

# ---------- Paths ----------
META_PQ = path("data_root", "2026_final_atlas", "output", "fig_secretory_polarization", "data", "meta.parquet")
USAGE_CSV = path("data_root", "2026_final_atlas", "output", "11d_epithelial_nmf", "11d_nmf_usage.csv")
OUT_SVG = path("figures_dir", "atlas_SecB_nmf_factor_umap.svg")
OUT_PNG = path("figures_dir", "atlas_SecB_nmf_factor_umap.png")

assert os.path.exists(META_PQ), f"Missing: {META_PQ}"
assert os.path.exists(USAGE_CSV), f"Missing: {USAGE_CSV}"

# ---------- Config ----------
FACTOR = CFG["polarization"]["factor"]   # Factor_2 (SecB-defining)
LABEL = "SecB — adaptive"
GENES = "(KRT7, LCN2, TACSTD2)"
CMAP = LinearSegmentedColormap.from_list(
    "beige_secb_brown",
    [(0.00, "#DDD5CA"), (0.40, "#DDD5CA"), (0.55, "#C0A880"),
     (0.70, "#9A7D55"), (0.85, "#6E5535"), (1.00, "#3D2B15")], N=256,
)

# ---------- Style ----------
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

# ---------- Load ----------
meta = pd.read_parquet(META_PQ, columns=["UMAP1", "UMAP2"])
usage = pd.read_csv(USAGE_CSV, index_col=0)
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

# ---------- Plot ----------
fig, ax = plt.subplots(1, 1, figsize=(3.2, 3.2))
vals = usage[FACTOR].to_numpy()
vo = vals[SHUFFLE_ORDER]
vmax = float(np.quantile(vals, 0.99)) or 1.0
norm = Normalize(vmin=0.0, vmax=vmax)

ax.scatter(x_shuf, y_shuf, c=vo, cmap=CMAP, norm=norm, s=0.08, marker="o",
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

ax.text(0.03, 0.97, LABEL, transform=ax.transAxes, fontsize=FG,
        fontweight="bold", ha="left", va="top")
ax.text(0.03, 0.88, GENES, transform=ax.transAxes, fontsize=FG - 1,
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
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}")
