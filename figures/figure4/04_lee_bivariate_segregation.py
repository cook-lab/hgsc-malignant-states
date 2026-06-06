#!/usr/bin/env python3
"""
Figure 4E — Lee's bivariate L: spatial segregation of SecA and SecB
===================================================================
PURPOSE
    Strip/jitter plot of Lee's L per whole-tissue sample (purple) and per TMA
    patient (beige). Negative L = spatial segregation between SecA and SecB
    UCell scores. Median lines per group; "<- segregated" annotation.

INPUTS
    - output_root/44_spatial_autocorrelation/interpretation_summary.csv  (WT L)
    - output_root/44_spatial_autocorrelation/tma_patient_level_lee.csv    (TMA L)

OUTPUTS
    - figures_dir/lee_bivariate_segregation.{png,svg}

MANUSCRIPT PANEL(S): Fig 4E.

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

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

# ---------- Paths ----------
INTERP = path("output_root", "44_spatial_autocorrelation", "interpretation_summary.csv")
TMA_PT = path("output_root", "44_spatial_autocorrelation", "tma_patient_level_lee.csv")
OUT_PNG = path("figures_dir", "lee_bivariate_segregation.png")
OUT_SVG = path("figures_dir", "lee_bivariate_segregation.svg")

# ---------- Style ----------
FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       FA,
    "axes.labelsize":  FA,
    "axes.titlesize":  FA,
    "xtick.labelsize": FK,
    "ytick.labelsize": FK,
    "legend.fontsize": FN,
    "pdf.fonttype":    42,
    "svg.fonttype":    "none",
    "savefig.dpi":     450,
    "figure.dpi":      150,
})

WT_COLOR = "#8A5DAF"
TMA_COLOR = "#D4A574"

# ---------- Load ----------
interp = pd.read_csv(INTERP)
wt_rows = interp[~interp["sample"].str.contains("tma")].copy()
wt_rows["lee_L"] = (wt_rows["SecA_SecB_segregation"]
                    .str.extract(r"L=([-\d.]+)")[0].astype(float))
wt_lee = wt_rows["lee_L"].values

tma_pt = pd.read_csv(TMA_PT)
tma_lee = tma_pt["mean_lee_L"].values
n_neg = (tma_lee < 0).sum()
n_sig = (tma_pt["n_sig"] > 0).sum()
print(f"  WT n={len(wt_lee)}; TMA n={len(tma_lee)}, neg={n_neg}, sig={n_sig}")

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(2.8, 1.65))
ax.axhline(0, color="black", lw=0.4, ls="-", alpha=0.3, zorder=0)
rng = np.random.default_rng(SEED)

tma_x = 2 + rng.uniform(-0.3, 0.3, size=len(tma_lee))
ax.scatter(tma_x, tma_lee, s=4, color=TMA_COLOR, alpha=0.5, edgecolors="none",
           zorder=2, rasterized=True)
wt_x = 1 + rng.uniform(-0.15, 0.15, size=len(wt_lee))
ax.scatter(wt_x, wt_lee, s=12, color=WT_COLOR, alpha=0.8, edgecolors="white",
           linewidths=0.3, zorder=3)

for pos, vals, col in [(1, wt_lee, WT_COLOR), (2, tma_lee, TMA_COLOR)]:
    med = np.median(vals)
    ax.plot([pos - 0.35, pos + 0.35], [med, med], color=col, lw=1.0, zorder=4)

ax.set_xticks([1, 2])
ax.set_xticklabels(["Whole\ntissue", "TMA\npatients"], fontsize=FK)
ax.set_xlim(0.4, 2.6)
ax.set_ylabel("Lee's L\n(SecA vs SecB)")
ax.text(0.02, 0.02, "← segregated", fontsize=3.5, color="grey", style="italic",
        transform=ax.transAxes, va="bottom")

for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
for spine in ("left", "bottom"):
    ax.spines[spine].set_linewidth(0.4)
    ax.spines[spine].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

fig.tight_layout(pad=0.3)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {os.path.basename(OUT_PNG)}\n  Saved: {os.path.basename(OUT_SVG)}\nDone.")
