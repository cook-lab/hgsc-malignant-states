#!/usr/bin/env python3
"""
Figure 4 (panel f) — GAM microenvironment vs polarization
=========================================================
PURPOSE
    Draw the pooled-GAM curves (distance to vasculature, matrix remodeling,
    cell density, hypoxia vs focal polarization) fitted by 09_fig4f_gam_fit.R,
    with matplotlib so the SVG is fully inlined and opens correctly in Affinity.

INPUTS  (written by 09_fig4f_gam_fit.R)
    - figures_dir/figure4/gam_curves.csv   (polarization, fitted, lower, upper, feature)
    - figures_dir/figure4/gam_labels.csv   (feature, polarization, fitted, fitted_orig, color)

OUTPUTS
    - figures_dir/figure4/gam_microenvironment_polarization.{png,svg}

MANUSCRIPT PANEL(S): Fig 4f. RUNTIME TIER: fast.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

CURVES = path("figures_dir", "figure4", "gam_curves.csv")
LABELS = path("figures_dir", "figure4", "gam_labels.csv")
OUT_PNG = path("figures_dir", "figure4", "gam_microenvironment_polarization.png")
OUT_SVG = path("figures_dir", "figure4", "gam_microenvironment_polarization.svg")

FEATURE_ORDER = ["Distance to vasculature", "Matrix remodeling", "Cell density", "Hypoxia"]
COLORS = {"Distance to vasculature": "#B87A7A", "Matrix remodeling": "#E07850",
          "Cell density": "#E6A141", "Hypoxia": "#D14E6C"}

FA, FK, FL = 7, 6, 8   # axis title, ticks, feature labels
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.dpi": 450, "figure.dpi": 150,
})

curves = pd.read_csv(CURVES)
labels = pd.read_csv(LABELS)

fig, ax = plt.subplots(figsize=(4.25, 2.5))
ax.axhline(0, ls=":", color="grey60" if False else "#999999", lw=0.4, zorder=0)

for feat in FEATURE_ORDER:
    d = curves[curves["feature"] == feat].sort_values("polarization")
    if d.empty:
        continue
    c = COLORS[feat]
    ax.fill_between(d["polarization"], d["lower"], d["upper"], color=c, alpha=0.12,
                    linewidth=0, zorder=1)
    ax.plot(d["polarization"], d["fitted"], color=c, lw=1.4, zorder=2, solid_capstyle="round")

# right-edge feature labels at the de-overlapped y (fitted), connector to fitted_orig
xmax = curves["polarization"].max()
for _, r in labels.iterrows():
    c = COLORS.get(r["feature"], "#555555")
    if abs(r["fitted"] - r["fitted_orig"]) > 0.02:
        ax.plot([xmax, xmax + 0.006], [r["fitted_orig"], r["fitted"]], color=c, lw=0.4, zorder=2)
    ax.text(xmax + 0.012, r["fitted"], r["feature"], color=c, fontsize=FL,
            fontweight="bold", va="center", ha="left", zorder=3)

ax.set_xlabel("Focal cell polarization score (SecA to SecB)", fontsize=FA)
ax.set_ylabel("50µm neighborhood\nz-score", fontsize=FA)
ax.set_xticks(np.arange(-0.75, 0.76, 0.25))
ax.tick_params(labelsize=FK, width=0.5, length=2.5, colors="black")
xr = xmax - curves["polarization"].min()
ax.set_xlim(curves["polarization"].min() - xr * 0.02, xmax + xr * 0.42)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.5); ax.spines[s].set_color("black")

fig.tight_layout(pad=0.3)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450, facecolor="white")
fig.savefig(OUT_SVG, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved (clean matplotlib SVG):\n  {OUT_PNG}\n  {OUT_SVG}")
