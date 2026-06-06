#!/usr/bin/env python3
"""
Figure 3D — OPTO-098 SecB UCell decay across days in standard culture
=====================================================================
PURPOSE
    Four violins (Day02 -> Day12, "Growth" arm) showing SecB UCell scores
    decrease over time in standard Matrigel culture, with atlas Intermediate /
    SecB medians as dashed reference lines. Manuscript-ready.

INPUTS
    - OPTO-098 scored cells : <organoids_root>/output/08_OPTO98_SecB_Conditions/
        opto98_ucell_scores_aligned.csv  (experiment == "Growth"; sample_id day)
        EXTERNAL DEPENDENCY (override ORGANOIDS_ROOT).
    - atlas reference : output_root/18_ucell_atlas/{atlas_ucell_scores.csv,
        atlas_secretory_metadata.csv}

OUTPUTS
    - figures_dir/organoids_secB_timecourse.{png,svg}

MANUSCRIPT PANEL(S): Fig 3D.

RUNTIME TIER: fast.

NOTE: atlas NMF label standardized "Transitioning epithelium" -> "Intermediate
epithelium".
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
from config.config import obj, path, SEED  # noqa: E402

# ---------- Paths ----------
OPTO_CSV = os.path.join(obj("organoids_root"),
                        "output/08_OPTO98_SecB_Conditions/opto98_ucell_scores_aligned.csv")
ATLAS_SCORES = path("output_root", "18_ucell_atlas", "atlas_ucell_scores.csv")
ATLAS_META = path("output_root", "18_ucell_atlas", "atlas_secretory_metadata.csv")
OUT_PNG = path("figures_dir", "organoids_secB_timecourse.png")
OUT_SVG = path("figures_dir", "organoids_secB_timecourse.svg")

# ---------- Style ----------
FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       FA,
    "axes.labelsize":  FA,
    "axes.titlesize":  0,
    "xtick.labelsize": FK,
    "ytick.labelsize": FK,
    "legend.fontsize": FN,
    "pdf.fonttype":    42,
    "svg.fonttype":    "none",
    "savefig.dpi":     450,
    "figure.dpi":      150,
})

DAY_ORDER = ["Day02", "Day06", "Day10", "Day12"]
DAY_PAL = {"Day02": "#2E5D2E", "Day06": "#4F8F4F", "Day10": "#8FBC8F", "Day12": "#C8DEC8"}
ATLAS_PAL = {"Intermediate epithelium": "#C08E48", "SecB epithelium": "#6B4D2E"}
ATLAS_SAMPLE_N = 20000

# ---------- Load ----------
opto = pd.read_csv(OPTO_CSV)
opto = opto[opto["experiment"] == "Growth"].copy()
day_data = [opto.loc[opto["sample_id"] == d, "SecB_UCell"].values for d in DAY_ORDER]

atlas_scores = pd.read_csv(ATLAS_SCORES)
atlas_meta = pd.read_csv(ATLAS_META, index_col=0)
atlas = atlas_meta.join(atlas_scores.set_index("barcode"), how="inner")
rng = np.random.default_rng(SEED)
atlas_medians = {}
for label in ATLAS_PAL.keys():
    vals = atlas.loc[atlas["celltype_nmf"] == label, "SecB_UCell"].values
    if len(vals) > ATLAS_SAMPLE_N:
        vals = vals[rng.choice(len(vals), ATLAS_SAMPLE_N, replace=False)]
    atlas_medians[label] = float(np.median(vals))

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(64 / 25.4, 60 / 25.4))
x_days = np.arange(len(DAY_ORDER)) + 1


def _draw_violin(data, positions, fills, w=0.90, alpha=0.85):
    vp = ax.violinplot(data, positions=positions, showmeans=False,
                       showmedians=False, showextrema=False, widths=w)
    for body, c in zip(vp["bodies"], fills):
        body.set_facecolor(c); body.set_edgecolor("none")
        body.set_linewidth(0); body.set_alpha(alpha)
    ax.boxplot(data, positions=positions, widths=0.18, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="white", edgecolor="black", lw=0.5),
               whiskerprops=dict(color="black", lw=0.4),
               capprops=dict(color="black", lw=0.4),
               medianprops=dict(color="black", lw=0.8))


ax.axhline(atlas_medians["Intermediate epithelium"],
           color=ATLAS_PAL["Intermediate epithelium"], lw=0.6, ls=(0, (4, 2)),
           alpha=0.85, zorder=0)
ax.axhline(atlas_medians["SecB epithelium"], color=ATLAS_PAL["SecB epithelium"],
           lw=0.6, ls=(0, (4, 2)), alpha=0.85, zorder=0)
ax.text(len(DAY_ORDER) + 0.55, atlas_medians["Intermediate epithelium"], "Atlas Int.",
        ha="right", va="bottom", fontsize=FN, color=ATLAS_PAL["Intermediate epithelium"])
ax.text(len(DAY_ORDER) + 0.55, atlas_medians["SecB epithelium"], "Atlas SecB",
        ha="right", va="bottom", fontsize=FN, color=ATLAS_PAL["SecB epithelium"])

_draw_violin(day_data, x_days, [DAY_PAL[d] for d in DAY_ORDER])

ax.set_xticks(x_days)
ax.set_xticklabels(list(DAY_ORDER), rotation=30, ha="right")
ax.set_xlim(0.4, x_days[-1] + 0.6)
ax.set_ylabel("SecB UCell score")
ax.set_ylim(-0.02, 1.0)
ax.text(0.02, 0.97, "OPTO-098", transform=ax.transAxes, ha="left", va="top",
        fontsize=FA, color="#2E5D2E", fontweight="bold")
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
for spine in ("left", "bottom"):
    ax.spines[spine].set_linewidth(0.4); ax.spines[spine].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

fig.tight_layout(pad=0.3)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved:\n  {os.path.basename(OUT_PNG)}\n  {os.path.basename(OUT_SVG)}")
