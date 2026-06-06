#!/usr/bin/env python3
"""
Figure 3E — OPTO-098 SecB UCell across microenvironmental perturbations
======================================================================
PURPOSE
    Five perturbation violins (PBS control, IFNg, TGFb, TNFa, WNT7a; "Treatment"
    arm), with atlas Intermediate / SecB medians as dashed reference lines.
    IFNg is the only factor that induces SecB expression. Manuscript-ready.

INPUTS
    - OPTO-098 scored cells : <organoids_root>/output/08_OPTO98_SecB_Conditions/
        opto98_ucell_scores_aligned.csv  (experiment == "Treatment")
        EXTERNAL DEPENDENCY (override ORGANOIDS_ROOT).
    - atlas reference : output_root/18_ucell_atlas/{atlas_ucell_scores.csv,
        atlas_secretory_metadata.csv}  (subsampled 20,000 per NMF label, seed=SEED)

OUTPUTS
    - figures_dir/organoids_secB_perturbations.{png,svg}

MANUSCRIPT PANEL(S): Fig 3E.

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
OUT_PNG = path("figures_dir", "organoids_secB_perturbations.png")
OUT_SVG = path("figures_dir", "organoids_secB_perturbations.svg")

ATLAS_SAMPLE_N = 20000

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

ATLAS_PAL = {"SecA epithelium": "#E6A141", "Intermediate epithelium": "#C08E48",
             "SecB epithelium": "#6B4D2E"}
ATLAS_ORDER = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
TREAT_ORDER = ["PBS", "IFNg", "TGFb", "TNFa", "WNT7a"]
TREAT_PAL = {"PBS": "#B2CCAC", "IFNg": "#2E7D32", "TGFb": "#81C784",
             "TNFa": "#4CAF50", "WNT7a": "#A5D6A7"}

# ---------- Load ----------
opto = pd.read_csv(OPTO_CSV)
opto = opto[opto["experiment"] == "Treatment"].copy()

atlas_scores = pd.read_csv(ATLAS_SCORES)
atlas_meta = pd.read_csv(ATLAS_META, index_col=0)
atlas = atlas_meta.join(atlas_scores.set_index("barcode"), how="inner")
rng = np.random.default_rng(SEED)
atlas_medians = {}
for label in ATLAS_ORDER:
    vals = atlas.loc[atlas["celltype_nmf"] == label, "SecB_UCell"].values
    if len(vals) > ATLAS_SAMPLE_N:
        vals = vals[rng.choice(len(vals), ATLAS_SAMPLE_N, replace=False)]
    atlas_medians[label] = float(np.median(vals))

treat_data = [opto.loc[opto["sample_id"] == s, "SecB_UCell"].values for s in TREAT_ORDER]

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(76 / 25.4, 60 / 25.4))
n_treat = len(TREAT_ORDER)
x_treat = np.arange(n_treat) + 1


def _draw_violin(data, positions, fills, w=0.90, alpha=0.85):
    vp = ax.violinplot(data, positions=positions, showmeans=False,
                       showmedians=False, showextrema=False, widths=w)
    for body, c in zip(vp["bodies"], fills):
        body.set_facecolor(c); body.set_edgecolor("none")
        body.set_linewidth(0); body.set_alpha(alpha)
    ax.boxplot(data, positions=positions, widths=0.22, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="white", edgecolor="black", lw=0.5),
               whiskerprops=dict(color="black", lw=0.4),
               capprops=dict(color="black", lw=0.4),
               medianprops=dict(color="black", lw=0.8))


ax.axhline(atlas_medians["Intermediate epithelium"],
           color=ATLAS_PAL["Intermediate epithelium"], lw=0.6, ls=(0, (4, 2)),
           alpha=0.85, zorder=0)
ax.axhline(atlas_medians["SecB epithelium"], color=ATLAS_PAL["SecB epithelium"],
           lw=0.6, ls=(0, (4, 2)), alpha=0.85, zorder=0)
ax.text(n_treat + 0.55, atlas_medians["Intermediate epithelium"], "Atlas Int.",
        ha="right", va="bottom", fontsize=FN, color=ATLAS_PAL["Intermediate epithelium"])
ax.text(n_treat + 0.55, atlas_medians["SecB epithelium"], "Atlas SecB",
        ha="right", va="bottom", fontsize=FN, color=ATLAS_PAL["SecB epithelium"])

_draw_violin(treat_data, x_treat, [TREAT_PAL[s] for s in TREAT_ORDER])

ax.set_xticks(x_treat)
ax.set_xticklabels(list(TREAT_ORDER), rotation=30, ha="right")
ax.set_xlim(0.4, n_treat + 0.6)
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
