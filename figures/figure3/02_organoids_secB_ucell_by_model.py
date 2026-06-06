#!/usr/bin/env python3
"""
Figure 3B — Per-PDO-model SecB UCell distributions with atlas reference medians
==============================================================================
PURPOSE
    Per-model SecB UCell score violins (ordered low->high SecB median), with the
    atlas Intermediate / SecB medians overlaid as dashed reference lines.
    Manuscript-ready (no in-figure title / panel letter).

INPUTS
    - organoid scores : <organoids_root>/output/02_Secretory_Polarization_v5_SecB_only/
        organoid_secB_classified_v5.csv   (per-cell SecB_UCell x model)
        EXTERNAL DEPENDENCY (override ORGANOIDS_ROOT).
    - atlas reference : output_root/18_ucell_atlas/{atlas_ucell_scores.csv,
        atlas_secretory_metadata.csv}     (per-cell SecB_UCell x celltype_nmf)

OUTPUTS
    - figures_dir/organoids_secB_ucell_by_model.{png,svg}

MANUSCRIPT PANEL(S): Fig 3B.

RUNTIME TIER: moderate.

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
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

# ---------- Paths ----------
ORG_CSV = os.path.join(obj("organoids_root"),
                       "output/02_Secretory_Polarization_v5_SecB_only/"
                       "organoid_secB_classified_v5.csv")
ATLAS_SCORES = path("output_root", "18_ucell_atlas", "atlas_ucell_scores.csv")
ATLAS_META = path("output_root", "18_ucell_atlas", "atlas_secretory_metadata.csv")
OUT_PNG = path("figures_dir", "organoids_secB_ucell_by_model.png")
OUT_SVG = path("figures_dir", "organoids_secB_ucell_by_model.svg")

ATLAS_SAMPLE_N = 20000

# ---------- Style ----------
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
    "pdf.fonttype":     42,
    "svg.fonttype":     "none",
    "savefig.dpi":      450,
    "figure.dpi":       150,
})

ATLAS_PAL = {"SecA epithelium": "#E6A141",
             "Intermediate epithelium": "#C08E48",
             "SecB epithelium": "#6B4D2E"}
ATLAS_ORDER = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]

ASCITES_FILL = "#5665B6"
PRIMARY_FILL = "#8FBC8F"
ASCITES_MODELS = {"OCAD106", "OCAD93", "OCAD96", "OCAD97"}
PRIMARY_MODELS = {"OPTO98", "OPTO112", "OPTO129", "PDO66"}


def _model_fill(m):
    return ASCITES_FILL if m in ASCITES_MODELS else PRIMARY_FILL


# ---------- Load ----------
org = pd.read_csv(ORG_CSV)
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

model_order = (org.groupby("model")["SecB_UCell"].median()
               .sort_values(ascending=True).index.tolist())
pdo_data = [org.loc[org["model"] == m, "SecB_UCell"].values for m in model_order]

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(140 / 25.4, 60 / 25.4))
n_pdo = len(model_order)
x_pdo = np.arange(n_pdo) + 1


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
ax.text(n_pdo + 0.55, atlas_medians["Intermediate epithelium"], "Atlas Int.",
        ha="right", va="bottom", fontsize=FN, color=ATLAS_PAL["Intermediate epithelium"])
ax.text(n_pdo + 0.55, atlas_medians["SecB epithelium"], "Atlas SecB",
        ha="right", va="bottom", fontsize=FN, color=ATLAS_PAL["SecB epithelium"])

_draw_violin(pdo_data, x_pdo, [_model_fill(m) for m in model_order])

ax.set_xticks(x_pdo)
ax.set_xticklabels(list(model_order), rotation=30, ha="right")
ax.set_xlim(0.4, n_pdo + 0.6)
ax.set_ylabel("SecB UCell score")
ax.set_ylim(-0.02, 1.0)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
for spine in ("left", "bottom"):
    ax.spines[spine].set_linewidth(0.4)
    ax.spines[spine].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

ax.legend(handles=[Patch(facecolor=ASCITES_FILL, edgecolor="black", lw=0.4, label="PDO — ascites"),
                   Patch(facecolor=PRIMARY_FILL, edgecolor="black", lw=0.4, label="PDO — primary")],
          loc="upper right", frameon=False, fontsize=FN, handlelength=1.0,
          handletextpad=0.4, borderpad=0.2, borderaxespad=0.3, labelspacing=0.25)

fig.tight_layout(pad=0.3)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved:\n  {os.path.basename(OUT_PNG)}\n  {os.path.basename(OUT_SVG)}")
