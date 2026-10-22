#!/usr/bin/env python3
"""
Figure 3b — per-model SecB UCell violins
========================================
PURPOSE
    One violin per PDO model (boxplot inside) of the per-cell SecB UCell score,
    coloured by tissue origin (primary = atlas T/NK light blue, ascites = atlas
    B-cell dark blue), with the atlas Intermediate and SecB epithelium median
    UCell as dashed reference lines. Model names are HORIZONTAL on x.
    Manuscript-ready (no title / panel letter).

    NOTE ON READABILITY: the per-cell SecB score is zero-inflated and the
    inflation differs sharply by model (ascites/OCAD ~9-35% exact zeros -> real
    violins; primary/OPTO ~62-86% zeros, median 0 -> the body collapses onto the
    floor). The boxplot inside makes each model's median/IQR readable even where
    the violin bulges at 0.

    MODEL-LEVEL TEST: ascites- vs primary-derived PDOs compared on each
    model's median per-cell SecB UCell score, two-sided exact Mann-Whitney U test with the model as
    the unit (4 vs 4). The medians separate completely (ascites 0.073-0.199, primary 0-0.035), so
    U = 16 and p = 2/70 = 0.029, the smallest p a 4 vs 4 test can give. Per-cell p values are not
    used (cells are not independent units).

INPUTS
    - organoid scores : <organoids_root>/output/02_Secretory_Polarization_v5_SecB_only/
        organoid_secB_classified_v5.csv   (per-cell SecB_UCell x model)
    - atlas reference : data_root/2026_final_atlas/output/18_ucell_atlas/
        {atlas_ucell_scores.csv, atlas_secretory_metadata.csv}

OUTPUTS
    - figures_dir/figure3/secB_by_model_violin_letter.{png,svg}
    - figures_dir/figure3/secB_by_model_violin_letter_models.csv (per model: origin, cells, median, mean, % zero)
    - figures_dir/figure3/secB_by_model_violin_letter_stats.csv  (the model-level ascites vs primary test)

MANUSCRIPT PANEL(S): Fig 3b. Footprint 4.25 x 2.0 in (half a letter page). Legend sits in one
row above the axes. FONTS 8 pt / 6 pt ticks.
Drawn at true footprint (no tight-bbox crop). RUNTIME TIER: moderate.
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
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

ORG_CSV = path("organoids_root",
               "output/02_Secretory_Polarization_v5_SecB_only/organoid_secB_classified_v5.csv")
ATLAS_SCORES = path("data_root", "2026_final_atlas", "output", "18_ucell_atlas", "atlas_ucell_scores.csv")
ATLAS_META = path("data_root", "2026_final_atlas", "output", "18_ucell_atlas", "atlas_secretory_metadata.csv")
OUT_PNG = path("figures_dir", "figure3", "secB_by_model_violin_letter.png")
OUT_SVG = path("figures_dir", "figure3", "secB_by_model_violin_letter.svg")
OUT_MODELS = path("figures_dir", "figure3", "secB_by_model_violin_letter_models.csv")
OUT_STATS = path("figures_dir", "figure3", "secB_by_model_violin_letter_stats.csv")

ATLAS_SAMPLE_N = 20000

FA, FK, FN = 8, 6, 6
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "axes.titlesize": 0,
    "xtick.labelsize": FK, "ytick.labelsize": FK, "legend.fontsize": FN,
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 450, "figure.dpi": 150,
})

# atlas-lineage recolour: primary = T/NK light blue, ascites = B-cell dark blue
ASCITES_FILL = "#5665B6"
PRIMARY_FILL = "#87CEFA"
ASCITES_MODELS = {"OCAD106", "OCAD93", "OCAD96", "OCAD97"}
ATLAS_PAL = {"SecA epithelium": "#E6A141", "Intermediate epithelium": "#C08E48",
             "SecB epithelium": "#6B4D2E"}
ATLAS_ORDER = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
ATLAS_SHORT = {"SecA epithelium": "Atlas SecA", "Intermediate epithelium": "Atlas Int.",
               "SecB epithelium": "Atlas SecB"}


def _model_fill(m):
    return ASCITES_FILL if m in ASCITES_MODELS else PRIMARY_FILL


# ---------- Load ----------
org = pd.read_csv(ORG_CSV)
atlas_scores = pd.read_csv(ATLAS_SCORES)
atlas_meta = pd.read_csv(ATLAS_META, index_col=0)
atlas = atlas_meta.join(atlas_scores.set_index("barcode"), how="inner")
atlas["celltype_nmf"] = atlas["celltype_nmf"].replace(
    {"Transitioning epithelium": "Intermediate epithelium"})
rng = np.random.default_rng(SEED)
atlas_medians = {}
for lab in ATLAS_ORDER:
    vals = atlas.loc[atlas["celltype_nmf"] == lab, "SecB_UCell"].values
    if len(vals) > ATLAS_SAMPLE_N:
        vals = vals[rng.choice(len(vals), ATLAS_SAMPLE_N, replace=False)]
    atlas_medians[lab] = float(np.median(vals))

stat = org.groupby("model")["SecB_UCell"].agg(["median", "mean"])
model_order = stat.sort_values(["median", "mean"]).index.tolist()
pdo_data = [org.loc[org["model"] == m, "SecB_UCell"].values for m in model_order]

# ---------- Model-level test: ascites vs primary PDOs (unit = model) ----------
models = org.groupby("model")["SecB_UCell"].agg(
    n_cells="size", median_SecB_UCell="median", mean_SecB_UCell="mean",
    pct_zero=lambda v: 100 * float((v == 0).mean())).loc[model_order]
models.insert(0, "origin", ["ascites" if m in ASCITES_MODELS else "primary" for m in models.index])
asc = models.loc[models.origin == "ascites", "median_SecB_UCell"]
pri = models.loc[models.origin == "primary", "median_SecB_UCell"]
assert len(asc) == 4 and len(pri) == 4, (len(asc), len(pri))
mw = mannwhitneyu(asc, pri, alternative="two-sided", method="exact")
test = pd.DataFrame([{
    "panel": "b", "comparison": "ascites vs primary PDO models", "unit": "model (median per-cell SecB UCell)",
    "test": "Mann-Whitney U, two-sided, exact", "n_ascites": len(asc), "n_primary": len(pri),
    "U": float(mw.statistic), "p": float(mw.pvalue),
    "ascites_median_range": f"{asc.min():.3f}-{asc.max():.3f}", "primary_median_range": f"{pri.min():.3f}-{pri.max():.3f}"}])
os.makedirs(os.path.dirname(OUT_STATS), exist_ok=True)
models.round(4).reset_index().to_csv(OUT_MODELS, index=False)
test.to_csv(OUT_STATS, index=False)

# ---------- Plot ----------
FIG_W, FIG_H = 4.25, 2.0   # half of an 8.5 in letter page; 2 in tall (author, 2026-10-08)
fig = plt.figure(figsize=(FIG_W, FIG_H), constrained_layout=True)
ax = fig.add_subplot(111)
n = len(model_order)
xpos = np.arange(n) + 1


def _draw_violin(data, positions, fills, w=0.9, alpha=0.85):
    vp = ax.violinplot(data, positions=positions, showmeans=False,
                       showmedians=False, showextrema=False, widths=w)
    for body, c in zip(vp["bodies"], fills):
        body.set_facecolor(c); body.set_edgecolor("black")
        body.set_linewidth(0.3); body.set_alpha(alpha)
    ax.boxplot(data, positions=positions, widths=0.2, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="white", edgecolor="black", lw=0.5),
               whiskerprops=dict(color="black", lw=0.4),
               capprops=dict(color="black", lw=0.4),
               medianprops=dict(color="black", lw=0.9))


# only Intermediate + SecB reference lines (SecA removed to focus on Int->SecB)
for lab in ("Intermediate epithelium", "SecB epithelium"):
    ax.axhline(atlas_medians[lab], color=ATLAS_PAL[lab], lw=0.6, ls=(0, (4, 2)),
               alpha=0.85, zorder=0)
    ax.text(0.5, atlas_medians[lab], ATLAS_SHORT[lab], ha="left", va="bottom",   # left end: clear of OCAD96
            fontsize=FN, color=ATLAS_PAL[lab])

_draw_violin(pdo_data, xpos, [_model_fill(m) for m in model_order])

ax.set_xticks(xpos)
ax.set_xticklabels(list(model_order), rotation=0, ha="center")
ax.set_xlim(0.45, n + 0.55)
ax.set_ylabel("SecB UCell score")
ax.set_ylim(-0.02, 0.78)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
for spine in ("left", "bottom"):
    ax.spines[spine].set_linewidth(0.4); ax.spines[spine].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

ax.legend(handles=[Patch(facecolor=PRIMARY_FILL, edgecolor="black", lw=0.4, label="PDO — primary"),
                   Patch(facecolor=ASCITES_FILL, edgecolor="black", lw=0.4, label="PDO — ascites")],
          loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, fontsize=FN,
          handlelength=1.0, handletextpad=0.4, borderpad=0.1, columnspacing=1.2)

fig.savefig(OUT_PNG, dpi=450)
fig.savefig(OUT_SVG)
plt.close(fig)
print(f"\nSaved:\n  {os.path.basename(OUT_PNG)}\n  {os.path.basename(OUT_SVG)}")
print("order:", model_order)
print("atlas:", {k: round(v, 3) for k, v in atlas_medians.items()})
print(f"ascites vs primary (per-model medians, 4 vs 4): U = {mw.statistic:.0f}, exact two-sided p = {mw.pvalue:.4f}")
