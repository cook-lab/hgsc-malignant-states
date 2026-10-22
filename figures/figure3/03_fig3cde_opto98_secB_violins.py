#!/usr/bin/env python3
"""
Figure 3c-e — OPTO-098 SecB UCell across days / perturbations / media
=====================================================================
PURPOSE
    One wide figure with THREE sub-panels side by side sharing the y-axis
    ("SecB UCell score", labelled once on the left); in the code the panels are
    named D, E and F (also in the output name):
        c (D) = days in standard culture   (Day02, Day06, Day10, Day12)
        d (E) = perturbations              (PBS control leftmost, then WNT7A/TGFb/TNFa/IFNg)
        e (F) = media conditions           (the 8 media; PBS is shown in d)
    No sub-panel titles; "OPTO-098" labelled once (top-left). Shared atlas
    Intermediate + SecB dashed reference lines. Horizontal x labels.
    Light-blue violins with a white inner boxplot (black edge and median), as in
    panel b. The `mode` argument has the single value "violin".

INPUTS
    - <organoids_root>/output/08_OPTO98_SecB_Conditions/opto98_ucell_scores_aligned.csv
    - atlas reference : data_root/2026_final_atlas/output/18_ucell_atlas/
        {atlas_ucell_scores.csv, atlas_secretory_metadata.csv}

OUTPUTS
    - figures_dir/figure3/secB_DEF_combined_violin.{png,svg}

MANUSCRIPT PANEL(S): Fig 3c, 3d, 3e (one render). Footprint ~7.5 x 2.2 in. FONTS 8/6 pt.
Drawn at true footprint (no tight-bbox crop). RUNTIME TIER: fast.
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
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

OPTO_CSV = path("organoids_root",
                "output/08_OPTO98_SecB_Conditions/opto98_ucell_scores_aligned.csv")
ATLAS_SCORES = path("data_root", "2026_final_atlas", "output", "18_ucell_atlas", "atlas_ucell_scores.csv")
ATLAS_META = path("data_root", "2026_final_atlas", "output", "18_ucell_atlas", "atlas_secretory_metadata.csv")
OUT = {"violin": (path("figures_dir", "figure3", "secB_DEF_combined_violin.png"),
                  path("figures_dir", "figure3", "secB_DEF_combined_violin.svg")),
       }

FA, FK = 8, 6
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "axes.titlesize": 0,
    "xtick.labelsize": FK, "ytick.labelsize": FK, "legend.fontsize": FK,
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 450, "figure.dpi": 150,
})

VIOLIN_FILL = "#87CEFA"     # lighter blue (atlas T/NK)
BOX_EDGE = "#1F5C8B"
ATLAS_PAL = {"Intermediate epithelium": "#C08E48", "SecB epithelium": "#6B4D2E"}
YLIM = (-0.02, 0.90)
ATLAS_SAMPLE_N = 20000

# group orders (exactly as specified)
D_ORDER = [("Day02", "2"), ("Day06", "6"), ("Day10", "10"), ("Day12", "12")]
E_ORDER = [("PBS", "PBS"), ("WNT7a", "WNT7A"), ("TGFb", "TGFβ"), ("TNFa", "TNFα"), ("IFNg", "IFNγ")]
# media order: basal media first, then CGM, CGM + ascites, HPLM (the leading PBS entry is not drawn in e)
F_ORDER = [("PBS", "PBS"), ("BM", "BM"), ("BM1", "BM1"), ("BM2", "BM2"), ("BM02", "BM02"), ("BM03", "BM03"),
           ("CGM", "CGM"), ("CGM_Ascites", "CGM\nAsc"), ("HPLM", "HPLM")]

# ---------- Load ----------
full = pd.read_csv(OPTO_CSV)
growth = full[full["experiment"] == "Growth"]
treat = full[full["experiment"] == "Treatment"]
media = full[full["experiment"] == "Media"]


def _vals(df, sid):
    v = df.loc[df["sample_id"] == sid, "SecB_UCell"].values
    return v[np.isfinite(v)]


def build_panels(mode):
    """D = days, E = perturbations (PBS leftmost). F = media conditions only
    (no PBS: the E panel already carries PBS)."""
    D = [(lab, _vals(growth, sid)) for sid, lab in D_ORDER]
    E = [(lab, _vals(treat, sid)) for sid, lab in E_ORDER]
    media_groups = [(lab, _vals(media, sid)) for sid, lab in F_ORDER[1:]]
    F = media_groups
    return [("D", D), ("E", E), ("F", F)]

atlas_scores = pd.read_csv(ATLAS_SCORES)
atlas_meta = pd.read_csv(ATLAS_META, index_col=0)
atlas = atlas_meta.join(atlas_scores.set_index("barcode"), how="inner")
atlas["celltype_nmf"] = atlas["celltype_nmf"].replace(
    {"Transitioning epithelium": "Intermediate epithelium"})
rng0 = np.random.default_rng(SEED)
atlas_med = {}
for lab in ATLAS_PAL:
    v = atlas.loc[atlas["celltype_nmf"] == lab, "SecB_UCell"].values
    if len(v) > ATLAS_SAMPLE_N:
        v = v[rng0.choice(len(v), ATLAS_SAMPLE_N, replace=False)]
    atlas_med[lab] = float(np.median(v))


def _box(ax, data, positions, mode):
    if mode == "violin":   # matches panel b (02_fig3b_pdo_secB_by_model.py): white box, black edge/median
        ax.boxplot(data, positions=positions, widths=0.18, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="white", edgecolor="black", lw=0.5),
                   whiskerprops=dict(color="black", lw=0.4),
                   capprops=dict(color="black", lw=0.4),
                   medianprops=dict(color="black", lw=0.9), zorder=6)


def make_figure(mode):
    panels = build_panels(mode)
    n = [len(p[1]) for p in panels]
    fig = plt.figure(figsize=(7.5, 2.2), constrained_layout=False)
    gs = GridSpec(1, 3, width_ratios=n, wspace=0.06, left=0.075, right=0.985,
                  top=0.9, bottom=0.26)
    axes = []
    for j, (name, groups) in enumerate(panels):
        ax = fig.add_subplot(gs[0, j], sharey=axes[0] if axes else None)
        axes.append(ax)
        x = np.arange(len(groups)) + 1
        # shared reference lines
        for lab in ("Intermediate epithelium", "SecB epithelium"):
            ax.axhline(atlas_med[lab], color=ATLAS_PAL[lab], lw=0.6, ls=(0, (4, 2)),
                       alpha=0.85, zorder=1)
        data = [g[1] for g in groups]
        if mode == "violin":
            vp = ax.violinplot(data, positions=x, showmeans=False, showmedians=False,
                               showextrema=False, widths=0.85)
            for body in vp["bodies"]:   # matches panel b: black edge, alpha 0.85
                body.set_facecolor(VIOLIN_FILL); body.set_edgecolor("black")
                body.set_linewidth(0.3); body.set_alpha(0.85)
        _box(ax, data, x, mode)

        ax.set_xticks(x)
        ax.set_xticklabels([g[0] for g in groups], rotation=0, ha="center")
        ax.set_xlim(0.4, len(groups) + 0.6)
        ax.set_ylim(*YLIM)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.4); ax.spines[sp].set_color("black")
        ax.tick_params(width=0.4, length=2, colors="black")
        if j == 0:
            ax.set_ylabel("SecB UCell score")
        else:
            ax.tick_params(labelleft=False)
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", length=0)

    # atlas labels once, at the right of the last panel
    axF = axes[-1]
    nF = len(panels[-1][1])
    # white backing so the last violin's tail (HPLM) passes behind the labels; lifted
    # 1.5 pt so the backing stays clear of the dashed line
    LBL_BOX = dict(facecolor="white", edgecolor="none", pad=0.4)
    for lab, txt in (("SecB epithelium", "Atlas SecB"), ("Intermediate epithelium", "Atlas Int.")):
        axF.annotate(txt, xy=(nF + 0.55, atlas_med[lab]), xytext=(0, 1.5), textcoords="offset points",
                     ha="right", va="bottom", fontsize=FK, color=ATLAS_PAL[lab], bbox=LBL_BOX, zorder=8)
    # one shared OPTO-098 label, top-left of the figure
    fig.text(0.012, 0.96, "OPTO-098", ha="left", va="top", fontsize=FA,
             color=BOX_EDGE, fontweight="bold")

    png, svg = OUT[mode]
    fig.savefig(png, dpi=450)
    fig.savefig(svg)
    plt.close(fig)
    print(f"  {os.path.basename(png)}\n  {os.path.basename(svg)}")


print("Saved:")
make_figure("violin")
print("atlas:", {k: round(v, 3) for k, v in atlas_med.items()})
