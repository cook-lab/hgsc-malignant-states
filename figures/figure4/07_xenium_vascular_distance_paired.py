#!/usr/bin/env python3
"""
Figure 4H — Distance-to-vasculature across the secretory gradient (paired)
==========================================================================
PURPOSE
    Paired-sample comparison of median distance-to-nearest-vascular-cell across
    the SecA -> Intermediate -> SecB gradient, for the 8 whole-tissue samples
    (3-position paired traces) and the TMA cohort (SecA->SecB with violins).
    Wilcoxon paired SecA-vs-SecB bracket per panel.

INPUTS
    - output_root/22_vascular_proximity/vascular_distance_summary.csv
      (per-sample/per-patient median um to nearest Pericyte/Endothelial;
       tissue == "whole_tissue" | "TMA")

OUTPUTS
    - figures_dir/figure4/xenium_vascular_distance_paired.{png,svg}

MANUSCRIPT PANEL(S): Fig 4H.

RUNTIME TIER: fast.

NOTE: epitype label standardized "Transitioning" -> "Intermediate".
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
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
SUMMARY = path("data_root", "2026_final_xenium_analysis", "output", "22_vascular_proximity", "vascular_distance_summary.csv")
OUT_PNG = path("figures_dir", "figure4", "xenium_vascular_distance_paired.png")
OUT_SVG = path("figures_dir", "figure4", "xenium_vascular_distance_paired.svg")

# ---------- Style ----------
FA, FK, FN = 8, 7, 6.5
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

EPI_ORDER = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
EPI_LABEL = ["SecA", "Int.", "SecB"]
PAL = {"SecA epithelium": "#E6A141", "Intermediate epithelium": "#C08E48",
       "SecB epithelium": "#9A7D55"}

# ---------- Load + pivot ----------
d = pd.read_csv(SUMMARY)
d["cell_label"] = d["cell_label"].replace({"Transitioning epithelium": "Intermediate epithelium"})


def pivot_wide(sub):
    return (sub.pivot_table(index="sample_id", columns="cell_label", values="median_dist")
              .reindex(columns=EPI_ORDER).dropna())


wt_wide = pivot_wide(d[d["tissue"] == "whole_tissue"])
tma_wide = pivot_wide(d[d["tissue"] == "TMA"])
print(f"  WT: n={len(wt_wide)}; TMA: n={len(tma_wide)}")


def paired_p(df, a, b):
    sub = df[[a, b]].dropna()
    if len(sub) < 5:
        return np.nan
    return wilcoxon(sub[a], sub[b]).pvalue


def fmt_p(p):
    if np.isnan(p):
        return "NA"
    if p < 0.0001:
        return "p < 0.0001"
    if p < 0.001:
        return f"p = {p:.4f}"
    return f"p = {p:.3f}"


wt_p_AB = paired_p(wt_wide, EPI_ORDER[0], EPI_ORDER[2])
tma_p_AB = paired_p(tma_wide, EPI_ORDER[0], EPI_ORDER[2])

# ---------- Plot ----------
fig, axes = plt.subplots(1, 2, figsize=(65 / 25.4, 35 / 25.4), sharey=True,
                         gridspec_kw={"width_ratios": [1.1, 1.0]})
TMA_PAIR = [EPI_ORDER[0], EPI_ORDER[2]]
TMA_LABEL = ["SecA", "SecB"]


def _draw_paired_traces(ax, wide, x_pos, cols, line_alpha, line_lw, dot_size,
                        dot_alpha=0.95, dot_edge=0.3):
    vals = wide[cols].values
    for row in vals:
        ax.plot(x_pos, row, color="grey", lw=line_lw, alpha=line_alpha, zorder=2)
    for j, ct in enumerate(cols):
        ax.scatter(np.full(len(vals), x_pos[j]), vals[:, j], s=dot_size, color=PAL[ct],
                   edgecolor="white", linewidth=dot_edge, alpha=dot_alpha, zorder=3)


ax = axes[0]
x_pos_wt = np.arange(3) + 1
_draw_paired_traces(ax, wt_wide, x_pos_wt, EPI_ORDER, line_alpha=0.55, line_lw=0.7, dot_size=14)
ax.set_xticks(x_pos_wt); ax.set_xticklabels(EPI_LABEL); ax.set_xlim(0.55, 3.45)
ax.set_ylabel("Distance to nearest\nvascular cell (µm)")
ax.text(0.97, 0.04, f"WT (n={len(wt_wide)})", transform=ax.transAxes, ha="right",
        va="bottom", fontsize=FA, fontweight="bold", color="black")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

ax = axes[1]
x_pos_tma = np.arange(2) + 1
tma_pair_vals = tma_wide[TMA_PAIR].dropna()
vp = ax.violinplot([tma_pair_vals[c].values for c in TMA_PAIR], positions=x_pos_tma,
                   showmeans=False, showmedians=False, showextrema=False, widths=0.85)
for body, ct in zip(vp["bodies"], TMA_PAIR):
    body.set_facecolor(PAL[ct]); body.set_edgecolor("none")
    body.set_alpha(0.30); body.set_linewidth(0)
_draw_paired_traces(ax, tma_pair_vals, x_pos_tma, TMA_PAIR, line_alpha=0.18,
                    line_lw=0.3, dot_size=5, dot_edge=0.2)
ax.set_xticks(x_pos_tma); ax.set_xticklabels(TMA_LABEL); ax.set_xlim(0.45, 2.55)
ax.text(0.97, 0.04, f"TMA (n={len(tma_pair_vals)})", transform=ax.transAxes, ha="right",
        va="bottom", fontsize=FA, fontweight="bold", color="black")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")


def _bracket(ax, x1, x2, y, text, h_frac=0.04):
    yr = ax.get_ylim()
    h = (yr[1] - yr[0]) * h_frac
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", lw=0.4, clip_on=False)
    ax.text((x1 + x2) / 2, y + h * 1.2, text, ha="center", va="bottom",
            fontsize=FN, color="black", clip_on=False)


all_vals = np.concatenate([wt_wide[EPI_ORDER].values.ravel(),
                           tma_wide[EPI_ORDER].dropna().values.ravel()])
ymin, ymax = float(np.min(all_vals)), float(np.max(all_vals))
span = ymax - ymin
axes[0].set_ylim(ymin - span * 0.05, ymax + span * 0.18)
axes[1].set_ylim(ymin - span * 0.05, ymax + span * 0.18)
y_bracket = ymax + span * 0.04
_bracket(axes[0], 1, 3, y_bracket, fmt_p(wt_p_AB))
_bracket(axes[1], 1, 2, y_bracket, fmt_p(tma_p_AB))

fig.tight_layout(pad=0.3, w_pad=0.6)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved:\n  {os.path.basename(OUT_PNG)}\n  {os.path.basename(OUT_SVG)}")
