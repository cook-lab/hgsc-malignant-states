#!/usr/bin/env python3
"""
Figure 4J — Secretory epithelial density across the gradient (paired)
=====================================================================
PURPOSE
    Paired-sample comparison of secretory epithelial density (cells/mm^2) across
    SecA -> Intermediate -> SecB, for the 8 whole-tissue samples (3-position
    paired traces) and the TMA cohort (SecA->SecB with violins). Log y-axis;
    Wilcoxon paired SecA-vs-SecB bracket per panel.

INPUTS
    - output_root/10_clinical_v2/per_patient_features_v2.csv
      (dens_<epitype>_epithelium columns; "source" == "WT" | "TMA")

OUTPUTS
    - figures_dir/figure4/xenium_epi_density_paired.{png,svg}

MANUSCRIPT PANEL(S): Fig 4J.

RUNTIME TIER: fast.

NOTE: epitype label standardized "Transitioning" -> "Intermediate"; the upstream
per-patient cache column name (dens_Transitioning_epithelium) is read as-is.
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
PER_PATIENT = path("data_root", "2026_final_xenium_analysis", "output", "10_clinical_v2", "per_patient_features_v2.csv")
OUT_PNG = path("figures_dir", "figure4", "xenium_epi_density_paired.png")
OUT_SVG = path("figures_dir", "figure4", "xenium_epi_density_paired.svg")

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

EPI_COLS = ["dens_SecA_epithelium", "dens_Transitioning_epithelium", "dens_SecB_epithelium"]
EPI_LABEL = ["SecA", "Int.", "SecB"]
PAL = {"dens_SecA_epithelium": "#E6A141", "dens_Transitioning_epithelium": "#C08E48",
       "dens_SecB_epithelium": "#9A7D55"}

# ---------- Load ----------
d = pd.read_csv(PER_PATIENT)
wt = d[d["source"] == "WT"][EPI_COLS].dropna().copy()
tma = d[d["source"] == "TMA"][EPI_COLS].dropna().copy()
print(f"  WT: n={len(wt)}; TMA: n={len(tma)}")


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


wt_p_AB = paired_p(wt, EPI_COLS[0], EPI_COLS[2])
tma_p_AB = paired_p(tma, EPI_COLS[0], EPI_COLS[2])

# ---------- Plot ----------
fig, axes = plt.subplots(1, 2, figsize=(65 / 25.4, 35 / 25.4), sharey=True,
                         gridspec_kw={"width_ratios": [1.1, 1.0]})
TMA_PAIR_COLS = [EPI_COLS[0], EPI_COLS[2]]
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
_draw_paired_traces(ax, wt, x_pos_wt, EPI_COLS, line_alpha=0.55, line_lw=0.7, dot_size=14)
ax.set_xticks(x_pos_wt); ax.set_xticklabels(EPI_LABEL); ax.set_xlim(0.55, 3.45)
ax.set_yscale("log")
ax.set_ylabel("Density (cells / mm²)")
ax.text(0.97, 0.04, f"WT (n={len(wt)})", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=FA, fontweight="bold", color="black")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

ax = axes[1]
x_pos_tma = np.arange(2) + 1
tma_pair = tma[TMA_PAIR_COLS].dropna()
vp = ax.violinplot([tma_pair[c].values for c in TMA_PAIR_COLS], positions=x_pos_tma,
                   showmeans=False, showmedians=False, showextrema=False, widths=0.85)
for body, ct in zip(vp["bodies"], TMA_PAIR_COLS):
    body.set_facecolor(PAL[ct]); body.set_edgecolor("none")
    body.set_alpha(0.30); body.set_linewidth(0)
_draw_paired_traces(ax, tma_pair, x_pos_tma, TMA_PAIR_COLS, line_alpha=0.18,
                    line_lw=0.3, dot_size=5, dot_edge=0.2)
ax.set_xticks(x_pos_tma); ax.set_xticklabels(TMA_LABEL); ax.set_xlim(0.45, 2.55)
ax.text(0.97, 0.04, f"TMA (n={len(tma_pair)})", transform=ax.transAxes, ha="right",
        va="bottom", fontsize=FA, fontweight="bold", color="black")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")


def _bracket(ax, x1, x2, y, text, h_frac=0.04):
    yr = ax.get_ylim()
    if ax.get_yscale() == "log":
        log_h = (np.log10(yr[1]) - np.log10(yr[0])) * h_frac
        y_top = 10 ** (np.log10(y) + log_h)
        y_text = 10 ** (np.log10(y) + log_h * 1.4)
    else:
        h = (yr[1] - yr[0]) * h_frac
        y_top = y + h
        y_text = y + h * 1.4
    ax.plot([x1, x1, x2, x2], [y, y_top, y_top, y], color="black", lw=0.4, clip_on=False)
    ax.text((x1 + x2) / 2, y_text, text, ha="center", va="bottom", fontsize=FN,
            color="black", clip_on=False)


all_vals = np.concatenate([wt[EPI_COLS].values.ravel(), tma[EPI_COLS].values.ravel()])
ymin = float(np.min(all_vals[all_vals > 0]))
ymax = float(np.max(all_vals))
log_span = np.log10(ymax) - np.log10(ymin)
axes[0].set_ylim(10 ** (np.log10(ymin) - log_span * 0.05), 10 ** (np.log10(ymax) + log_span * 0.20))
axes[1].set_ylim(10 ** (np.log10(ymin) - log_span * 0.05), 10 ** (np.log10(ymax) + log_span * 0.20))
y_bracket = 10 ** (np.log10(ymax) + log_span * 0.04)
_bracket(axes[0], 1, 3, y_bracket, fmt_p(wt_p_AB))
_bracket(axes[1], 1, 2, y_bracket, fmt_p(tma_p_AB))

fig.tight_layout(pad=0.3, w_pad=0.6)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved:\n  {os.path.basename(OUT_PNG)}\n  {os.path.basename(OUT_SVG)}")
