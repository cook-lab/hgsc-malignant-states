#!/usr/bin/env python3
"""
Figure 4i — Local cell density around SecA, Intermediate and SecB cells (paired)
=================================================================================
PURPOSE
    Per-sample test of the Fig 4f density curve, with the TMA as validation. Local cell
    density = all cells within 50 um of a secretory cell, per mm^2 (the Fig 4f
    definition); each sample's value is the median over that class's cells. Left: the 8
    whole tissues (SecA -> Intermediate -> SecB paired traces). Right: TMA patients (SecA
    vs SecB boxes). Linear y axes, each fitted to its own points. Paired Wilcoxon SecA vs
    SecB per arm.

INPUTS
    - output_root/singler_recompute/per_patient_local_density_SINGLER.csv
      (figure4/14_fig4i_local_density.py; ldens_{SecA,Int,SecB} plotted, med_nb_* (the median
      neighbour counts they are made from) tested, as their ties are exact; "source" == "WT" | "TMA");
      env PER_PATIENT_CSV overrides the path.

OUTPUTS
    - figures_dir/figure4/xenium_epi_density_paired.{png,svg}
    - figures_dir/figure4/xenium_epi_density_paired_stats.csv (n, medians, units lower around SecB, p per arm)
    - figures_dir/figure4/xenium_epi_density_paired_units.csv (the plotted values and cell counts per sample;
      units keep the names of the input table, TMA patients as '<TMA>|<patient number>', not study codes)

MANUSCRIPT PANEL(S): Fig 4i.

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
import matplotlib.ticker
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
PER_PATIENT = path("output_root", "singler_recompute", "per_patient_local_density_SINGLER.csv")
_pp_override = os.environ.get("PER_PATIENT_CSV", "")
if _pp_override and os.path.exists(_pp_override):
    PER_PATIENT = _pp_override
OUT_PNG = path("figures_dir", "figure4", "xenium_epi_density_paired.png")
OUT_SVG = path("figures_dir", "figure4", "xenium_epi_density_paired.svg")
OUT_STATS = path("figures_dir", "figure4", "xenium_epi_density_paired_stats.csv")
OUT_UNITS = path("figures_dir", "figure4", "xenium_epi_density_paired_units.csv")

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

EPI_COLS = ["ldens_SecA", "ldens_Int", "ldens_SecB"]
EPI_LABEL = ["SecA", "Int.", "SecB"]
PAL = {"ldens_SecA": "#E6A141", "ldens_Int": "#C08E48", "ldens_SecB": "#9A7D55"}

# ---------- Load ----------
d = pd.read_csv(PER_PATIENT)
NB = {c: c.replace("ldens_", "med_nb_") for c in EPI_COLS}   # tested columns (density = med_nb / circle area)
wt = d[d["source"] == "WT"][EPI_COLS + list(NB.values())].dropna().copy()
tma = d[d["source"] == "TMA"][[EPI_COLS[0], EPI_COLS[2], NB[EPI_COLS[0]], NB[EPI_COLS[2]]]].dropna().copy()   # SecA and SecB only
print(f"  WT: n={len(wt)}; TMA: n={len(tma)}")


def paired_p(df, a, b):
    a, b = NB[a], NB[b]
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
fig, axes = plt.subplots(1, 2, figsize=(65 / 25.4, 35 / 25.4), sharey=False,
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
ax.set_ylabel("Local cell density\n(cells / mm², 50 µm)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

ax = axes[1]
x_pos_tma = np.arange(2) + 1
tma_pair = tma[TMA_PAIR_COLS]
bp = ax.boxplot([tma_pair[c].values for c in TMA_PAIR_COLS], positions=x_pos_tma,
                widths=0.55, showfliers=False, patch_artist=True,
                medianprops=dict(color="black", lw=0.9),
                whiskerprops=dict(color="black", lw=0.6),
                capprops=dict(color="black", lw=0.6), boxprops=dict(lw=0.5))
for patch, ct in zip(bp["boxes"], TMA_PAIR_COLS):
    patch.set_facecolor(PAL[ct]); patch.set_alpha(0.55); patch.set_edgecolor("black")
_rng = np.random.default_rng(42)
for j, ct in enumerate(TMA_PAIR_COLS):
    _v = tma_pair[ct].values
    ax.scatter(x_pos_tma[j] + _rng.uniform(-0.13, 0.13, len(_v)), _v, s=2, color=PAL[ct],
               alpha=0.25, edgecolors="none", zorder=1)
ax.set_xticks(x_pos_tma); ax.set_xticklabels(TMA_LABEL); ax.set_xlim(0.45, 2.55)
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


# linear axes, each fitted to its own points; every point stays inside the frame
for ax_, vals, pad_lo, pad_hi, x1, x2, p_ in ((axes[0], wt[EPI_COLS].values.ravel(), 0.10, 0.32, 1, 3, wt_p_AB),
                                              (axes[1], tma_pair.values.ravel(), 0.04, 0.16, 1, 2, tma_p_AB)):
    lo, hi = float(vals.min()), float(vals.max()); span = hi - lo
    ax_.set_ylim(lo - span * pad_lo, hi + span * pad_hi)
    ax_.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    _bracket(ax_, x1, x2, hi + span * (0.06 if ax_ is axes[0] else 0.03), fmt_p(p_))

stats = []
for arm, df in (("whole tissue", wt), ("TMA patients", tma)):
    a, b = df[EPI_COLS[0]], df[EPI_COLS[2]]
    stats.append({"arm": arm, "n": len(df), "median_SecA": a.median(), "median_SecB": b.median(),
                  "n_lower_around_SecB": int((b < a).sum()), "wilcoxon_p": paired_p(df, EPI_COLS[0], EPI_COLS[2])})
pd.DataFrame(stats).to_csv(OUT_STATS, index=False)
UNIT_COLS = ["source", "patient_id"] + EPI_COLS + ["n_SecA", "n_Int", "n_SecB"]
pd.concat([d.loc[wt.index, UNIT_COLS], d.loc[tma.index, UNIT_COLS]]).to_csv(OUT_UNITS, index=False)
print(pd.DataFrame(stats).to_string(index=False))

fig.tight_layout(pad=0.3, w_pad=0.6)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved:\n  {os.path.basename(OUT_PNG)}\n  {os.path.basename(OUT_SVG)}")
