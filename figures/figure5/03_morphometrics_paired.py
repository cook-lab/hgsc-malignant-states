#!/usr/bin/env python3
# ============================================================================
# Figure 5D,5E — Paired epithelial morphometrics across SecA/Intermediate/SecB
# ----------------------------------------------------------------------------
# PURPOSE
#   Paired-sample comparison of epithelial morphometric features (nuclear area,
#   nuclear perimeter, N:C ratio) across the SecA -> Intermediate -> SecB
#   polarisation axis.
#     WT panel:  SecA -> Intermediate -> SecB paired traces (whole-tissue
#                samples; sample-level medians).
#     TMA panel: SecA -> SecB endpoints with cell-type violins behind
#                per-patient paired traces (patient-level medians).
#   Eligibility uses the 06f-corrected polarisation labels baked into the 33a
#   morphometric cache; samples/patients missing an epitype or with <30 cells
#   in any compared epitype are dropped (eligible_paired flag).
#
# INPUTS
#   data_root/2026_final_xenium_analysis/output/33_morphometrics/
#     per_sample_summary_wt.csv      (WT, group_id = sample_key)
#     per_patient_summary_tma.csv    (TMA, group_id = patient_id)
#     wt_cell_counts.csv             (eligibility flag)
#     tma_cell_counts.csv            (eligibility flag)
#
# OUTPUTS
#   figures_dir/figure5/xenium_nuc_area_paired.{png,svg}
#   figures_dir/figure5/xenium_nuc_perimeter_paired.{png,svg}
#   figures_dir/figure5/xenium_nc_ratio_paired.{png,svg}
#
# MANUSCRIPT PANEL(S): Fig 5D (nuclear area), Fig 5E (N:C ratio)
#                      (perimeter panel is supporting / not in main figure)
# RUNTIME TIER: fast (reads summary CSVs only)
# ============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

# --- Config (script is 2 levels deep: figures/figure5/) ---------------------
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED

np.random.seed(SEED)

FIG_DIR  = path("figures_dir", "figure5")
MORPH_DIR = path("data_root", "2026_final_xenium_analysis",
                 "output", "33_morphometrics")
WT_CSV   = os.path.join(MORPH_DIR, "per_sample_summary_wt.csv")
TMA_CSV  = os.path.join(MORPH_DIR, "per_patient_summary_tma.csv")
WT_CNT   = os.path.join(MORPH_DIR, "wt_cell_counts.csv")
TMA_CNT  = os.path.join(MORPH_DIR, "tma_cell_counts.csv")

# ── Style ─────────────────────────────────────────────────────────────────
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
PAL = {"SecA epithelium":         "#E6A141",
       "Intermediate epithelium": "#C08E48",
       "SecB epithelium":         "#9A7D55"}

# Three features to plot: column name → y-axis label, filename slug
FEATURES = [
    ("nuc_area",      r"Nuclear area (μm$^2$)", "nuc_area"),
    ("nuc_perimeter", r"Nuclear perimeter (μm)", "nuc_perimeter"),
    ("nc_ratio",      "Nuclear : cytoplasmic ratio",  "nc_ratio"),
]

# ── Load + filter ─────────────────────────────────────────────────────────
print("Loading morphometric per-sample / per-patient summaries ...")
wt_long  = pd.read_csv(WT_CSV)
tma_long = pd.read_csv(TMA_CSV)
wt_cnt   = pd.read_csv(WT_CNT)
tma_cnt  = pd.read_csv(TMA_CNT)
# Standardize cell_label "Transitioning epithelium" -> "Intermediate epithelium" (deposited cache predates the rename).
wt_long["cell_label"]  = wt_long["cell_label"].replace({"Transitioning epithelium": "Intermediate epithelium"})
tma_long["cell_label"] = tma_long["cell_label"].replace({"Transitioning epithelium": "Intermediate epithelium"})

wt_eligible  = set(wt_cnt.loc[wt_cnt["eligible_paired"] == True,
                               "sample_key"].astype(str))
tma_cnt = tma_cnt.dropna(subset=["group_id"])
tma_cnt["group_id"] = tma_cnt["group_id"].astype(str)
tma_eligible = set(tma_cnt.loc[tma_cnt["eligible_paired"] == True,
                                "group_id"])

wt_long["group_id"]  = wt_long["group_id"].astype(str)
tma_long["group_id"] = tma_long["group_id"].astype(str)
wt_long  = wt_long[wt_long["group_id"].isin(wt_eligible)]
tma_long = tma_long[tma_long["group_id"].isin(tma_eligible)]
print(f"  WT eligible samples:  {len(wt_eligible)}")
print(f"  TMA eligible patients: {len(tma_eligible)}")


def pivot_wide(sub, value_col):
    return (sub.pivot_table(index="group_id", columns="cell_label",
                            values=value_col)
              .reindex(columns=EPI_ORDER)
              .dropna())


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


def _draw_paired_traces(ax, wide, x_pos, cols, line_alpha, line_lw,
                        dot_size, dot_alpha=0.95, dot_edge=0.3):
    vals = wide[cols].values
    for row in vals:
        ax.plot(x_pos, row, color="grey",
                lw=line_lw, alpha=line_alpha, zorder=2)
    for j, ct in enumerate(cols):
        ax.scatter(np.full(len(vals), x_pos[j]), vals[:, j],
                   s=dot_size, color=PAL[ct],
                   edgecolor="white", linewidth=dot_edge,
                   alpha=dot_alpha, zorder=3)


def _bracket(ax, x1, x2, y, text, h_frac=0.04):
    yr = ax.get_ylim()
    h = (yr[1] - yr[0]) * h_frac
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y],
            color="black", lw=0.4, clip_on=False)
    ax.text((x1 + x2) / 2, y + h * 1.4, text,
            ha="center", va="bottom", fontsize=FN, color="black",
            clip_on=False)


def make_panel(value_col, ylab, slug):
    wt_w  = pivot_wide(wt_long,  value_col)
    tma_w = pivot_wide(tma_long, value_col)
    print(f"\n[{value_col}] WT n={len(wt_w)}, TMA n={len(tma_w)}")

    wt_p_AB  = paired_p(wt_w,  EPI_ORDER[0], EPI_ORDER[2])
    tma_p_AB = paired_p(tma_w, EPI_ORDER[0], EPI_ORDER[2])
    print(f"  WT  paired SecA vs SecB: {fmt_p(wt_p_AB)}")
    print(f"  TMA paired SecA vs SecB: {fmt_p(tma_p_AB)}")

    fig_w_mm, fig_h_mm = 65, 35
    fig, axes = plt.subplots(1, 2,
                             figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4),
                             sharey=True,
                             gridspec_kw={"width_ratios": [1.1, 1.0]})

    TMA_PAIR  = [EPI_ORDER[0], EPI_ORDER[2]]
    TMA_LABEL = ["SecA", "SecB"]

    # ── WT
    ax = axes[0]
    x_pos_wt = np.arange(3) + 1
    _draw_paired_traces(ax, wt_w, x_pos_wt, EPI_ORDER,
                        line_alpha=0.55, line_lw=0.7, dot_size=14)
    ax.set_xticks(x_pos_wt)
    ax.set_xticklabels(EPI_LABEL)
    ax.set_xlim(0.55, 3.45)
    ax.set_ylabel(ylab)
    ax.text(0.97, 0.04, f"WT (n={len(wt_w)})",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=FA, fontweight="bold", color="black")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
    ax.tick_params(width=0.4, length=2, colors="black")

    # ── TMA: violin + paired traces (SecA → SecB)
    ax = axes[1]
    x_pos_tma = np.arange(2) + 1
    tma_pair = tma_w[TMA_PAIR].dropna()

    vp_data = [tma_pair[c].values for c in TMA_PAIR]
    vp = ax.violinplot(vp_data, positions=x_pos_tma,
                       showmeans=False, showmedians=False, showextrema=False,
                       widths=0.85)
    for body, ct in zip(vp["bodies"], TMA_PAIR):
        body.set_facecolor(PAL[ct]); body.set_edgecolor("none")
        body.set_alpha(0.30); body.set_linewidth(0)

    _draw_paired_traces(ax, tma_pair, x_pos_tma, TMA_PAIR,
                        line_alpha=0.18, line_lw=0.3,
                        dot_size=5, dot_edge=0.2)

    ax.set_xticks(x_pos_tma)
    ax.set_xticklabels(TMA_LABEL)
    ax.set_xlim(0.45, 2.55)
    ax.text(0.97, 0.04, f"TMA (n={len(tma_pair)})",
            transform=ax.transAxes, ha="right", va="bottom",
            fontsize=FA, fontweight="bold", color="black")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
    ax.tick_params(width=0.4, length=2, colors="black")

    # Shared y-axis sized to combined data + headroom for bracket
    all_vals = np.concatenate([wt_w[EPI_ORDER].values.ravel(),
                               tma_w[EPI_ORDER].values.ravel()])
    ymin = float(np.min(all_vals))
    ymax = float(np.max(all_vals))
    span = ymax - ymin
    y_lo = ymin - span * 0.05
    y_hi = ymax + span * 0.20
    axes[0].set_ylim(y_lo, y_hi)
    axes[1].set_ylim(y_lo, y_hi)

    y_bracket = ymax + span * 0.04
    _bracket(axes[0], 1, 3, y_bracket, fmt_p(wt_p_AB))
    _bracket(axes[1], 1, 2, y_bracket, fmt_p(tma_p_AB))

    fig.tight_layout(pad=0.3, w_pad=0.6)

    out_png = os.path.join(FIG_DIR, f"xenium_{slug}_paired.png")
    out_svg = os.path.join(FIG_DIR, f"xenium_{slug}_paired.svg")
    fig.savefig(out_png, bbox_inches="tight", dpi=450)
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {os.path.basename(out_png)}")
    print(f"         {os.path.basename(out_svg)}")


for col, ylab, slug in FEATURES:
    make_panel(col, ylab, slug)

print("\nDone.")
