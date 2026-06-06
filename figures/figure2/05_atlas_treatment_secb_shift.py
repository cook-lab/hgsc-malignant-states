#!/usr/bin/env python3
"""
Figure 2F,G — Primary epitype composition by treatment + matched pre/post SecB shift
====================================================================================
PURPOSE
    Left: stacked epitype composition (%) for primary tumours, pre-treatment vs
    post-chemo. Right: per-patient mean SecB NMF Factor-2 loading, pre vs post
    (>=100 epithelial cells each timepoint), with a paired Wilcoxon signed-rank
    significance bracket.

INPUTS
    - data_fig1i_treatment_proportions.csv :
        output_root/fig_data_fig1/data_fig1i_treatment_proportions.csv (bars)
    - meta.parquet :
        output_root/fig_secretory_polarization/data/meta.parquet
        (patient_id, treatment_status, Factor_2; produced by the prep helper)

OUTPUTS
    - figures_dir/atlas_treatment_secb_shift.{svg,png}

MANUSCRIPT PANEL(S): Fig 2F (left), Fig 2G (right).

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
from matplotlib.patches import Patch
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

# ---------- Paths ----------
BARS_CSV = path("output_root", "fig_data_fig1", "data_fig1i_treatment_proportions.csv")
META_PATH = path("output_root", "fig_secretory_polarization", "data", "meta.parquet")
OUT_SVG = path("figures_dir", "atlas_treatment_secb_shift.svg")
OUT_PNG = path("figures_dir", "atlas_treatment_secb_shift.png")

MIN_CELLS_PAIRED = 100

# ---------- Style ----------
PALETTE = {"Ciliated": "#E07850", "SecA": "#E6A141",
           "Intermediate": "#C08E48", "SecB": "#9A7D55"}
ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]
INCREASE = "#D65146"
DECREASE = "#4575B4"
MEDIAN_C = "#111111"

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         7,
    "axes.titlesize":    9,
    "axes.labelsize":    8,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "legend.fontsize":   7,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "svg.fonttype":      "none",
    "axes.linewidth":    0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size":  2.5,
    "ytick.major.size":  2.5,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.facecolor":    "white",
    "figure.facecolor":  "white",
})


def p_to_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def add_bracket(ax, x1, x2, y, h, stars, fontsize=7):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.8, color="black", clip_on=False)
    ax.text((x1 + x2) / 2, y + h, stars, ha="center", va="bottom",
            fontsize=fontsize, fontweight="bold")


# ---------- Load ----------
bars = pd.read_csv(BARS_CSV).set_index("treatment_status")
meta = pd.read_parquet(META_PATH, columns=["patient_id", "treatment_status", "Factor_2"])
meta = meta[meta["treatment_status"].isin(["pre-treatment", "post-chemotherapy"])]
counts = meta.groupby(["patient_id", "treatment_status"], observed=True).size()
means = meta.groupby(["patient_id", "treatment_status"], observed=True)["Factor_2"].mean()
per_pt = pd.concat([counts.rename("n_cells"), means.rename("mean_f2")], axis=1).reset_index()
per_pt = per_pt[per_pt["n_cells"] >= MIN_CELLS_PAIRED]

paired = per_pt.pivot_table(index="patient_id", columns="treatment_status",
                            values="mean_f2").dropna()
paired.columns.name = None
pre = paired["pre-treatment"].values
post = paired["post-chemotherapy"].values
n_paired = len(paired)
n_inc = int((post > pre).sum())
wilc_p = stats.wilcoxon(pre, post).pvalue
median_pre, median_post = float(np.median(pre)), float(np.median(post))

N_PRE, N_POST = 116, 24

# ---------- Figure ----------
fig, (ax_b, ax_p) = plt.subplots(1, 2, figsize=(7.0, 3.4),
                                 gridspec_kw=dict(width_ratios=[1.0, 1.3], wspace=0.45))

grp_labels = ["pre-treatment", "post-chemotherapy"]
grp_display = {"pre-treatment": f"Pre-treatment\n(n={N_PRE})",
               "post-chemotherapy": f"Post-chemo\n(n={N_POST})"}
bar_w = 0.7
for xi, grp in enumerate(grp_labels):
    bottom = 0.0
    for ct in ORDER:
        val = float(bars.loc[grp, ct])
        ax_b.bar(xi, val, bottom=bottom, width=bar_w, color=PALETTE[ct],
                 edgecolor="#333333", linewidth=0.3)
        if val > 6:
            ax_b.text(xi, bottom + val / 2, f"{val:.0f}%", ha="center", va="center",
                      fontsize=7, fontweight="bold", color="white" if val > 12 else "#333333")
        bottom += val
ax_b.set_xticks(range(len(grp_labels)))
ax_b.set_xticklabels([grp_display[g] for g in grp_labels])
ax_b.set_ylim(0, 100)
ax_b.set_ylabel("Composition (%)")
ax_b.set_title("Epitype composition", fontweight="bold")

x_pre, x_post = 0.0, 1.0
jitter = np.random.default_rng(SEED).uniform(-0.03, 0.03, n_paired)
for i, (pid, row) in enumerate(paired.iterrows()):
    p0, p1 = row["pre-treatment"], row["post-chemotherapy"]
    color = INCREASE if p1 > p0 else DECREASE
    ax_p.plot([x_pre + jitter[i], x_post + jitter[i]], [p0, p1], color=color,
              linewidth=0.7, alpha=0.75, zorder=2)
    ax_p.scatter([x_pre + jitter[i], x_post + jitter[i]], [p0, p1], s=16, color=color,
                 alpha=0.85, zorder=3, edgecolors="white", linewidths=0.3)

ax_p.plot([x_pre, x_post], [median_pre, median_post], color=MEDIAN_C, linewidth=1.1,
          linestyle=(0, (3, 2)), zorder=4)
ax_p.scatter([x_pre, x_post], [median_pre, median_post], s=30, color=MEDIAN_C, zorder=5)
ax_p.text(x_pre - 0.08, median_pre, f"{median_pre:.3f}", ha="right", va="center",
          fontsize=7, fontweight="bold")
ax_p.text(x_post + 0.08, median_post, f"{median_post:.3f}", ha="left", va="center",
          fontsize=7, fontweight="bold")

all_vals = np.concatenate([pre, post])
y_max = np.max(all_vals)
bracket_y, bracket_h = y_max + 0.02, 0.01
add_bracket(ax_p, x_pre, x_post, bracket_y, bracket_h, p_to_stars(wilc_p), fontsize=9)
ax_p.text((x_pre + x_post) / 2, bracket_y + bracket_h + 0.02,
          f"n={n_paired} paired, {n_inc}/{n_paired} increase",
          ha="center", va="bottom", fontsize=6, color="#555555")
ax_p.set_ylim(ax_p.get_ylim()[0], bracket_y + bracket_h + 0.06)
ax_p.set_xticks([x_pre, x_post])
ax_p.set_xticklabels(["Pre-treatment", "Post-chemo"])
ax_p.set_xlim(-0.5, 1.5)
ax_p.set_ylabel("Mean SecB NMF factor loading")
ax_p.set_title("Patient-matched SecB shift", fontweight="bold")

handles = [Patch(facecolor=PALETTE[c], edgecolor="#333333", linewidth=0.3, label=c)
           for c in ["Ciliated", "SecA", "Intermediate", "SecB"]]
fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02), frameon=False)
fig.tight_layout(rect=(0, 0.05, 1, 1))

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight")
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight")
plt.close(fig)
print("Done.")
