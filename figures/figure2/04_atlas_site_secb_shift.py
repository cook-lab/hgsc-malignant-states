#!/usr/bin/env python3
"""
Figure 2D,E — Epitype composition by site + patient-matched SecB across sites
=============================================================================
PURPOSE
    Left: stacked epitype composition (%) by metastatic site (primary / ascites /
    metastasis). Right: per-patient SecB proportion across sites (ascites-matched
    patients), with paired Wilcoxon signed-rank significance brackets.

INPUTS
    - panel_b_site_proportions.csv   (Fig 2D bars)
    - panel_c_paired_site.csv        (Fig 2E paired)
        both under output_root/fig_secretory_polarization/data/
        (produced by figures/_prep/fig_secretory_polarization_00_prepare_data.py)

OUTPUTS
    - figures_dir/atlas_site_secb_shift.{svg,png}

MANUSCRIPT PANEL(S): Fig 2D (left), Fig 2E (right).

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
DATA_DIR = path("output_root", "fig_secretory_polarization", "data")
BARS_CSV = os.path.join(DATA_DIR, "panel_b_site_proportions.csv")
PAIRED_CSV = os.path.join(DATA_DIR, "panel_c_paired_site.csv")
OUT_SVG = path("figures_dir", "atlas_site_secb_shift.svg")
OUT_PNG = path("figures_dir", "atlas_site_secb_shift.png")

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
bars = pd.read_csv(BARS_CSV).set_index("metastatic_site")
paired_long = pd.read_csv(PAIRED_CSV)

n_per_site = paired_long.groupby("metastatic_site")["patient_id"].nunique()
N_PRI = int(n_per_site.get("primary", 0))
N_ASC = int(n_per_site.get("ascites", 0))
N_META = int(n_per_site.get("metastasis", 0))

SITE_ORDER = ["primary", "ascites", "metastasis"]
asc_pts = paired_long.loc[paired_long["metastatic_site"] == "ascites", "patient_id"].unique()
paired_sub = paired_long[paired_long["patient_id"].isin(asc_pts)].copy()
n_paired_pts = paired_sub["patient_id"].nunique()
paired_wide = (paired_sub.pivot_table(index="patient_id", columns="metastatic_site",
                                      values="secb_prop").reindex(columns=SITE_ORDER))


def paired_wilcoxon(df, a, b):
    d = df[[a, b]].dropna()
    if len(d) < 2:
        return np.nan, 0
    return float(stats.wilcoxon(d[a].values, d[b].values).pvalue), int(len(d))


p_pa, _ = paired_wilcoxon(paired_wide, "primary", "ascites")
p_am, _ = paired_wilcoxon(paired_wide, "ascites", "metastasis")
p_pm, _ = paired_wilcoxon(paired_wide, "primary", "metastasis")
medians = np.array([paired_wide[s].median() for s in SITE_ORDER])

# ---------- Figure ----------
fig, (ax_b, ax_p) = plt.subplots(1, 2, figsize=(7.6, 3.4),
                                 gridspec_kw=dict(width_ratios=[1.0, 1.3], wspace=0.45))

grp_display = {"primary": f"Primary\n(n={N_PRI})", "ascites": f"Ascites\n(n={N_ASC})",
               "metastasis": f"Metastasis\n(n={N_META})"}
bar_w = 0.7
for xi, site in enumerate(SITE_ORDER):
    bottom = 0.0
    for ct in ORDER:
        val = float(bars.loc[site, ct])
        ax_b.bar(xi, val, bottom=bottom, width=bar_w, color=PALETTE[ct],
                 edgecolor="#333333", linewidth=0.3)
        if val > 6:
            ax_b.text(xi, bottom + val / 2, f"{val:.0f}%", ha="center", va="center",
                      fontsize=7, fontweight="bold", color="white" if val > 12 else "#333333")
        bottom += val
ax_b.set_xticks(range(len(SITE_ORDER)))
ax_b.set_xticklabels([grp_display[s] for s in SITE_ORDER])
ax_b.set_ylim(0, 100)
ax_b.set_ylabel("Composition (%)")
ax_b.set_title("Epitype composition", fontweight="bold")

x_pos = np.array([0.0, 1.0, 2.0])
rng = np.random.default_rng(SEED)
jitter = rng.uniform(-0.04, 0.04, n_paired_pts)
for i, (pid, row) in enumerate(paired_wide.iterrows()):
    vals = row.values
    for k in range(len(SITE_ORDER) - 1):
        v1, v2 = vals[k], vals[k + 1]
        if np.isnan(v1) or np.isnan(v2):
            continue
        seg_color = INCREASE if v2 > v1 else DECREASE
        ax_p.plot([x_pos[k] + jitter[i], x_pos[k + 1] + jitter[i]],
                  [v1 * 100, v2 * 100], color=seg_color, linewidth=0.7, alpha=0.75, zorder=2)
    for k, v in enumerate(vals):
        if np.isnan(v):
            continue
        if k < len(SITE_ORDER) - 1 and not np.isnan(vals[k + 1]):
            ref = INCREASE if vals[k + 1] > v else DECREASE
        elif k > 0 and not np.isnan(vals[k - 1]):
            ref = INCREASE if v > vals[k - 1] else DECREASE
        else:
            ref = MEDIAN_C
        ax_p.scatter(x_pos[k] + jitter[i], v * 100, s=16, color=ref, alpha=0.85,
                     zorder=3, edgecolors="white", linewidths=0.3)

ax_p.plot(x_pos, medians * 100, color=MEDIAN_C, linewidth=1.1, linestyle=(0, (3, 2)), zorder=4)
ax_p.scatter(x_pos, medians * 100, s=30, color=MEDIAN_C, zorder=5)
for k, (off, ha) in enumerate([(-0.10, "right"), (0.0, "center"), (0.10, "left")]):
    if ha == "center":
        ax_p.text(x_pos[k] + off, medians[k] * 100 + 3.5, f"{medians[k]*100:.1f}%",
                  ha=ha, va="bottom", fontsize=7, fontweight="bold")
    else:
        ax_p.text(x_pos[k] + off, medians[k] * 100, f"{medians[k]*100:.1f}%",
                  ha=ha, va="center", fontsize=7, fontweight="bold")

all_secb = paired_wide.values.flatten()
all_secb = all_secb[~np.isnan(all_secb)] * 100
y_max = np.max(all_secb)
bracket_h, gap = 2.0, 2.5
y1 = y_max + 4
add_bracket(ax_p, x_pos[0], x_pos[1], y1, bracket_h, p_to_stars(p_pa), fontsize=8)
y2 = y1 + bracket_h + gap + 2
add_bracket(ax_p, x_pos[1], x_pos[2], y2, bracket_h, p_to_stars(p_am), fontsize=8)
y3 = y2 + bracket_h + gap + 2
add_bracket(ax_p, x_pos[0], x_pos[2], y3, bracket_h, p_to_stars(p_pm), fontsize=8)
ax_p.text(1.0, y3 + bracket_h + 4, f"n={n_paired_pts} ascites-matched patients",
          ha="center", va="bottom", fontsize=6, color="#555555")
ax_p.set_ylim(ax_p.get_ylim()[0], y3 + bracket_h + 10)
ax_p.set_xticks(x_pos)
ax_p.set_xticklabels(["Primary", "Ascites", "Metastasis"])
ax_p.set_xlim(-0.5, 2.5)
ax_p.set_ylabel("Mean SecB proportion (%)")
ax_p.set_title("Patient-matched SecB across sites", fontweight="bold")

handles = [Patch(facecolor=PALETTE[c], edgecolor="#333333", linewidth=0.3, label=c)
           for c in ["Ciliated", "SecA", "Intermediate", "SecB"]]
fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02), frameon=False)
fig.tight_layout(rect=(0, 0.05, 1, 1))

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight")
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight")
plt.close(fig)
print("Done.")
