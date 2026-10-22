#!/usr/bin/env python3
"""
Figure 2f — Epithelial state composition by consensusOV subtype (grouped boxes)
==============================================================================
PURPOSE
    For each consensusOV subtype (x, ordered by median SecB), side-by-side boxes
    show the per-sample % of the epithelial compartment that is SecA and SecB
    (with --no-intermediate, the published panel; without it Intermediate boxes
    are drawn too), with every sample drawn as a point. Kruskal-Wallis p across
    subtypes is given per state; the statistics table tests SecA, Intermediate
    and SecB either way.

    Samples: 96 pre-treatment adnexal samples (77 patients) with >= 500 cells and
    >= 100 epithelial cells.

    DATA SOURCE: whole-sample RAW-count pseudobulks of the final atlas and the
    per-sample state composition (load_epithelial()['state']) from
    figures/_prep/fig2f_consensusov.py (edgeR TMM -> log2(CPM+1) -> consensusOV in
    fig2f_consensusov.R). The subtype of each sample is the 25-seed consensus call
    of figures/_prep/fig2f_consensusov_multiseed.R, mapped onto the per-sample table
    here (it overrides the single-seed column of that table).

    Published run: --no-intermediate with the default --cohort primary_pretreatment;
    --cohort all_tumours (every non-healthy sample) is not in the paper.

INPUTS
    output_root/figure2/fig2f_consensusov/fig2f_per_sample_joined.csv
    output_root/figure2/fig2f_consensusov/consensusov_calls_multiseed.csv

OUTPUTS (with --no-intermediate; without it the _SecA_SecB suffix is dropped)
    figures_dir/figure2/atlas_subtype_epitype_grouped_boxes_SecA_SecB.{svg,png}
    figures_dir/figure2/atlas_subtype_epitype_grouped_boxes_stats_SecA_SecB.csv

MANUSCRIPT PANEL(S): Fig 2f.
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
from matplotlib.patches import Patch
from scipy.stats import kruskal

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
import argparse  # noqa: E402
_ap = argparse.ArgumentParser()
_ap.add_argument("--cohort", choices=["primary_pretreatment", "all_tumours"], default="primary_pretreatment")
_ap.add_argument("--no-intermediate", action="store_true",
                 help="plot SecA and SecB only (variant; Intermediate still tested in the stats table)")
_args = _ap.parse_args()
COHORT, NO_INT = _args.cohort, _args.no_intermediate
_SUB = "fig2f_consensusov" if COHORT == "primary_pretreatment" else "fig2f_consensusov_all"
_SFX = ("" if COHORT == "primary_pretreatment" else "_all_tumours") + ("_SecA_SecB" if NO_INT else "")
SAMPLE_CSV = path("output_root", "figure2", _SUB, "fig2f_per_sample_joined.csv")
OUT_SVG = path("figures_dir", "figure2", f"atlas_subtype_epitype_grouped_boxes{_SFX}.svg")
OUT_PNG = path("figures_dir", "figure2", f"atlas_subtype_epitype_grouped_boxes{_SFX}.png")
OUT_STATS = path("figures_dir", "figure2", f"atlas_subtype_epitype_grouped_boxes_stats{_SFX}.csv")
assert os.path.exists(SAMPLE_CSV), f"Missing: {SAMPLE_CSV} (run figures/_prep/fig2f_consensusov.py)"
# Subtype calls: 25-seed consensus (figures/_prep/fig2f_consensusov_multiseed.R).
CALLS_CSV = path("output_root", "figure2", _SUB, "consensusov_calls_multiseed.csv")
assert os.path.exists(CALLS_CSV), f"Missing: {CALLS_CSV} (run figures/_prep/fig2f_consensusov_multiseed.R)"

# ---------- Style (Figure 2 rule: 8pt text; 6pt ticks/categories) ----------
FA, FK = 8, 6
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "xtick.labelsize": FK, "ytick.labelsize": FK,
    "legend.fontsize": FK, "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "figure.dpi": 450, "savefig.dpi": 450,
})
# House palette shared with Fig 1k and Fig 2c.
PALETTE = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55"}
STATES = ["SecA", "Intermediate", "SecB"]
COLS = {"SecA": "pct_secA", "Intermediate": "pct_trans", "SecB": "pct_secB"}
SUBTYPE_ORDER = ["PRO", "IMR", "DIF", "MES"]          # ascending median SecB
FULL = {"PRO": "Proliferative", "IMR": "Immunoreactive",
        "DIF": "Differentiated", "MES": "Mesenchymal"}
FIGSIZE = (4.1, 2.6)                                    # row 3, beside the volcano (2g)
SEED = 42


def fmt_p(p: float) -> str:
    return "< 0.001" if p < 0.0005 else f"{p:.2g}" if p >= 0.01 else f"{p:.1g}" if p >= 0.00095 else f"{p:.2f}"


# ---------- Load ----------
samp = pd.read_csv(SAMPLE_CSV)
_calls = pd.read_csv(CALLS_CSV)
_key = lambda x: x.astype(str).str.replace(r"\.0$", "", regex=True)  # noqa: E731
samp["subtype_bulk"] = _key(samp["sample_id"]).map(dict(zip(_key(_calls["sample_id"]), _calls["subtype"])))
assert samp["subtype_bulk"].notna().all(), "samples without a consensus call"
n = samp["subtype_bulk"].value_counts()

# ---------- Statistics ----------
rows = []
for st in STATES:
    groups = [samp.loc[samp.subtype_bulk == s, COLS[st]].values for s in SUBTYPE_ORDER]
    h, p = kruskal(*groups)
    pt = samp.sort_values("n_epi", ascending=False).drop_duplicates("patient_id")  # one sample per patient
    p_pt = kruskal(*[pt.loc[pt.subtype_bulk == s_, COLS[st]].values for s_ in SUBTYPE_ORDER]).pvalue
    row = {"state": st, "kruskal_H": h, "kruskal_p": p, "kruskal_p_patient_level": p_pt,
           "n_patients": pt.patient_id.nunique()}
    for s, g in zip(SUBTYPE_ORDER, groups):
        row[f"median_{s}"] = np.median(g)
        row[f"n_{s}"] = len(g)
    rows.append(row)
stats = pd.DataFrame(rows)
stats.to_csv(OUT_STATS, index=False)
kw = dict(zip(stats.state, stats.kruskal_p))

# ---------- Plot ----------
rng = np.random.default_rng(SEED)
fig, ax = plt.subplots(figsize=FIGSIZE)
PLOT_STATES = ["SecA", "SecB"] if NO_INT else STATES
w = 0.26 if not NO_INT else 0.34
for i, s in enumerate(SUBTYPE_ORDER):
    for k, st in enumerate(PLOT_STATES):
        v = samp.loc[samp.subtype_bulk == s, COLS[st]].values
        x = i + (k - (len(PLOT_STATES) - 1) / 2) * w * 1.05
        ax.boxplot(v, positions=[x], widths=w * 0.9, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor=PALETTE[st], edgecolor="k", lw=0.5),
                   medianprops=dict(color="k", lw=0.8),
                   whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
        ax.scatter(x + rng.uniform(-w * 0.3, w * 0.3, len(v)), v, s=3, color="k",
                   lw=0, alpha=0.6, zorder=3)

ax.set_xticks(range(len(SUBTYPE_ORDER)))
ax.set_xticklabels([f"{FULL[s]}\n(n = {n[s]})" for s in SUBTYPE_ORDER])
ax.set_xlim(-0.6, len(SUBTYPE_ORDER) - 0.4)
ax.set_ylim(-2, 104)
ax.set_ylabel("% of epithelium per sample")
ax.tick_params(axis="x", length=0)
ax.legend(handles=[Patch(facecolor=PALETTE[st], edgecolor="k", lw=0.5,
                         label=f"{st} (p = {fmt_p(kw[st])})") for st in PLOT_STATES],
          frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0),
          handlelength=1, columnspacing=1.2)

fig.savefig(OUT_SVG, format="svg", bbox_inches="tight", pad_inches=0.04)
fig.savefig(OUT_PNG, format="png", bbox_inches="tight", pad_inches=0.04)
plt.close(fig)
print(stats.round(4).to_string(index=False))
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}\n  Saved: {OUT_STATS}")
