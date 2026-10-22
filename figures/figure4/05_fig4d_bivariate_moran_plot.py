#!/usr/bin/env python3
"""
Fig 4d — global bivariate Moran's I between SecA and neighbour SecB scores (plot + source tables).

PURPOSE
    Strip plot of the global bivariate Moran's I per whole tissue (purple) and per TMA patient (beige,
    mean over cores), with the median per group and the across-unit test above each group. Negative =
    SecA-high cells sit among SecB-low neighbours. Binary 50 um weights, as Fig 4e.

INPUTS (output_root/singler_recompute/, from 04_fig4d_bivariate_moran.R)
    bivmoran_wt_SINGLER.csv, bivmoran_tma_patient_SINGLER.csv, bivmoran_tests_SINGLER.csv

OUTPUTS (output_root/figures/figure4/)
    Fig4d_bivariate_moran.{png,svg}
    Fig4d_bivariate_moran_wholetissue.csv  per tissue, study codes (source data)
    Fig4d_bivariate_moran_tests.csv        across-unit tests

MANUSCRIPT PANEL(S): Fig 4d
RUNTIME TIER: fast
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402
from config.study_codes import code_of  # noqa: E402

SG = path("output_root", "singler_recompute")
OUT = path("output_root", "figures", "figure4")
os.makedirs(OUT, exist_ok=True)

FA, FK, FN = 8, 8, 7
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "xtick.labelsize": FK, "ytick.labelsize": FK,
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 450,
})
WT_COLOR, TMA_COLOR = "#8A5DAF", "#D4A574"

wt = pd.read_csv(f"{SG}/bivmoran_wt_SINGLER.csv")
pt = pd.read_csv(f"{SG}/bivmoran_tma_patient_SINGLER.csv")
tests = pd.read_csv(f"{SG}/bivmoran_tests_SINGLER.csv").set_index("cohort")


def p_text(p):
    if p >= 0.001:
        return f"p = {p:.2g}" if p < 0.01 else f"p = {p:.2f}"
    e = int(np.floor(np.log10(p)))
    return f"p < 10$^{{{e + 1}}}$"


# ---- source data (study codes, code order) ----
src = wt.assign(study_code=[code_of("Whole tissue", t) for t in wt.tissue])
src = src[["study_code", "n_cells", "mean_neighbours", "I", "p_perm", "within_cell_r"]].sort_values("study_code")
src.rename(columns={"I": "bivariate_morans_I", "p_perm": "p_permutation_999",
                    "within_cell_r": "within_cell_pearson_r_SecA_SecB"}).round(5).to_csv(
    f"{OUT}/Fig4d_bivariate_moran_wholetissue.csv", index=False)
tt = tests.reset_index().round({"median_I": 6, "min_I": 6, "max_I": 6})
for c in ["sign_test_p", "wilcoxon_p"]:   # 3 significant digits (round(6) wrote the TMA p values as 0.0)
    tt[c] = tt[c].map(lambda v: float(f"{v:.3g}"))
tt.to_csv(f"{OUT}/Fig4d_bivariate_moran_tests.csv", index=False)

# ---- plot ----
fig, ax = plt.subplots(figsize=(2.4, 1.35))
ax.axhline(0, color="black", lw=0.4, alpha=0.3, zorder=0)
rng = np.random.default_rng(SEED)
tma_v, wt_v = pt.mean_I.values, wt.I.values
ax.scatter(2 + rng.uniform(-0.3, 0.3, len(tma_v)), tma_v, s=4, color=TMA_COLOR, alpha=0.5,
           edgecolors="none", zorder=2, rasterized=True)
ax.scatter(1 + rng.uniform(-0.15, 0.15, len(wt_v)), wt_v, s=12, color=WT_COLOR, alpha=0.8,
           edgecolors="white", linewidths=0.3, zorder=3)
for pos, vals, col in [(1, wt_v, WT_COLOR), (2, tma_v, TMA_COLOR)]:
    med = np.median(vals)
    ax.plot([pos - 0.35, pos + 0.35], [med, med], color=col, lw=1.0, zorder=4)
top = max(np.max(wt_v), np.max(tma_v))
w, t = tests.loc["whole tissue"], tests.loc["TMA"]
ax.text(1, top, f"{int(w.n_negative)}/{int(w.n_units)} < 0\nsign test {p_text(w.sign_test_p)}",
        ha="center", va="bottom", fontsize=6)
ax.text(2, top, f"{int(t.n_negative)}/{int(t.n_units)} < 0\nWilcoxon {p_text(t.wilcoxon_p)}",
        ha="center", va="bottom", fontsize=6)
ax.set_xticks([1, 2])
ax.set_xticklabels(["Whole\ntissue", "TMA\npatients"], fontsize=FK)
ax.set_xlim(0.4, 2.6)
ax.set_ylabel("Bivariate Moran's I\n(SecA vs SecB)", fontsize=FA)
ax.text(0.02, 0.02, "← segregated", fontsize=6, color="grey", style="italic", transform=ax.transAxes, va="bottom")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.4)
ax.tick_params(width=0.4, length=2, colors="black")
fig.tight_layout(pad=0.3)
fig.savefig(f"{OUT}/Fig4d_bivariate_moran.png", bbox_inches="tight", dpi=450)
fig.savefig(f"{OUT}/Fig4d_bivariate_moran.svg", bbox_inches="tight")
plt.close(fig)
print(src.to_string(index=False))
print(tests.to_string())
print("saved", f"{OUT}/Fig4d_bivariate_moran.{{png,svg}}")
