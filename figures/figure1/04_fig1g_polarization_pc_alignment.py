#!/usr/bin/env python3
"""
Figure 1g — Polarization axis is the leading coordinated axis of variation
==========================================================================
PURPOSE
    Per tumour, the rank of the principal component of malignant cell-state
    variation that is most correlated with the SecA<->SecB secretory-polarization
    score (PC1 in 54%; top-3 in 82% of 250 tumours). Histogram of that rank
    across atlas tumours (ranks 1-10; >10 clipped to 10), dashed marker at PC1.
    Manuscript-ready.

INPUTS
    - Per-sample PC alignment computed from the final atlas by
      figures/_prep/fig1g_polarization_pc_alignment.py (malignant set =
      CopyKAT-aneuploid epithelial cells of the final atlas, counts from the
      epithelial subset h5ad, >= 150 cells per sample; polarization = organoid
      SecB - SecA score):
        output_root/_prep_caches/fig1g_polarization_pc_alignment.csv
      Column `pol_top_pc` = rank of the PC most correlated with the polarization
      score for that tumour.

OUTPUTS
    - figures_dir/figure1/atlas_polarization_pc_alignment.{svg,png}
    - output_root/figure1/atlas_polarization_pc_alignment_data.csv (rank->n)
    - printed: PC1 % and top-3 % of tumours (Results text)

MANUSCRIPT PANEL(S): Fig 1g.

RUNTIME TIER: fast (reads one CSV; no h5ad). Run the _prep script first (moderate).
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

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
RESULTS_CSV = path("output_root", "_prep_caches", "fig1g_polarization_pc_alignment.csv")
OUT_SVG = path("figures_dir", "figure1", "atlas_polarization_pc_alignment.svg")
OUT_PNG = path("figures_dir", "figure1", "atlas_polarization_pc_alignment.png")
OUT_CSV = path("output_root", "figure1", "atlas_polarization_pc_alignment_data.csv")

assert os.path.exists(RESULTS_CSV), (
    f"Missing: {RESULTS_CSV} (run figures/_prep/fig1g_polarization_pc_alignment.py first)")

# SecB colour (the shared epithelial-state palette).
SECB_COLOR = "#9A7D55"

# ---------- Style ----------
FA, FK = 8, 6   # standardized 8/6; kept tight for a small (~1.55 in) square panel
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        FA,
    "axes.labelsize":   FA,
    "axes.titlesize":   0,
    "xtick.labelsize":  FK,
    "ytick.labelsize":  FK,
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,
    "figure.dpi":       450,
    "savefig.dpi":      450,
})

# ---------- Load & bin ----------
fr = pd.read_csv(RESULTS_CSV)
top = fr["pol_top_pc"].dropna().astype(int)
n_tumours = int(top.shape[0])

pc1_pct = round(100 * float((top == 1).mean()))
top3_pct = round(100 * float((top <= 3).mean()))

rank = np.minimum(top, 10)                                  # >10 clipped to 10
counts = rank.value_counts().reindex(range(1, 11), fill_value=0).sort_index()

pd.DataFrame({"rank": counts.index, "n": counts.values}).to_csv(OUT_CSV, index=False)
print(f"PC1 {pc1_pct}%, top-3 {top3_pct}% (n={n_tumours} tumours)")

# ---------- Plot ----------
# Horizontal orientation: rank on Y running top->bottom
# (1 at top), Tumours (n) on X extending rightward.
fig, ax = plt.subplots(figsize=(1.30, 1.30))
ax.barh(counts.index, counts.values, height=1.0, color=SECB_COLOR,
        edgecolor="white", linewidth=0.2, align="center")
ax.axhline(1, linestyle="--", linewidth=0.4, color="#404040")
ax.text(0.97, 0.45, f"PC1: {pc1_pct}%", transform=ax.transAxes,
        ha="right", va="center", fontsize=FA, color="#4D4D4D")

ax.set_yticks([1, 2, 4, 6, 8, 10])
ax.set_yticklabels(["1", "2", "4", "6", "8", "10+"])
ax.set_ylim(0.4, 10.6)
ax.invert_yaxis()                          # rank 1 at top
ax.set_xlim(0, counts.values.max() * 1.12)
ax.set_ylabel("Rank of polarization PC", fontsize=FA, labelpad=3, linespacing=0.9)
ax.set_xlabel("Tumours (n)", fontsize=FA, labelpad=3)
ax.tick_params(axis="both", length=2, width=0.4, pad=1)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["left"].set_linewidth(0.5)

fig.savefig(OUT_SVG, format="svg", bbox_inches="tight", pad_inches=0.04)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}\n  Saved: {OUT_CSV}")
