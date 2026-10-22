#!/usr/bin/env python3
"""
Figure 1b-e, k — Study-contribution 5-panel column with epithelial composition
==============================================================================
PURPOSE
    Vertical 5-panel column under the Fig 1a UMAP: cells per study (log; b),
    metastatic-site % (c), treatment-status % (d), cell-type % (level 1; e), and
    epithelial state composition % (Ciliated/SecA/Intermediate/SecB; k).
    Manuscript-ready (no titles / panel letters). Studies ordered descending by
    total cell count.

INPUTS
    All values are computed here from the final atlas via the shared loader
    figures/_prep/atlas_final_tables.py (caches in output_root/_prep_caches/):
    - load_obs()        (1,980,703 cells: study, patient_id, metastatic_site,
        treatment_status, celltype_level1) -> cells/study, patients/study,
        site %, treatment %, level-1 cell-type %.
    - load_epithelial() (575,366 epithelial cells: study, state
        SecA/Intermediate/SecB/Ciliated) -> epithelial state % per study.

OUTPUTS
    - figures_dir/figure1/atlas_study_contribution_with_epi_composition.{svg,png}
    - output_root/figure1/atlas_study_contribution_values.csv (plotted values, one row per study)
    - printed: cells, studies and patients of the atlas (Results text)

MANUSCRIPT PANEL(S): Fig 1b, 1c, 1d, 1e, 1k.

RUNTIME TIER: fast.

NOTE: epithelial labels come from the loader's `state` column
    (SecA / Intermediate / SecB / Ciliated).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

from figures._prep.atlas_final_tables import load_obs, load_epithelial  # noqa: E402

# ---------- Paths ----------
OUT_SVG = path("figures_dir", "figure1", "atlas_study_contribution_with_epi_composition.svg")
OUT_PNG = path("figures_dir", "figure1", "atlas_study_contribution_with_epi_composition.png")
OUT_CSV = path("output_root", "figure1", "atlas_study_contribution_values.csv")

# ---------- Constants ----------
LEVEL1_ORDER = [
    "Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte",
    "Endothelial", "T/NK cell", "B cell", "Plasma cell", "Macrophage",
    "DC", "Neutrophil", "Mast cell",
]
LEVEL1_PALETTE = {
    "Epithelial":     "#E6A141",
    "Mesothelial":    "#D4A574",
    "Fibroblast":     "#C4B9A8",
    "Smooth muscle":  "#D14E6C",
    "Pericyte":       "#B87A7A",
    "Endothelial":    "#7D4E4E",
    "T/NK cell":      "#87CEFA",
    "B cell":         "#5665B6",
    "Plasma cell":    "#8A5DAF",
    "Macrophage":     "#8FBC8F",
    "DC":             "#2E8B57",
    "Neutrophil":     "#6B8E23",
    "Mast cell":      "#8B9B6B",
}
STUDY_DISPLAY = {
    "denisenko_2024": "Denisenko 2024", "geistlinger_2020": "Geistlinger 2020",
    "hornburg_2021": "Hornburg 2021", "loret_2022": "Loret 2022",
    "luo_2024": "Luo 2024", "nath_2021": "Nath 2021",
    "olalekan_2021": "Olalekan 2021", "olbrecht_2021": "Olbrecht 2021",
    "regner_2021": "Regner 2021", "vazquez_garcia_2022": "Vázquez-García 2022",
    "xu_2022": "Xu 2022", "zhang_2022": "Zhang 2022", "zheng_2023": "Zheng 2023",
}

EPI_ORDER = ["Ciliated", "SecA", "Intermediate", "SecB"]
EPI_PALETTE = {
    "Ciliated":     "#E07850",
    "SecA":         "#E6A141",
    "Intermediate": "#C08E48",
    "SecB":         "#9A7D55",
}

SITE_ORDER = ["primary", "metastasis", "ascites", "healthy"]
SITE_DISPLAY = {"primary": "Primary", "metastasis": "Metastasis",
                "ascites": "Ascites", "healthy": "Healthy"}
SITE_PALETTE = {"primary": "#7A9EBF", "metastasis": "#B07AA1",
                "ascites": "#8FAC8C", "healthy": "#C2956B"}

TREATMENT_ORDER = [
    "pre-treatment", "post-chemotherapy", "post-chemotherapy_niraparib",
    "post-chemotherapy_olaparib", "post-chemotherapy_pembro",
    "post-niraparib", "NA",
]
TREATMENT_DISPLAY = {
    "pre-treatment": "Pre-treatment", "post-chemotherapy": "Post-chemo",
    "post-chemotherapy_niraparib": "Post-chemo + niraparib",
    "post-chemotherapy_olaparib": "Post-chemo + olaparib",
    "post-chemotherapy_pembro": "Post-chemo + pembro",
    "post-niraparib": "Post-niraparib", "NA": "Unknown",
}
TREATMENT_PALETTE = {
    "pre-treatment": "#DDD5CA", "post-chemotherapy": "#009E73",
    "post-chemotherapy_niraparib": "#E69F00", "post-chemotherapy_olaparib": "#0072B2",
    "post-chemotherapy_pembro": "#D55E00", "post-niraparib": "#CC79A7", "NA": "#666666",
}
# Display grouping (Fig 1d): lump the three post-chemo maintenance combos into a
# single "Other" stacked segment (sum of their per-study %). Data values unchanged;
# Pre-treatment / Post-chemo / Post-niraparib / Unknown kept as-is.
TREATMENT_LUMP = ["post-chemotherapy_niraparib", "post-chemotherapy_olaparib", "post-chemotherapy_pembro"]
TREATMENT_ORDER_DISPLAY = ["pre-treatment", "post-chemotherapy", "Other", "post-niraparib", "NA"]
TREATMENT_DISPLAY["Other"] = "Other"
TREATMENT_PALETTE["Other"] = "#9B8AA6"  # muted mauve-grey; distinct from Unknown (#666666)

FA, FK, FN = 8, 6, 8
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        FA,
    "axes.labelsize":   FA,
    "axes.titlesize":   0,
    "xtick.labelsize":  FK,
    "ytick.labelsize":  FK,
    "legend.fontsize":  FN,
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
})

# ---------- Load inputs (final atlas) ----------
study_meta = load_obs(columns=["study", "patient_id", "metastatic_site", "treatment_status", "celltype_level1"])
epi_meta = load_epithelial(columns=["study", "state"])
assert len(study_meta) == 1_980_703 and len(epi_meta) == 575_366, (len(study_meta), len(epi_meta))

cells_per_study = study_meta.groupby("study", observed=True).size()
patients_per_study = study_meta.groupby("study", observed=True)["patient_id"].nunique()
comp = pd.crosstab(study_meta["study"], study_meta["celltype_level1"])
comp = comp.div(comp.sum(axis=1), axis=0) * 100
assert set(cells_per_study.index) == set(STUDY_DISPLAY), set(cells_per_study.index) ^ set(STUDY_DISPLAY)
print(f"[atlas] {len(study_meta):,} cells, {cells_per_study.size} studies, "
      f"{study_meta['patient_id'].nunique()} patients; {len(epi_meta):,} epithelial cells")

site_counts = (
    study_meta.groupby(["study", "metastatic_site"], observed=True)
    .size().unstack(fill_value=0).reindex(columns=SITE_ORDER, fill_value=0)
)
site_prop = site_counts.div(site_counts.sum(axis=1), axis=0) * 100

tx_counts = (
    study_meta.groupby(["study", "treatment_status"], observed=True)
    .size().unstack(fill_value=0).reindex(columns=TREATMENT_ORDER, fill_value=0)
)
tx_prop = tx_counts.div(tx_counts.sum(axis=1), axis=0) * 100

# Lump the three maintenance combos into "Other" for display (same denominator, so
# Other% == sum of the three component %; totals per study preserved at 100%).
tx_counts_disp = tx_counts.copy()
tx_counts_disp["Other"] = tx_counts_disp[TREATMENT_LUMP].sum(axis=1)
tx_counts_disp = tx_counts_disp.drop(columns=TREATMENT_LUMP).reindex(columns=TREATMENT_ORDER_DISPLAY, fill_value=0)
tx_prop_disp = tx_counts_disp.div(tx_counts_disp.sum(axis=1), axis=0) * 100
_verify = [s for s in tx_counts_disp.index if tx_counts_disp.loc[s, "Other"] > 0][:3] or list(tx_counts_disp.index[:3])
for _s in _verify:
    _comp = float(tx_prop.loc[_s, TREATMENT_LUMP].sum())
    _other = float(tx_prop_disp.loc[_s, "Other"])
    print(f"[verify treatment lump] {_s}: Other={_other:.4f}%  sum(3 components)={_comp:.4f}%  match={np.isclose(_other, _comp)}")

counts = epi_meta.groupby(["study", "state"], observed=True).size().unstack(fill_value=0)
assert set(counts.columns) <= set(EPI_ORDER), f"unexpected epithelial labels: {set(counts.columns) - set(EPI_ORDER)}"
counts = counts.reindex(columns=EPI_ORDER, fill_value=0)
epi_prop = counts.div(counts.sum(axis=1), axis=0) * 100

# Plotted values, one row per study (for the report / caption numbers).
pd.concat({
    "n": pd.DataFrame({"cells": cells_per_study, "patients": patients_per_study}),
    "site_pct": site_prop, "treatment_pct": tx_prop_disp,
    "celltype_pct": comp.reindex(columns=LEVEL1_ORDER, fill_value=0), "epithelial_pct": epi_prop,
}, axis=1).to_csv(OUT_CSV)
print(f"[save] {OUT_CSV}")

study_order = (
    cells_per_study.reindex(STUDY_DISPLAY.keys())
    .sort_values(ascending=False).index.tolist()
)
n_studies = len(study_order)
x_pos = np.arange(n_studies)
# Panel b labels are split: study NAME (no count) rotated above each bar (x-tick label),
# and the patient count "n=<count>" rendered white at the bar base (value-at-base style).
x_labels_full = [f"{STUDY_DISPLAY[s]} (n={int(patients_per_study[s])})" for s in study_order]
x_labels = [STUDY_DISPLAY[s] for s in study_order]
# Parse the count out of the existing "(n=…)" labels so the numbers are unchanged.
n_labels = [f"n={re.search(r'\(n=(\d+)\)', lbl).group(1)}" for lbl in x_labels_full]

# ---------- Figure ----------
fig = plt.figure(figsize=(4.50, 8.30))
gs = gridspec.GridSpec(
    5, 1, figure=fig, height_ratios=[1.35, 1.05, 1.05, 1.35, 1.05], hspace=0.24,
    left=0.10, right=0.62, top=0.82, bottom=0.04,
)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)
ax4 = fig.add_subplot(gs[3, 0], sharex=ax1)
ax5 = fig.add_subplot(gs[4, 0], sharex=ax1)
ALL_AXES = [ax1, ax2, ax3, ax4, ax5]

BAR_W = 0.78
EDGE_LW = 0.3
TICK_PARAMS = dict(length=2, width=0.4, pad=1)

# Panel 1 — cells per study
counts_sorted = np.array([cells_per_study[s] for s in study_order], dtype=float)
ax1.bar(x_pos, counts_sorted, width=BAR_W, color="#7A7A7A",
        edgecolor="#333333", linewidth=EDGE_LW)
ax1.set_yscale("log")
ax1.set_ylim(max(1.0, counts_sorted.min() * 0.7), counts_sorted.max() * 1.4)
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(
    lambda v, _: f"{v / 1e3:.0f}" if v >= 1000 else f"{v:.0f}"))
ax1.set_ylabel("Cells (×1,000)", fontsize=FA, labelpad=2)
ax1.tick_params(axis="y", which="both", **TICK_PARAMS)
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position("top")
ax1.set_xticks(x_pos)
ax1.set_xticklabels(x_labels, rotation=90, ha="center", va="bottom", fontsize=FK)
ax1.tick_params(axis="x", which="both", length=0, pad=2)
# Patient count as a vertical label per bar (centred), numbers parsed from the
# source labels. Default: WHITE inside the bar, sitting just above the x-axis line.
# For a bar too SHORT (on the log axis) to hold the label legibly, fall back to a
# DARK-GREY label placed just ABOVE the bar top so it reads on the white background.
# Threshold is general (compares each bar's measured top to the label's height), not
# hardcoded to any study.
_LAB_BASE = 0.02          # axes-fraction y where the white label's base sits
_LAB_MARGIN = 0.02        # headroom required above the label to keep it "inside"
fig.canvas.draw()         # realise a renderer so label heights can be measured
_rend = fig.canvas.get_renderer()
_ax_h_px = ax1.get_window_extent().height
_y0, _y1 = ax1.get_ylim()
_log_span = np.log10(_y1) - np.log10(_y0)
for xi, nlab, val in zip(x_pos, n_labels, counts_sorted):
    top_frac = (np.log10(val) - np.log10(_y0)) / _log_span   # bar top in axes fraction
    t = ax1.text(xi, _LAB_BASE, nlab, transform=ax1.get_xaxis_transform(),
                 rotation=90, ha="center", va="bottom", fontsize=FK)
    lab_frac = t.get_window_extent(renderer=_rend).height / _ax_h_px
    if top_frac >= _LAB_BASE + lab_frac + _LAB_MARGIN:
        t.set_color("white")                                 # fits inside the bar
    else:
        t.remove()                                           # too short -> dark, above top
        ax1.text(xi, top_frac + 0.01, nlab, transform=ax1.get_xaxis_transform(),
                 rotation=90, ha="center", va="bottom", fontsize=FK, color="#333333")

# Panel 2 — metastatic site
bottom_acc = np.zeros(n_studies)
for site in SITE_ORDER:
    vals = np.array([site_prop.loc[s, site] for s in study_order], dtype=float)
    if not vals.any():
        continue
    ax2.bar(x_pos, vals, bottom=bottom_acc, width=BAR_W,
            color=SITE_PALETTE[site], edgecolor="white", linewidth=EDGE_LW)
    bottom_acc = bottom_acc + vals
ax2.set_ylim(0, 100); ax2.set_yticks([0, 50, 100])
ax2.set_ylabel("Site (%)", fontsize=FA, labelpad=2)
ax2.tick_params(axis="y", which="both", **TICK_PARAMS)

# Panel 3 — treatment status (maintenance combos lumped into "Other")
bottom_acc = np.zeros(n_studies)
for tx in TREATMENT_ORDER_DISPLAY:
    vals = np.array([tx_prop_disp.loc[s, tx] for s in study_order], dtype=float)
    if not vals.any():
        continue
    ax3.bar(x_pos, vals, bottom=bottom_acc, width=BAR_W,
            color=TREATMENT_PALETTE[tx], edgecolor="white", linewidth=EDGE_LW)
    bottom_acc = bottom_acc + vals
ax3.set_ylim(0, 100); ax3.set_yticks([0, 50, 100])
ax3.set_ylabel("Treatment (%)", fontsize=FA, labelpad=2)
ax3.tick_params(axis="y", which="both", **TICK_PARAMS)

# Panel 4 — cell-type (level1)
bottom_acc = np.zeros(n_studies)
for ct in LEVEL1_ORDER:
    if ct not in comp.columns:
        continue
    vals = np.array([comp.loc[s, ct] for s in study_order], dtype=float)
    ax4.bar(x_pos, vals, bottom=bottom_acc, width=BAR_W,
            color=LEVEL1_PALETTE[ct], edgecolor="none", linewidth=0)
    bottom_acc = bottom_acc + vals
ax4.set_ylim(0, 100); ax4.set_yticks([0, 25, 50, 75, 100])
ax4.set_ylabel("Cell-type (%)", fontsize=FA, labelpad=2)
ax4.tick_params(axis="y", which="both", **TICK_PARAMS)

# Panel 5 — epithelial composition
bottom_acc = np.zeros(n_studies)
for cls in EPI_ORDER:
    vals = np.array([epi_prop.loc[s, cls] for s in study_order], dtype=float)
    ax5.bar(x_pos, vals, bottom=bottom_acc, width=BAR_W,
            color=EPI_PALETTE[cls], edgecolor="white", linewidth=EDGE_LW)
    bottom_acc = bottom_acc + vals
ax5.set_ylim(0, 100); ax5.set_yticks([0, 25, 50, 75, 100])
ax5.set_ylabel("Epithelial (%)", fontsize=FA, labelpad=2)
ax5.tick_params(axis="y", which="both", **TICK_PARAMS)

for ax in (ax2, ax3, ax4, ax5):
    ax.tick_params(axis="x", which="both", length=0, labelbottom=False, labeltop=False)
ax5.set_xlim(-0.6, n_studies - 0.4)

for ax in ALL_AXES:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.spines["left"].set_linewidth(0.5)

LEGEND_KW = dict(
    loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=FN,
    handlelength=1.0, handleheight=0.9, handletextpad=0.4, columnspacing=0.7,
    borderaxespad=0.0,
)
ax2.legend(handles=[Patch(facecolor=SITE_PALETTE[s], edgecolor="white", linewidth=0.3,
                          label=SITE_DISPLAY[s])
                    for s in SITE_ORDER if site_counts[s].sum() > 0], ncol=1, **LEGEND_KW)
ax3.legend(handles=[Patch(facecolor=TREATMENT_PALETTE[t], edgecolor="white", linewidth=0.3,
                          label=TREATMENT_DISPLAY[t])
                    for t in TREATMENT_ORDER_DISPLAY if tx_counts_disp[t].sum() > 0], ncol=1, **LEGEND_KW)
ax4.legend(handles=[Patch(facecolor=LEVEL1_PALETTE[ct], edgecolor="none", label=ct)
                    for ct in LEVEL1_ORDER], ncol=2, **LEGEND_KW)  # 2 cols: 1-col overlapped the treatment legend
ax5.legend(handles=[Patch(facecolor=EPI_PALETTE[cls], edgecolor="white", linewidth=0.3,
                          label=cls) for cls in EPI_ORDER], ncol=1, **LEGEND_KW)

fig.savefig(OUT_SVG, format="svg", bbox_inches="tight", pad_inches=0.02)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print(f"[save] {OUT_SVG}\n[save] {OUT_PNG}\nDone.")
