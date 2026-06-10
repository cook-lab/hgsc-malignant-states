#!/usr/bin/env python3
"""
Figure 1B-E,I — Study-contribution 5-panel column with epithelial composition
=============================================================================
PURPOSE
    Vertical 5-panel column under the Fig 1A UMAP: cells per study (log),
    metastatic-site %, treatment-status %, cell-type % (level1), and epithelial
    composition % (Ciliated/SecA/Intermediate/SecB). Manuscript-ready (no titles
    / panel letters). Studies ordered descending by total cell count.

INPUTS
    - fig_data_dir/meta.parquet     (schema_nmf 4-class per epithelial cell;
        produced by figures/_prep/fig_secretory_polarization_00_prepare_data.py)
    - fig_data_fig1/meta.parquet    (whole-atlas obs: study/site/treatment)
    - fig_data_fig1/panel_b_cells_per_study.csv
    - fig_data_fig1/panel_b_patients_per_study.csv
    - fig_data_fig1/panel_g_composition_by_study.csv

OUTPUTS
    - figures_dir/figure1/atlas_study_contribution_with_epi_composition.{svg,png}

MANUSCRIPT PANEL(S): Fig 1B, 1C, 1D, 1E, 1I.

RUNTIME TIER: fast.

NOTE: epithelial label standardized "Transitioning" -> "Intermediate".
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
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
SCHEMA_PQ = path("data_root", "2026_final_atlas", "output", "fig_secretory_polarization", "data", "meta.parquet")
STUDY_META = path("data_root", "202605_epitype_manuscript", "final_publication_figures", "data_fig1", "meta.parquet")
CELLS_CSV = path("data_root", "202605_epitype_manuscript", "final_publication_figures", "data_fig1", "panel_b_cells_per_study.csv")
PATIENTS_CSV = path("data_root", "202605_epitype_manuscript", "final_publication_figures", "data_fig1", "panel_b_patients_per_study.csv")
COMP_CSV = path("data_root", "202605_epitype_manuscript", "final_publication_figures", "data_fig1", "panel_g_composition_by_study.csv")

OUT_SVG = path("figures_dir", "figure1", "atlas_study_contribution_with_epi_composition.svg")
OUT_PNG = path("figures_dir", "figure1", "atlas_study_contribution_with_epi_composition.png")

for p in (SCHEMA_PQ, STUDY_META, CELLS_CSV, PATIENTS_CSV, COMP_CSV):
    assert os.path.exists(p), f"missing input: {p}"

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
    "regner_2021": "Regner 2021", "vazquez_garcia_2022": "Vazquez-Garcia 2022",
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

FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        FA,
    "axes.labelsize":   FA,
    "axes.titlesize":   0,
    "xtick.labelsize":  FK,
    "ytick.labelsize":  FN,
    "legend.fontsize":  FN,
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,
    "ps.fonttype":      42,
})

# ---------- Load inputs ----------
cells_per_study = pd.read_csv(CELLS_CSV, index_col=0).squeeze()
patients_per_study = pd.read_csv(PATIENTS_CSV, index_col=0).squeeze()
comp = pd.read_csv(COMP_CSV, index_col=0)

schema = pd.read_parquet(SCHEMA_PQ, columns=["schema_nmf"])
study_meta = pd.read_parquet(STUDY_META, columns=["study", "metastatic_site", "treatment_status"])

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

df = schema.join(study_meta[["study"]], how="left")
df = df.dropna(subset=["study"]).copy()
counts = df.groupby(["study", "schema_nmf"], observed=True).size().unstack(fill_value=0)
counts = counts.reindex(columns=EPI_ORDER, fill_value=0)
epi_prop = counts.div(counts.sum(axis=1), axis=0) * 100

study_order = (
    cells_per_study.reindex(STUDY_DISPLAY.keys())
    .sort_values(ascending=False).index.tolist()
)
n_studies = len(study_order)
x_pos = np.arange(n_studies)
x_labels = [f"{STUDY_DISPLAY[s]} (n={int(patients_per_study[s])})" for s in study_order]

# ---------- Figure ----------
fig = plt.figure(figsize=(5.4, 9.0))
gs = gridspec.GridSpec(
    5, 1, figure=fig, height_ratios=[1.4, 0.9, 0.9, 1.4, 0.9], hspace=0.18,
    left=0.10, right=0.62, top=0.78, bottom=0.05,
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
ax1.set_xticklabels(x_labels, rotation=90, ha="center", va="bottom", fontsize=FN)
ax1.tick_params(axis="x", which="both", length=0, pad=2)

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

# Panel 3 — treatment status
bottom_acc = np.zeros(n_studies)
for tx in TREATMENT_ORDER:
    vals = np.array([tx_prop.loc[s, tx] for s in study_order], dtype=float)
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
                    for t in TREATMENT_ORDER if tx_counts[t].sum() > 0], ncol=1, **LEGEND_KW)
ax4.legend(handles=[Patch(facecolor=LEVEL1_PALETTE[ct], edgecolor="none", label=ct)
                    for ct in LEVEL1_ORDER], ncol=2, **LEGEND_KW)
ax5.legend(handles=[Patch(facecolor=EPI_PALETTE[cls], edgecolor="white", linewidth=0.3,
                          label=cls) for cls in EPI_ORDER], ncol=1, **LEGEND_KW)

fig.savefig(OUT_SVG, format="svg", bbox_inches="tight", pad_inches=0.02)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print(f"[save] {OUT_SVG}\n[save] {OUT_PNG}\nDone.")
