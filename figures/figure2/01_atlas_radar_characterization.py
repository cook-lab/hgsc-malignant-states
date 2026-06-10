#!/usr/bin/env python3
"""
Figure 2A — Functional radar plots, SecA vs SecB (PROGENy / Hallmark / DoRothEA / Flux)
=======================================================================================
PURPOSE
    Four radar panels comparing SecA vs SecB functional characterization using
    curated features aligned with the manuscript text. Emits a 2x2, a 1x4, and
    individual single-panel SVGs.

INPUTS
    - 21 functional characterization z-scores:
        output_root/04_functional/21_{progeny,hallmark,dorothea,flux}_zscored.csv
        (rows indexed by epitype: SecA / Intermediate / SecB / Ciliated)

OUTPUTS
    - figures_dir/figure2/atlas_radar_characterization.{svg,png}              (2x2)
    - figures_dir/figure2/atlas_radar_characterization_1x4.{svg,png}          (1x4)
    - figures_dir/figure2/atlas_radar_characterization_{progeny,hallmark,dorothea,flux}.svg

MANUSCRIPT PANEL(S): Fig 2A.

RUNTIME TIER: fast.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
DATA_DIR = path("data_root", "2026_final_atlas", "output", "21_epitype_functional_characterization")
OUT_STEM = "atlas_radar_characterization"

# ---------- Style ----------
FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FA,
    "axes.labelsize":    FA,
    "xtick.labelsize":   FK,
    "ytick.labelsize":   FK,
    "legend.fontsize":   FN,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
})

GROUPS = ["SecA", "SecB"]
DISPLAY = {"SecA": "SecA", "SecB": "SecB"}
PALETTE = {"SecA": "#E6A141", "SecB": "#9A7D55"}

PROGENY_FEATURES = ["MAPK", "p53", "TGFb", "TNFa", "VEGF", "Hypoxia"]
HALLMARK_FEATURES = ["E2F_TARGETS", "MYC_TARGETS_V1", "P53_PATHWAY",
                     "HYPOXIA", "PROTEIN_SECRETION", "TGF_BETA_SIGNALING"]
DOROTHEA_FEATURES = ["E2F4", "MYC", "HIF1A", "FOXO1", "TP53", "TEAD4"]
FLUX_FEATURES = [
    ("M_163", "Thymidylate\nsynthesis"), ("M_7", "Citrate\nsynthesis"),
    ("M_23", "Glutathione\nsynthesis"), ("M_5", "TCA entry"),
    ("M_169", "Steroid hormone\nsynthesis"), ("M_4", "Glycolysis"),
]

LABEL_CLEAN = {
    "E2F_TARGETS": "E2F targets", "MYC_TARGETS_V1": "MYC targets",
    "P53_PATHWAY": "p53 pathway", "HYPOXIA": "Hypoxia",
    "PROTEIN_SECRETION": "Protein secretion", "TGF_BETA_SIGNALING": "TGF-β signaling",
    "TGFb": "TGF-β", "TNFa": "TNF-α",
}


def clean_labels(features):
    return [LABEL_CLEAN.get(f, f) for f in features]


# ---------- Load ----------
prog_z = pd.read_csv(os.path.join(DATA_DIR, "21_progeny_zscored.csv"), index_col=0)
hall_z = pd.read_csv(os.path.join(DATA_DIR, "21_hallmark_zscored.csv"), index_col=0)
doro_z = pd.read_csv(os.path.join(DATA_DIR, "21_dorothea_zscored.csv"), index_col=0)
flux_z = pd.read_csv(os.path.join(DATA_DIR, "21_flux_zscored.csv"), index_col=0)

prog_z = prog_z.loc[GROUPS, PROGENY_FEATURES]
hall_z = hall_z.loc[GROUPS, HALLMARK_FEATURES]
doro_z = doro_z.loc[GROUPS, DOROTHEA_FEATURES]
flux_ids = [f[0] for f in FLUX_FEATURES]
flux_labels = [f[1] for f in FLUX_FEATURES]
flux_z = flux_z.loc[GROUPS, flux_ids].copy()
flux_z.columns = flux_labels


def radar_plot(ax, df_z, label_fontsize=None):
    if label_fontsize is None:
        label_fontsize = FN
    features = df_z.columns.tolist()
    labels = clean_labels(features)
    n = len(features)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles += angles[:1]
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    for grp in GROUPS:
        values = df_z.loc[grp].values.tolist()
        values += values[:1]
        ax.plot(angles, values, linewidth=1.2, color=PALETTE[grp])
        ax.fill(angles, values, alpha=0.18, color=PALETTE[grp])
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=label_fontsize, linespacing=0.9)
    ax.set_rlabel_position(30)
    ax.tick_params(axis="y", labelsize=label_fontsize - 0.5, pad=1)
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax.grid(linewidth=0.3, alpha=0.5)
    ax.spines["polar"].set_linewidth(0.3)


# ---------- 2x2 ----------
fig = plt.figure(figsize=(8.5, 11))
gs = gridspec.GridSpec(2, 2, figure=fig, left=0.08, right=0.92, bottom=0.06,
                       top=0.94, wspace=0.45, hspace=0.40)
for spec, df_z in [(gs[0, 0], prog_z), (gs[0, 1], hall_z), (gs[1, 0], doro_z), (gs[1, 1], flux_z)]:
    radar_plot(fig.add_subplot(spec, polar=True), df_z)
handles = [Line2D([0], [0], color=PALETTE[g], linewidth=1.5, label=DISPLAY[g]) for g in GROUPS]
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=FN,
           bbox_to_anchor=(0.5, 0.01))
fig.savefig(path("figures_dir", "figure2", f"{OUT_STEM}.svg"), format="svg")
fig.savefig(path("figures_dir", "figure2", f"{OUT_STEM}.png"), format="png")
plt.close(fig)

# ---------- 1x4 ----------
fig_h = plt.figure(figsize=(8.5, 3.5))
gs_h = gridspec.GridSpec(1, 4, figure=fig_h, left=0.04, right=0.96, bottom=0.10,
                         top=0.88, wspace=0.55)
FL = 7
for spec, df_z in [(gs_h[0, 0], prog_z), (gs_h[0, 1], hall_z), (gs_h[0, 2], doro_z), (gs_h[0, 3], flux_z)]:
    radar_plot(fig_h.add_subplot(spec, polar=True), df_z, label_fontsize=FL)
fig_h.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=FL,
             bbox_to_anchor=(0.5, 0.01))
fig_h.savefig(path("figures_dir", "figure2", f"{OUT_STEM}_1x4.svg"), format="svg")
fig_h.savefig(path("figures_dir", "figure2", f"{OUT_STEM}_1x4.png"), format="png")
plt.close(fig_h)

# ---------- Individual ----------
for df_z, stem in [(prog_z, f"{OUT_STEM}_progeny"), (hall_z, f"{OUT_STEM}_hallmark"),
                   (doro_z, f"{OUT_STEM}_dorothea"), (flux_z, f"{OUT_STEM}_flux")]:
    fig, ax = plt.subplots(figsize=(80 / 25.4, 80 / 25.4), subplot_kw=dict(polar=True))
    radar_plot(ax, df_z)
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.35, 1.1),
              frameon=False, fontsize=FN)
    fig.tight_layout()
    fig.savefig(path("figures_dir", "figure2", f"{stem}.svg"), format="svg")
    plt.close(fig)

print("Done.")
