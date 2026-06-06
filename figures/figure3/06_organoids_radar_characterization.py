#!/usr/bin/env python3
"""
Figure 3G — Organoid functional radar plots, SecA vs SecB (atlas-matched features)
==================================================================================
PURPOSE
    Three radar panels (PROGENy / Hallmark / DoRothEA) using the EXACT atlas
    curated feature sets, for organoid SecA (SecB-low) vs SecB (SecB-high).
    Emits a 1x3 and individual single-panel SVGs. Flux omitted (no 09b flux).

INPUTS
    - <organoids_root>/output/09_organoid_secB_characterization/
        09b_{progeny,hallmark,dorothea}_zscored.csv  (rows SecB-low / SecB-high)
        EXTERNAL DEPENDENCY (override ORGANOIDS_ROOT).

OUTPUTS
    - figures_dir/organoids_radar_characterization_secAB.{svg,png}        (1x3)
    - figures_dir/organoids_radar_characterization_secAB_{progeny,hallmark,dorothea}.svg

MANUSCRIPT PANEL(S): Fig 3G.

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
from config.config import obj, path  # noqa: E402

# ---------- Paths ----------
DATA_DIR = os.path.join(obj("organoids_root"),
                        "output/09_organoid_secB_characterization")
OUT_STEM = "organoids_radar_characterization_secAB"

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

PLOT_GROUPS = ["SecA", "SecB"]
DISPLAY = {"SecA": "SecA", "SecB": "SecB"}
PALETTE = {"SecA": "#E6A141", "SecB": "#9A7D55"}
RENAME = {"SecB-low": "SecA", "SecB-high": "SecB"}

PROGENY_FEATURES = ["MAPK", "p53", "TGFb", "TNFa", "VEGF", "Hypoxia"]
HALLMARK_FEATURES = ["E2F_TARGETS", "MYC_TARGETS_V1", "P53_PATHWAY",
                     "HYPOXIA", "PROTEIN_SECRETION", "TGF_BETA_SIGNALING"]
DOROTHEA_FEATURES = ["E2F4", "MYC", "HIF1A", "FOXO1", "TP53", "TEAD4"]

LABEL_CLEAN = {
    "E2F_TARGETS": "E2F targets", "MYC_TARGETS_V1": "MYC targets",
    "P53_PATHWAY": "p53 pathway", "HYPOXIA": "Hypoxia",
    "PROTEIN_SECRETION": "Protein secretion", "TGF_BETA_SIGNALING": "TGF-β signaling",
    "TGFb": "TGF-β", "TNFa": "TNF-α",
}


def clean_labels(features):
    return [LABEL_CLEAN.get(f, f) for f in features]


def load_subset(csv_name, features):
    df = pd.read_csv(os.path.join(DATA_DIR, csv_name), index_col=0)
    df = df.loc[["SecB-low", "SecB-high"], features]
    df.index = df.index.map(RENAME)
    return df


prog_z = load_subset("09b_progeny_zscored.csv", PROGENY_FEATURES)
hall_z = load_subset("09b_hallmark_zscored.csv", HALLMARK_FEATURES)
doro_z = load_subset("09b_dorothea_zscored.csv", DOROTHEA_FEATURES)


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
    for grp in PLOT_GROUPS:
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


handles = [Line2D([0], [0], color=PALETTE[g], linewidth=1.5, label=DISPLAY[g]) for g in PLOT_GROUPS]

# ---------- 1x3 ----------
fig = plt.figure(figsize=(8.5, 3.5))
gs = gridspec.GridSpec(1, 3, figure=fig, left=0.04, right=0.96, bottom=0.10, top=0.88, wspace=0.55)
FL = 7
for spec, df_z in [(gs[0, 0], prog_z), (gs[0, 1], hall_z), (gs[0, 2], doro_z)]:
    radar_plot(fig.add_subplot(spec, polar=True), df_z, label_fontsize=FL)
fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=FL,
           bbox_to_anchor=(0.5, 0.01))
fig.savefig(path("figures_dir", f"{OUT_STEM}.svg"), format="svg")
fig.savefig(path("figures_dir", f"{OUT_STEM}.png"), format="png")
plt.close(fig)

# ---------- Individual ----------
for df_z, stem in [(prog_z, f"{OUT_STEM}_progeny"), (hall_z, f"{OUT_STEM}_hallmark"),
                   (doro_z, f"{OUT_STEM}_dorothea")]:
    fig, ax = plt.subplots(figsize=(80 / 25.4, 80 / 25.4), subplot_kw=dict(polar=True))
    radar_plot(ax, df_z)
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(1.35, 1.1),
              frameon=False, fontsize=FN)
    fig.tight_layout()
    fig.savefig(path("figures_dir", f"{stem}.svg"), format="svg")
    plt.close(fig)

print("Done.")
