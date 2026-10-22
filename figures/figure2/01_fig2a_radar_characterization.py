#!/usr/bin/env python3
"""
Figure 2a — Functional radar plots, SecA vs SecB (PROGENy / Hallmark / DoRothEA / Flux)
=======================================================================================
PURPOSE
    Four radar panels in one row comparing SecA vs SecB functional
    characterization, using curated features (six per resource) aligned with the
    manuscript text.

INPUTS
    - Per-state z-scored functional summaries rebuilt from the deposited per-cell
      decoupler / scFEA outputs joined to the final-atlas epithelial state
      (load_epithelial()['state']) by figures/_prep/fig2a_functional_state_summaries.py:
        output_root/figure2/fig2a_functional/21_{progeny,hallmark,dorothea,flux}_zscored.csv
        (rows: SecA / Intermediate / SecB / Ciliated; z across the four states, ddof=1)

OUTPUTS
    - figures_dir/figure2/atlas_radar_characterization_1x4.{svg,png}          (1x4 row)

MANUSCRIPT PANEL(S): Fig 2a.

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
DATA_DIR = path("output_root", "figure2", "fig2a_functional")  # figures/_prep/fig2a_functional_state_summaries.py
OUT_STEM = "atlas_radar_characterization"

# ---------- Style ----------
# Enlarged for 2A: legend/titles = 9pt (FA/FN); radar spoke + radial-tick
# (axis-type) labels = 7.5pt (FK).
FA, FK, FN = 9, 7.5, 9
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
for _n in ("progeny", "hallmark", "dorothea", "flux"):
    assert os.path.exists(os.path.join(DATA_DIR, f"21_{_n}_zscored.csv")), \
        f"run figures/_prep/fig2a_functional_state_summaries.py first ({DATA_DIR})"
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
        label_fontsize = FK  # radar spoke labels are axis-type labels -> 6pt
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
    ax.tick_params(axis="y", labelsize=FK, pad=1)  # radial ticks -> 6pt
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax.grid(linewidth=0.3, alpha=0.5)
    ax.spines["polar"].set_linewidth(0.3)


# ---------- Legend handles ----------
handles = [Line2D([0], [0], color=PALETTE[g], linewidth=1.5, label=DISPLAY[g]) for g in GROUPS]

# ---------- 1x4 (Fig 2a: letter-page footprint 7.6 x 1.95, wide-and-short) ----------
fig_h = plt.figure(figsize=(7.6, 2.05))
gs_h = gridspec.GridSpec(1, 4, figure=fig_h, left=0.02, right=0.87, bottom=0.14,
                         top=0.88, wspace=1.05)
for spec, df_z in [(gs_h[0, 0], prog_z), (gs_h[0, 1], hall_z), (gs_h[0, 2], doro_z), (gs_h[0, 3], flux_z)]:
    radar_plot(fig_h.add_subplot(spec, polar=True), df_z, label_fontsize=FK)
fig_h.legend(handles=handles, loc="center right", ncol=1, frameon=False, fontsize=FN,
             bbox_to_anchor=(0.995, 0.5))
fig_h.savefig(path("figures_dir", "figure2", f"{OUT_STEM}_1x4.svg"), format="svg")
fig_h.savefig(path("figures_dir", "figure2", f"{OUT_STEM}_1x4.png"), format="png")
plt.close(fig_h)

print("Done.")
