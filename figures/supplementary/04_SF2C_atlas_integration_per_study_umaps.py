#!/usr/bin/env python3
"""
SF2C — Per-study split UMAPs
============================

Purpose
    Small multiples: one UMAP per study, highlighting that study's cells in
    colour against a grey background. Shows each study contributes to the major
    clusters (no study-specific islands).

INPUTS
    output_root/figures/data/atlas_final_umap.parquet
        (columns used: UMAP1, UMAP2, study; from 00b_extract_integration_umaps.py)

OUTPUTS
    output_root/figures/supplementary/SF2C_atlas_integration_per_study_umaps.{svg,png}

MANUSCRIPT PANEL(S)
    SF2C.

RUNTIME TIER
    moderate (rasterized background + per-study scatter, subsampled).
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

DATA_PQ = path("output_root", "figures", "data", "atlas_final_umap.parquet")
OUT_SVG = path("output_root", "figures", "supplementary", "SF2C_atlas_integration_per_study_umaps.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF2C_atlas_integration_per_study_umaps.png")

# ============================================================================
# STYLE
# ============================================================================

FA, FK, FN = 6, 5.5, 5

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FA,
    "axes.labelsize":    FA,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
})

STUDY_PALETTE = {
    "denisenko_2024":      "#875692",
    "geistlinger_2020":    "#F38400",
    "hornburg_2021":       "#A1CAF1",
    "loret_2022":          "#BE0032",
    "luo_2024":            "#C2B280",
    "nath_2021":           "#848482",
    "olalekan_2021":       "#008856",
    "olbrecht_2021":       "#E68FAC",
    "regner_2021":         "#0067A5",
    "vazquez_garcia_2022": "#F99379",
    "xu_2022":             "#604E97",
    "zhang_2022":          "#F6A600",
    "zheng_2023":          "#B3446C",
}
STUDY_DISPLAY = {
    "denisenko_2024": "Denisenko 2024", "geistlinger_2020": "Geistlinger 2020",
    "hornburg_2021": "Hornburg 2021", "loret_2022": "Loret 2022", "luo_2024": "Luo 2024",
    "nath_2021": "Nath 2021", "olalekan_2021": "Olalekan 2021",
    "olbrecht_2021": "Olbrecht 2021", "regner_2021": "Regner 2021",
    "vazquez_garcia_2022": "Vazquez-Garcia 2022", "xu_2022": "Xu 2022",
    "zhang_2022": "Zhang 2022", "zheng_2023": "Zheng 2023",
}
STUDY_ORDER = list(STUDY_PALETTE.keys())

# ============================================================================
# DATA
# ============================================================================

print("Loading data...", flush=True)
df = pd.read_parquet(DATA_PQ, columns=["UMAP1", "UMAP2", "study"])
print(f"  {len(df):,} cells, {df['study'].nunique()} studies")

rng = np.random.default_rng(SEED)
bg_idx = rng.choice(len(df), size=min(500_000, len(df)), replace=False)
bg_umap1 = df["UMAP1"].values[bg_idx]
bg_umap2 = df["UMAP2"].values[bg_idx]

# ============================================================================
# PLOT
# ============================================================================

print("Plotting...", flush=True)

n_studies = len(STUDY_ORDER)
ncols = 4
nrows = int(np.ceil(n_studies / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(180 / 25.4, nrows * 42 / 25.4),
                         gridspec_kw={"wspace": 0.02, "hspace": 0.15})
axes = np.array(axes).flatten()

for i, study in enumerate(STUDY_ORDER):
    ax = axes[i]
    ax.scatter(bg_umap1, bg_umap2, c="#E8E8E8", s=0.01, alpha=0.3,
               linewidths=0, rasterized=True)

    mask = df["study"] == study
    study_df = df.loc[mask]
    if len(study_df) > 200_000:
        study_df = study_df.sample(n=200_000, random_state=SEED)

    ax.scatter(study_df["UMAP1"], study_df["UMAP2"], c=STUDY_PALETTE[study],
               s=0.05, alpha=0.6, linewidths=0, rasterized=True)
    ax.set_axis_off()
    ax.text(0.5, -0.01, f"{STUDY_DISPLAY[study]}\nn = {mask.sum():,}",
            transform=ax.transAxes, ha="center", va="top", fontsize=FN, color="#333333")

for j in range(n_studies, len(axes)):
    axes[j].set_visible(False)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight")
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_SVG}")
print(f"Saved: {OUT_PNG}")
