#!/usr/bin/env python3
"""
SF1 — Atlas QC metrics by study with filtering cutoffs
======================================================

Purpose
    Three stacked violin panels, each showing the unfiltered distribution of a
    QC metric per study with the filtering threshold overlaid:
      1. Total UMI counts (log scale) — cutoff < 500
      2. Genes detected — cutoff < 300
      3. Doublet score — cutoff > 0.25
    UMI / genes use the raw pre-filter parquet (2.73M cells). Doublet score
    uses the post-filter parquet (scores computed before doublet removal).

INPUTS
    output_root/figures/data/atlas_obs_prefilter.parquet
    output_root/figures/data/atlas_obs_postfilter.parquet
    (produced by 00_extract_atlas_obs.py — run that first)

OUTPUTS
    output_root/figures/supplementary/SF1_atlas_qc_metrics_by_study.{svg,png}

MANUSCRIPT PANEL(S)
    SF1A-C.

RUNTIME TIER
    fast (reads cached parquets; subsamples to 50k cells per study for violins).
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

obs_pre = pd.read_parquet(path("output_root", "figures", "data", "atlas_obs_prefilter.parquet"))
obs_post = pd.read_parquet(path("output_root", "figures", "data", "atlas_obs_postfilter.parquet"))

OUT_SVG = path("output_root", "figures", "supplementary", "SF1_atlas_qc_metrics_by_study.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF1_atlas_qc_metrics_by_study.png")

# ============================================================================
# STYLE
# ============================================================================

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

# ============================================================================
# STUDY PALETTE & ORDER
# ============================================================================

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
STUDY_ORDER = list(STUDY_PALETTE.keys())

# METRICS — (column, ylabel, log_scale, threshold, direction, label, source_df)
METRICS = [
    ("total_counts",           "Total UMI counts", True,  500,  "min", "< 500 UMI",  obs_pre),
    ("n_genes_by_counts",      "Genes detected",   False, 300,  "min", "< 300 genes", obs_pre),
    ("doublet_score_scrublet", "Doublet score",    False, 0.25, "max", "> 0.25",      obs_post),
]


def format_study_label(s):
    parts = s.split("_")
    year = parts[-1]
    name = "_".join(parts[:-1])
    nice = {
        "denisenko": "Denisenko", "geistlinger": "Geistlinger", "hornburg": "Hornburg",
        "loret": "Loret", "luo": "Luo", "nath": "Nath", "olalekan": "Olalekan",
        "olbrecht": "Olbrecht", "regner": "Regner", "vazquez_garcia": "Vazquez-Garcia",
        "xu": "Xu", "zhang": "Zhang", "zheng": "Zheng",
    }
    return f"{nice.get(name, name)} {year}"


def plot_violins(ax, obs_df, metric, studies, palette, log_scale,
                 thresh_val, thresh_dir, thresh_label):
    data_by_study = []
    colors = []
    for study in studies:
        vals = obs_df.loc[obs_df["study"] == study, metric].dropna().values
        if len(vals) > 50000:
            rng = np.random.default_rng(SEED)
            vals = rng.choice(vals, size=50000, replace=False)
        data_by_study.append(vals)
        colors.append(palette.get(study, "#A0A0A0"))

    if not data_by_study:
        return

    positions = list(range(len(studies)))
    parts = ax.violinplot(
        data_by_study, positions=positions, showmeans=False,
        showmedians=True, showextrema=False, widths=0.75,
    )
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(colors[i])
        body.set_edgecolor("none")
        body.set_alpha(0.85)
    parts["cmedians"].set_edgecolor("#333333")
    parts["cmedians"].set_linewidth(0.6)

    ax.axhline(thresh_val, color="#d73027", linewidth=0.8, linestyle="--", zorder=5)
    ax.text(len(studies) - 0.5, thresh_val, f"  {thresh_label}",
            va="bottom", ha="right", fontsize=FN, color="#d73027")

    if thresh_dir == "min":
        ax.axhspan(ax.get_ylim()[0] if not log_scale else 0.1,
                   thresh_val, color="#d73027", alpha=0.04, zorder=0)
    else:
        ax.axhspan(thresh_val, ax.get_ylim()[1], color="#d73027", alpha=0.04, zorder=0)

    labels = [format_study_label(s) for s in studies]
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=FK)

    if log_scale:
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(
            lambda x, _: f"{int(x):,}" if x >= 1 else f"{x:.2f}"
        ))
    ax.spines[["top", "right"]].set_visible(False)


# ============================================================================
# BUILD FIGURE — 3 rows × 1 column
# ============================================================================

print("Plotting SF1 — QC metrics by study...", flush=True)

fig, axes = plt.subplots(3, 1, figsize=(180 / 25.4, 180 / 25.4),
                         gridspec_kw={"hspace": 0.55})

for row, (metric, ylabel, log_scale, thresh_val, thresh_dir, thresh_label, obs_df) in enumerate(METRICS):
    ax = axes[row]
    studies = [s for s in STUDY_ORDER if s in obs_df["study"].unique()]
    plot_violins(ax, obs_df, metric, studies, STUDY_PALETTE, log_scale,
                 thresh_val, thresh_dir, thresh_label)
    ax.set_ylabel(ylabel, fontsize=FA)
    n = obs_df[metric].notna().sum()
    ax.text(0.02, 0.97, f"n = {n:,}", transform=ax.transAxes, va="top", ha="left",
            fontsize=FN, color="#555555")

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight")
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {OUT_SVG}")
print(f"Saved: {OUT_PNG}")
