#!/usr/bin/env python3
"""
Figure S1 — Atlas QC metrics by study with filtering cutoffs
============================================================

Purpose
    Three violin panels (one file each), each showing the unfiltered distribution of a
    QC metric per study with the filtering threshold overlaid:
      a. Total UMI counts (log scale) — cutoff < 500
      b. Genes detected — cutoff < 300
      c. Doublet score — cutoff > 0.25
    UMI / genes use the raw pre-filter parquet (2.73M cells). Doublet score
    uses the post-filter parquet (scores computed before doublet removal).

INPUTS
    output_root/_prep_caches/atlas_obs_prefilter.parquet
    output_root/_prep_caches/atlas_obs_postfilter.parquet
    (produced by 00_figS1_extract_atlas_obs.py — run that first)

OUTPUTS
    output_root/figures/supplementary/SF1A_atlas_qc_total_counts_by_study.{svg,png}
    output_root/figures/supplementary/SF1B_atlas_qc_genes_by_study.{svg,png}
    output_root/figures/supplementary/SF1C_atlas_qc_doublet_score_by_study.{svg,png}

MANUSCRIPT PANEL(S)
    Figure S1a-c.

RUNTIME TIER
    fast (reads cached parquets; every cell is drawn).
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

obs_pre = pd.read_parquet(path("output_root", "_prep_caches", "atlas_obs_prefilter.parquet"))
obs_post = pd.read_parquet(path("output_root", "_prep_caches", "atlas_obs_postfilter.parquet"))

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
        "olbrecht": "Olbrecht", "regner": "Regner", "vazquez_garcia": "Vázquez-García",
        "xu": "Xu", "zhang": "Zhang", "zheng": "Zheng",
    }
    return f"{nice.get(name, name)} {year}"


def plot_violins(ax, obs_df, metric, studies, palette, log_scale,
                 thresh_val, thresh_dir, thresh_label):
    """Returns the number of cells drawn."""
    data_by_study = []
    colors = []
    for study in studies:
        # Every cell is used (no subsampling) — final-render rule.
        vals = obs_df.loc[obs_df["study"] == study, metric].dropna().values.astype(float)
        if log_scale:
            # KDE on the log10 scale (a linear-space KDE drawn on a log axis distorts the
            # violins and runs to -inf for low-count cells). Zero-count cells clip to 1.
            vals = np.log10(np.clip(vals, 1, None))
        data_by_study.append(vals)
        colors.append(palette.get(study, "#A0A0A0"))

    if not data_by_study:
        return 0
    t = np.log10(thresh_val) if log_scale else thresh_val

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

    ax.axhline(t, color="#d73027", linewidth=0.8, linestyle="--", zorder=5)
    ax.text(1.005, t, thresh_label, transform=ax.get_yaxis_transform(),
            va="center", ha="left", fontsize=FN, color="#d73027", clip_on=False)

    lo, hi = ax.get_ylim()
    if thresh_dir == "min":
        ax.axhspan(lo, t, color="#d73027", alpha=0.06, zorder=0, linewidth=0)
    else:
        ax.axhspan(t, hi, color="#d73027", alpha=0.06, zorder=0, linewidth=0)
    ax.set_ylim(lo, hi)

    labels = [format_study_label(s) for s in studies]
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=FK)

    if log_scale:
        ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{10 ** x:,.0f}"))
    ax.spines[["top", "right"]].set_visible(False)
    return int(sum(len(v) for v in data_by_study))


def draw_panel(ax, metric, ylabel, log_scale, thresh_val, thresh_dir, thresh_label, obs_df):
    studies = [s for s in STUDY_ORDER if s in obs_df["study"].unique()]
    n = plot_violins(ax, obs_df, metric, studies, STUDY_PALETTE, log_scale,
                     thresh_val, thresh_dir, thresh_label)
    ax.set_ylabel(ylabel, fontsize=FA)
    ax.text(0.02, 0.97, f"n = {n:,} cells", transform=ax.transAxes, va="top", ha="left",
            fontsize=FN, color="#555555")


# ============================================================================
# BUILD FIGURES — one standalone file per panel (a, b, c)
# ============================================================================

print("Plotting SF1 — QC metrics by study...", flush=True)

PANEL_STEMS = ["SF1A_atlas_qc_total_counts_by_study",
               "SF1B_atlas_qc_genes_by_study",
               "SF1C_atlas_qc_doublet_score_by_study"]
for stem, spec in zip(PANEL_STEMS, METRICS):
    fig, ax = plt.subplots(figsize=(180 / 25.4, 60 / 25.4))
    draw_panel(ax, *spec)
    for ext in ("svg", "png"):
        out = path("output_root", "figures", "supplementary", f"{stem}.{ext}")
        fig.savefig(out, format=ext, dpi=450, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)
