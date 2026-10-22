#!/usr/bin/env python3
"""
Figure S3a-c — Atlas UMAPs by treatment status, anatomic site, metastatic site
==============================================================================

Purpose
    Three categorical UMAP panels of the full atlas coloured by clinical
    metadata. As in Figure S3, panel a is coloured by treatment status, b by
    anatomic site, c by metastatic site. (Panels d-f — TP53 / HRD / BRCA — are
    produced by the companion script 06_figS3def_genomic_umaps.py.) No titles /
    no panel letters.

INPUTS
    figures._prep.atlas_final_tables.load_obs()  (final atlas obj("atlas_final"): obs
        treatment_status, anatomic_site, metastatic_site (final, corrected labels) + obsm
        X_umap -> UMAP1/UMAP2; all 1,980,703 cells).

OUTPUTS
    output_root/figures/supplementary/SF3_treatment_status.{svg,png}
    output_root/figures/supplementary/SF3_anatomic_site.{svg,png}
    output_root/figures/supplementary/SF3_metastatic_site.{svg,png}

MANUSCRIPT PANEL(S)
    Figure S3a (treatment status), S3b (anatomic site), S3c (metastatic site).

RUNTIME TIER
    moderate (loads the obs cache; plots all 1,980,703 points, rasterized).
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402
from figures._prep.atlas_final_tables import load_obs  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

def out_path(stem, ext):
    return path("output_root", "figures", "supplementary", f"{stem}.{ext}")

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
# LOAD DATA
# ============================================================================

print("Loading final-atlas obs (load_obs)...", flush=True)
obs = load_obs(columns=["treatment_status", "anatomic_site", "metastatic_site", "UMAP1", "UMAP2"])
print(f"  {len(obs):,} cells")

# Every cell is plotted (no subsampling); the scatter layer is rasterized.

# ============================================================================
# COLOR PALETTES
# ============================================================================

ANATOMIC_PALETTE = {
    "adnexa": "#E6A141", "omentum": "#9A7D55", "ascites": "#87CEFA",
    "peritoneum": "#8FBC8F", "bowel": "#D14E6C", "upper_quadrant": "#5665B6",
    "pelvic_cavity": "#56AFC4", "colon": "#B87A7A", "bladder": "#8A5DAF",
    "liver": "#7D4E4E", "diaphragm": "#DDD5CA", "lymph_node": "#6B8E23",
}
METASTATIC_PALETTE = {
    "primary": "#E6A141", "metastasis": "#5665B6", "ascites": "#87CEFA", "healthy": "#DDD5CA",
}
# Treatment-status palette / labels (Figure S3a legend order):
#   Pre-treatment, Post-chemo, Post-niraparib, NA,
#   Post-chemo + niraparib, Post-chemo + pembro, Post-chemo + olaparib.
# Keys match the raw treatment_status categories in atlas_final.obs.
TREATMENT_PALETTE = {
    "pre-treatment":               "#56AFC4",
    "post-chemotherapy":           "#D14E6C",
    "post-niraparib":              "#E6A141",
    "NA":                          "#DDD5CA",
    "post-chemotherapy_niraparib": "#8A5DAF",
    "post-chemotherapy_pembro":    "#2E8B57",
    "post-chemotherapy_olaparib":  "#5665B6",
    # Nath 2021 P19 (chemo, niraparib and pembrolizumab); present since the Nath 2021 sample mapping fix (2026-10-20)
    "post-chemotherapy_niraparib_pembro": "#B87A7A",
}
TREATMENT_LABELS = {
    "pre-treatment":               "Pre-treatment",
    "post-chemotherapy":           "Post-chemo",
    "post-niraparib":              "Post-niraparib",
    "NA":                          "NA",
    "post-chemotherapy_niraparib": "Post-chemo + niraparib",
    "post-chemotherapy_pembro":    "Post-chemo + pembro",
    "post-chemotherapy_olaparib":  "Post-chemo + olaparib",
    "post-chemotherapy_niraparib_pembro": "Post-chemo + niraparib + pembro",
}
ANATOMIC_LABELS = {
    "adnexa": "Adnexa", "omentum": "Omentum", "ascites": "Ascites",
    "peritoneum": "Peritoneum", "bowel": "Bowel", "upper_quadrant": "Upper quadrant",
    "pelvic_cavity": "Pelvic cavity", "colon": "Colon", "bladder": "Bladder",
    "liver": "Liver", "diaphragm": "Diaphragm", "lymph_node": "Lymph node",
}
METASTATIC_LABELS = {
    "primary": "Primary", "metastasis": "Metastasis", "ascites": "Ascites", "healthy": "Healthy",
}


def plot_umap_categorical(df, col, palette, labels, out_stem):
    """Single UMAP panel with legend, no title, no panel label."""
    order = [k for k in palette if k in df[col].values]
    fig, ax = plt.subplots(figsize=(88 / 25.4, 75 / 25.4))
    plot_df = df.sample(frac=1, random_state=SEED)

    # Draw "NA"/unknown first as a background layer, then every remaining cell in ONE
    # shuffled scatter so no category is systematically painted over the others.
    is_na = (plot_df[col] == "NA").values
    for sel in (is_na, ~is_na):
        sub = plot_df.loc[sel]
        if len(sub):
            cols = sub[col].astype(str).map(palette).fillna("#DDD5CA").values
            ax.scatter(sub["UMAP1"].values, sub["UMAP2"].values, c=cols, s=0.02,
                       alpha=0.5, linewidths=0, rasterized=True)

    ax.set_axis_off()
    ax.set_aspect("equal")
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[k],
               markersize=3, linewidth=0, label=labels.get(k, k))
        for k in order
    ]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=False, fontsize=FN, handletextpad=0.3, borderpad=0.2,
              labelspacing=0.25, markerscale=1.2)

    out_svg = out_path(out_stem, "svg")
    out_png = out_path(out_stem, "png")
    fig.savefig(out_svg, format="svg", dpi=450, bbox_inches="tight")
    fig.savefig(out_png, format="png", dpi=450, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_svg}")
    print(f"  Saved: {out_png}")


print("\nPanel 1/3 (SF3A): treatment status", flush=True)
plot_umap_categorical(obs, "treatment_status", TREATMENT_PALETTE, TREATMENT_LABELS,
                      "SF3_treatment_status")

print("Panel 2/3 (SF3B): anatomic site", flush=True)
plot_umap_categorical(obs, "anatomic_site", ANATOMIC_PALETTE, ANATOMIC_LABELS,
                      "SF3_anatomic_site")

print("Panel 3/3 (SF3C): metastatic site", flush=True)
plot_umap_categorical(obs, "metastatic_site", METASTATIC_PALETTE, METASTATIC_LABELS,
                      "SF3_metastatic_site")

print("\nDone — 3 panels saved.")
