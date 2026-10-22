#!/usr/bin/env python3
"""
Figure S3d-f — Atlas UMAPs by TP53, HRD, and BRCA status
========================================================

Purpose
    Three categorical UMAP panels of the full atlas coloured by genomic status,
    matching the Figure S3a-c clinical-UMAP style.

INPUTS
    figures._prep.atlas_final_tables.load_obs()  (final atlas obj("atlas_final"): obs
        TP53_status, HRD_status, BRCA_status + obsm X_umap -> UMAP1/UMAP2; all 1,980,703
        cells). Values are plotted as held in the final atlas (no recode; "False"/"nan"
        are drawn as NA).

OUTPUTS
    output_root/figures/supplementary/SF3_tp53_status.{svg,png}
    output_root/figures/supplementary/SF3_hrd_status.{svg,png}
    output_root/figures/supplementary/SF3_brca_status.{svg,png}

MANUSCRIPT PANEL(S)
    Figure S3d (TP53), S3e (HRD), S3f (BRCA).

RUNTIME TIER
    moderate (loads the obs cache; plots all 1,980,703 points, rasterized).
"""

import os
import sys

import numpy as np
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
# COLOR PALETTES
# ============================================================================

TP53_PALETTE = {"mutated": "#D14E6C", "wildtype": "#87CEFA", "NA": "#DDD5CA"}
TP53_LABELS = {"mutated": "TP53 Mutated", "wildtype": "TP53 Wildtype", "NA": "Unknown"}
HRD_PALETTE = {"HRD": "#D14E6C", "HRP": "#87CEFA", "NA": "#DDD5CA"}
HRD_LABELS = {"HRD": "HRD", "HRP": "HR Proficient", "NA": "Unknown"}
BRCA_PALETTE = {"mutated": "#D14E6C", "wildtype": "#87CEFA", "NA": "#DDD5CA"}
BRCA_LABELS = {"mutated": "BRCA Mutated", "wildtype": "BRCA Wildtype", "NA": "Unknown"}

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading final-atlas obs (load_obs)...", flush=True)
obs = load_obs(columns=["TP53_status", "HRD_status", "BRCA_status", "UMAP1", "UMAP2"])
print(f"  {len(obs):,} cells")

for col in ["TP53_status", "HRD_status", "BRCA_status"]:
    print(f"  {col} categories in final atlas: {obs[col].astype(str).value_counts(dropna=False).to_dict()}")
    obs[col] = obs[col].astype(str)
    obs.loc[obs[col] == "False", col] = "NA"
    obs.loc[obs[col] == "nan", col] = "NA"

# Every cell is plotted (no subsampling); the scatter layer is rasterized.


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

    legend_order = [k for k in order if k != "NA"]
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=palette[k],
               markersize=3, linewidth=0, label=labels.get(k, k))
        for k in legend_order
    ]
    if "NA" in order:
        handles.append(
            Line2D([0], [0], marker="o", color="w", markerfacecolor=palette["NA"],
                   markersize=2.5, linewidth=0, label=labels.get("NA", "NA"))
        )
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


print("\nPanel 1/3: TP53 status", flush=True)
plot_umap_categorical(obs, "TP53_status", TP53_PALETTE, TP53_LABELS, "SF3_tp53_status")

print("Panel 2/3: HRD status", flush=True)
plot_umap_categorical(obs, "HRD_status", HRD_PALETTE, HRD_LABELS, "SF3_hrd_status")

print("Panel 3/3: BRCA status", flush=True)
plot_umap_categorical(obs, "BRCA_status", BRCA_PALETTE, BRCA_LABELS, "SF3_brca_status")

print("\nDone — 3 genomic UMAP panels saved.")
