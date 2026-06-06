#!/usr/bin/env python3
"""
SF3D-F — Atlas UMAPs by TP53, HRD, and BRCA status
==================================================

Purpose
    Three categorical UMAP panels of the full atlas coloured by genomic status,
    matching the SF3A-C clinical-UMAP style.

INPUTS
    obj("atlas_final")  (hgsc_atlas_final.h5ad; obsm['X_umap'];
        obs: TP53_status, HRD_status, BRCA_status)

OUTPUTS
    output_root/figures/supplementary/SF3_tp53_status.{svg,png}
    output_root/figures/supplementary/SF3_hrd_status.{svg,png}
    output_root/figures/supplementary/SF3_brca_status.{svg,png}

MANUSCRIPT PANEL(S)
    SF3D-F.

RUNTIME TIER
    moderate (loads atlas obs + UMAP, subsamples to 800k points).
"""

import os
import gc
import sys

import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

ATLAS_H5AD = obj("atlas_final")


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

print("Loading atlas h5ad (backed)...", flush=True)
adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
print(f"  {adata.n_obs:,} cells")

obs = adata.obs[["TP53_status", "HRD_status", "BRCA_status"]].copy()
umap = adata.obsm["X_umap"]
obs["UMAP1"] = umap[:, 0]
obs["UMAP2"] = umap[:, 1]
adata.file.close()
del adata
gc.collect()

for col in ["TP53_status", "HRD_status", "BRCA_status"]:
    obs[col] = obs[col].fillna("NA")
    obs.loc[obs[col] == "False", col] = "NA"
    obs.loc[obs[col] == "nan", col] = "NA"

MAX_CELLS = 800_000
rng = np.random.default_rng(SEED)
if len(obs) > MAX_CELLS:
    idx = rng.choice(len(obs), size=MAX_CELLS, replace=False)
    obs = obs.iloc[idx].copy()
    print(f"  Subsampled to {len(obs):,} cells")


def plot_umap_categorical(df, col, palette, labels, out_stem):
    """Single UMAP panel with legend, no title, no panel label."""
    order = [k for k in palette if k in df[col].values]
    fig, ax = plt.subplots(figsize=(88 / 25.4, 75 / 25.4))
    plot_df = df.sample(frac=1, random_state=SEED)

    for cat in reversed(order):
        mask = plot_df[col] == cat
        ax.scatter(plot_df.loc[mask, "UMAP1"], plot_df.loc[mask, "UMAP2"],
                   c=palette[cat], s=0.02, alpha=0.5, linewidths=0,
                   rasterized=True, label=labels.get(cat, cat))

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
