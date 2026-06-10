#!/usr/bin/env python3
"""
SF3A-C — Atlas UMAPs by treatment status, anatomic site, metastatic site
========================================================================

Purpose
    Three categorical UMAP panels of the full atlas coloured by clinical
    metadata. Per the published Supplemental Figure 3, panel A is coloured by
    treatment status, B by anatomic site, C by metastatic site. (Panels D-F of
    the published figure — TP53 / HRD / BRCA — are produced by the companion
    script 06_SF3D_atlas_genomic_umaps.py.) No titles / no panel letters.

INPUTS
    obj("atlas_final")  (hgsc_atlas_final.h5ad; obsm['X_umap'];
        obs: treatment_status, anatomic_site, metastatic_site)

OUTPUTS
    output_root/figures/supplementary/SF3_treatment_status.{svg,png}
    output_root/figures/supplementary/SF3_anatomic_site.{svg,png}
    output_root/figures/supplementary/SF3_metastatic_site.{svg,png}

MANUSCRIPT PANEL(S)
    SF3A (treatment status), SF3B (anatomic site), SF3C (metastatic site).

RUNTIME TIER
    moderate (loads atlas obs + UMAP, subsamples to 800k points).
"""

import os
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
# LOAD DATA
# ============================================================================

print("Loading atlas h5ad (backed)...", flush=True)
adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
print(f"  {adata.n_obs:,} cells")

obs = adata.obs[["treatment_status", "anatomic_site", "metastatic_site"]].copy()
umap = adata.obsm["X_umap"]
obs["UMAP1"] = umap[:, 0]
obs["UMAP2"] = umap[:, 1]
del adata

MAX_CELLS = 800_000
rng = np.random.default_rng(SEED)
if len(obs) > MAX_CELLS:
    idx = rng.choice(len(obs), size=MAX_CELLS, replace=False)
    obs = obs.iloc[idx].copy()
    print(f"  Subsampled to {len(obs):,} cells")

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
# Treatment-status palette / labels (published SF3A legend order):
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
}
TREATMENT_LABELS = {
    "pre-treatment":               "Pre-treatment",
    "post-chemotherapy":           "Post-chemo",
    "post-niraparib":              "Post-niraparib",
    "NA":                          "NA",
    "post-chemotherapy_niraparib": "Post-chemo + niraparib",
    "post-chemotherapy_pembro":    "Post-chemo + pembro",
    "post-chemotherapy_olaparib":  "Post-chemo + olaparib",
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

    for cat in reversed(order):
        mask = plot_df[col] == cat
        ax.scatter(plot_df.loc[mask, "UMAP1"], plot_df.loc[mask, "UMAP2"],
                   c=palette[cat], s=0.02, alpha=0.5, linewidths=0,
                   rasterized=True, label=labels.get(cat, cat))

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
