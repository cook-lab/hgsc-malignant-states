#!/usr/bin/env python3
"""
Atlas 01 — Step 05: publication UMAP suite

PURPOSE
    Generate individual publication-quality UMAP SVGs for all key metadata
    variables on the post-filter atlas (celltype_pred, study, leiden, anatomic
    site, treatment, metastatic site, doublet score, UMI, genes). Reads obs +
    UMAP only.

INPUTS
    obj("atlas_scanvi")  = hgsc_atlas_scanvi.h5ad  (written by 03b).

OUTPUTS
    output_root/01_preprocess_qc/05_umap_suite/05_umap_*.svg  (9 UMAPs)

MANUSCRIPT PANEL(S)
    Atlas overview UMAP substrate (Fig 1A and SF2/SF3 are rendered by per-panel
    figures/ scripts from the same object).

RUNTIME TIER
    moderate (backed read of obs + UMAP).
"""

import gc
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

warnings.filterwarnings("ignore")
np.random.seed(SEED)

# ============================================================================
# PATHS (resolved via central config)
# ============================================================================

H5AD_PATH = obj("atlas_scanvi")
FIG_DIR   = path("output_root", "01_preprocess_qc", "05_umap_suite")
os.makedirs(FIG_DIR, exist_ok=True)

# ============================================================================
# STYLE — Cook Lab v1.2
# ============================================================================

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "legend.fontsize": 6,
    "figure.dpi":      450,
    "savefig.dpi":     450,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
})

CELLTYPE_PALETTE = {
    "Epithelial":    "#E6A141",
    "Mesothelial":   "#A8A298",
    "Fibroblast":    "#DDD5CA",
    "Smooth_Muscle": "#D14E6C",
    "Pericyte":      "#B87A7A",
    "Endothelial":   "#7D4E4E",
    "T_cell":        "#87CEFA",
    "NK_cell":       "#56AFC4",
    "B_cell":        "#5665B6",
    "Plasma_cell":   "#8A5DAF",
    "Macrophage":    "#8FBC8F",
    "DC":            "#2E8B57",
    "Neutrophil":    "#6B8E23",
    "Mast":          "#8B9B6B",
    "Erythrocyte":   "#CD5C5C",
}

KELLY_22 = [
    "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
    "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379",
    "#604E97", "#F6A600", "#B3446C", "#882D17", "#8DB600",
    "#654522", "#E25822", "#2B3D26", "#CC79A7", "#56B4E9",
    "#009E73", "#D55E00",
]

TREATMENT_PALETTE = {
    "pre-treatment":               "#DDD5CA",
    "post-treatment":              "#56B4E9",
    "post-chemotherapy":           "#009E73",
    "post-chemotherapy_niraparib": "#E69F00",
    "post-chemotherapy_olaparib":  "#0072B2",
    "post-chemotherapy_pembro":    "#D55E00",
    "post-niraparib":              "#CC79A7",
    "NA":                          "#666666",
}

METASTATIC_PALETTE = {
    "primary":    "#7A9EBF",
    "metastasis": "#B07AA1",
    "ascites":    "#8FAC8C",
    "healthy":    "#C2956B",
}


def auto_palette(categories, base):
    return {cat: base[i % len(base)] for i, cat in enumerate(sorted(categories))}


# ============================================================================
# LOAD DATA (backed → lightweight extraction → close)
# ============================================================================

print("=" * 60)
print("Step 5 — UMAP Visualization Suite")
print("=" * 60)
print(f"\nLoading (backed): {H5AD_PATH}", flush=True)

adata = sc.read_h5ad(H5AD_PATH, backed="r")
print(f"  Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

OBS_COLS = [
    "celltype_pred", "study", "leiden",
    "anatomic_site", "treatment_status", "metastatic_site",
    "doublet_score_scrublet", "total_counts", "n_genes_by_counts",
]
obs = adata.obs[OBS_COLS].copy()
umap = np.array(adata.obsm["X_umap"])
obs["UMAP1"] = umap[:, 0]
obs["UMAP2"] = umap[:, 1]

adata.file.close()
del adata, umap
gc.collect()

for col in ["doublet_score_scrublet", "total_counts", "n_genes_by_counts"]:
    obs[col] = pd.to_numeric(obs[col], errors="coerce")


# ============================================================================
# PLOTTING HELPERS
# ============================================================================

def umap_categorical(obs_df, col, palette, title, fname,
                     figsize=(5.5, 4), pt_size=0.02, legend_ncol=1):
    fig, ax = plt.subplots(figsize=figsize)
    data = obs_df.dropna(subset=[col]).copy()
    data[col] = data[col].astype(str)
    data = data.sample(frac=1, random_state=SEED)
    categories = sorted(data[col].unique())
    for cat in categories:
        mask = data[col] == cat
        ax.scatter(data.loc[mask, "UMAP1"], data.loc[mask, "UMAP2"],
                   c=palette.get(cat, "#A0A0A0"), s=pt_size, alpha=0.6,
                   linewidths=0, rasterized=True, label=cat)
    ax.set_title(title, fontsize=9, pad=6)
    ax.set_axis_off()
    handles = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=palette.get(cat, "#A0A0A0"),
               markersize=4, linewidth=0, label=cat)
        for cat in categories
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              ncol=legend_ncol, fontsize=5.5, frameon=False,
              markerscale=1.2, handletextpad=0.3, borderpad=0.2)
    fig.savefig(os.path.join(FIG_DIR, fname), format="svg", dpi=450, bbox_inches="tight")
    plt.close(fig)


def umap_continuous(obs_df, col, cmap, title, fname,
                    figsize=(5, 4), pt_size=0.02, log_scale=False):
    fig, ax = plt.subplots(figsize=figsize)
    data = obs_df.dropna(subset=[col]).copy()
    data[col] = pd.to_numeric(data[col], errors="coerce").pipe(
        lambda s: np.log1p(s) if log_scale else s
    )
    data = data.dropna(subset=[col]).sort_values(col)
    sc_plot = ax.scatter(data["UMAP1"], data["UMAP2"], c=data[col], cmap=cmap,
                         s=pt_size, alpha=0.6, linewidths=0, rasterized=True)
    cbar = fig.colorbar(sc_plot, ax=ax, shrink=0.5, pad=0.02, aspect=20)
    cbar.ax.tick_params(labelsize=6)
    label = (f"log1p({col.replace('_', ' ').title()})"
             if log_scale else col.replace("_", " ").title())
    cbar.set_label(label, fontsize=6)
    ax.set_title(title, fontsize=9, pad=6)
    ax.set_axis_off()
    fig.savefig(os.path.join(FIG_DIR, fname), format="svg", dpi=450, bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# RENDER
# ============================================================================

umap_categorical(obs, "celltype_pred", CELLTYPE_PALETTE,
                 title="Cell type (CellAssign v3)",
                 fname="05_umap_celltype_pred.svg", figsize=(5.5, 4))

study_palette = auto_palette(obs["study"].dropna().unique(), KELLY_22)
umap_categorical(obs, "study", study_palette, title="Study",
                 fname="05_umap_study.svg", figsize=(6, 4))

leiden_cats = sorted(obs["leiden"].dropna().astype(str).unique(), key=lambda x: int(x))
leiden_palette = {cat: KELLY_22[i % len(KELLY_22)] for i, cat in enumerate(leiden_cats)}
umap_categorical(obs, "leiden", leiden_palette, title="Leiden clusters",
                 fname="05_umap_leiden.svg", figsize=(5.5, 4))

anatomic_palette = auto_palette(obs["anatomic_site"].dropna().unique(), KELLY_22)
umap_categorical(obs, "anatomic_site", anatomic_palette, title="Anatomic site",
                 fname="05_umap_anatomic_site.svg", figsize=(6, 4))

umap_categorical(obs, "treatment_status", TREATMENT_PALETTE, title="Treatment status",
                 fname="05_umap_treatment_status.svg", figsize=(6, 4))

umap_categorical(obs, "metastatic_site", METASTATIC_PALETTE, title="Metastatic site",
                 fname="05_umap_metastatic_site.svg", figsize=(5.5, 4))

umap_continuous(obs, "doublet_score_scrublet", cmap="viridis",
                title="Scrublet doublet score", fname="05_umap_doublet_score.svg")

umap_continuous(obs, "total_counts", cmap="magma", title="Total UMI counts (log1p)",
                fname="05_umap_total_counts.svg", log_scale=True)

umap_continuous(obs, "n_genes_by_counts", cmap="magma", title="Genes detected (log1p)",
                fname="05_umap_n_genes.svg", log_scale=True)

print("Done — 9 UMAP SVGs saved to output_root/01_preprocess_qc/05_umap_suite/")
