#!/usr/bin/env python3
"""
Atlas 01 — Step 08: re-filter doublets (threshold 0.25) + recompute UMAP

PURPOSE
    Apply the stricter post-integration Scrublet filter (score < 0.25, was 0.275),
    recompute neighbours/UMAP on the scANVI latent space (obsm["X_scanvi"]), and
    write the canonical downstream entry-point object. This consumes
    integrated_scanvi.h5ad produced by 07_finalize.py (--method scanvi), which emits
    the matching lowercase obsm key "X_scanvi" read here.

INPUTS
    integrated_scanvi.h5ad  (scANVI integration output; written by 07_finalize.py,
        or the deposited 20260213 run). Resolved under
        DATA_ROOT/2026_final_atlas/pre-processing and integration/...

OUTPUTS
    obj("atlas_scanvi")  = hgsc_atlas_scanvi.h5ad   (post-filter + UMAP; config entry-point)
    output_root/01_preprocess_qc/08_refilter/*.svg  (QC UMAP checks)

MANUSCRIPT PANEL(S)
    Upstream backend; underpins all atlas panels. No panel rendered directly.

RUNTIME TIER
    heavy (loads ~2.3M-cell h5ad; neighbours + UMAP).
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

# Prefer the freshly finalized integration output (07_finalize writes here); fall
# back to the deposited integrated object if a clean re-run has not produced it yet
# (audit H18 — the 07->08 chain is only continuous once 07 has run locally).
_RECOMPUTED = path("output_root", "01_preprocess_qc", "integration", "anndata", "integrated_scanvi.h5ad")
_DEPOSITED = path(
    "data_root", "2026_final_atlas", "pre-processing and integration",
    "20260213_integration", "output", "anndata", "integrated_scanvi.h5ad",
)
INPUT_H5AD = _RECOMPUTED if os.path.exists(_RECOMPUTED) else _DEPOSITED
OUTPUT_H5AD = obj("atlas_scanvi")
FIG_DIR     = path("output_root", "01_preprocess_qc", "08_refilter")
os.makedirs(FIG_DIR, exist_ok=True)

SCORE_THRESHOLD = 0.25  # stricter than previous 0.275

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

# ============================================================================
# LOAD
# ============================================================================

print("=" * 60)
print("Step 08 — Re-filter doublets (threshold 0.25) + UMAP")
print("=" * 60)

print(f"\n1. Loading: {INPUT_H5AD}", flush=True)
adata = sc.read_h5ad(INPUT_H5AD)
print(f"   Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# ============================================================================
# FILTER 1 — Remove doublet_scrublet == True
# ============================================================================

print("\n2. Removing doublet_scrublet == True...", flush=True)
n_before = adata.n_obs

doublet_col = adata.obs["doublet_scrublet"]
if doublet_col.dtype == object:
    mask_remove = doublet_col.astype(str).str.lower() == "true"
else:
    mask_remove = doublet_col == True
n_doublets = int(mask_remove.sum())
adata._inplace_subset_obs(~mask_remove)

print(f"   Before:  {n_before:,}")
print(f"   Removed: {n_doublets:,} (doublet_scrublet == True)")
print(f"   After:   {adata.n_obs:,}")

# ============================================================================
# FILTER 2 — Score threshold 0.25
# ============================================================================

print(f"\n3. Removing doublet_score_scrublet >= {SCORE_THRESHOLD}...", flush=True)
n_before2 = adata.n_obs
mask_keep2 = adata.obs["doublet_score_scrublet"].fillna(0) < SCORE_THRESHOLD
n_removed2 = int((~mask_keep2).sum())
adata._inplace_subset_obs(mask_keep2)

print(f"   Before:  {n_before2:,}")
print(f"   Removed: {n_removed2:,} (score >= {SCORE_THRESHOLD})")
print(f"   After:   {adata.n_obs:,}")

total_removed = n_doublets + n_removed2
print(f"\n   Total removed: {total_removed:,} ({total_removed / n_before * 100:.2f}%)")
print(f"   Final shape:   {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# ============================================================================
# NEIGHBORS + UMAP
# ============================================================================

print("\n4. Computing neighbors (n=10, X_scanvi)...", flush=True)
sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=10)

print("\n5. Computing UMAP (min_dist=0.2)...", flush=True)
sc.tl.umap(adata, min_dist=0.2, random_state=SEED)

# ============================================================================
# SAVE
# ============================================================================

print(f"\n6. Saving: {OUTPUT_H5AD}", flush=True)
adata.write_h5ad(OUTPUT_H5AD)
fsize = os.path.getsize(OUTPUT_H5AD) / (1024**3)
print(f"   File size: {fsize:.1f} GB")

# ============================================================================
# QUICK UMAP PLOTS FOR APPROVAL
# ============================================================================

print("\n7. Plotting UMAPs for approval...", flush=True)

umap = np.array(adata.obsm["X_umap"])
celltype = adata.obs["celltype_pred"].astype(str).values

rng = np.random.default_rng(SEED)
idx = rng.permutation(len(umap))

# --- UMAP by celltype_pred ---
cats = sorted(set(celltype))

fig, ax = plt.subplots(figsize=(5.5, 4))
for cat in cats:
    mask = np.array([celltype[i] == cat for i in idx])
    ax.scatter(
        umap[idx[mask], 0], umap[idx[mask], 1],
        c=CELLTYPE_PALETTE.get(cat, "#A0A0A0"),
        s=0.5, alpha=0.6, linewidths=0, rasterized=True,
    )
handles = [
    Line2D([0], [0], marker="o", color="w",
           markerfacecolor=CELLTYPE_PALETTE.get(c, "#A0A0A0"),
           markersize=4, linewidth=0, label=c)
    for c in cats
]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
          ncol=1, fontsize=5.5, frameon=False, markerscale=1.2,
          handletextpad=0.3, borderpad=0.2)
ax.set_title(f"celltype_pred — {adata.n_obs:,} cells (doublet < {SCORE_THRESHOLD})",
             fontsize=9, pad=6)
ax.set_axis_off()
fig.savefig(os.path.join(FIG_DIR, "08_umap_celltype_pred.svg"),
            format="svg", dpi=450, bbox_inches="tight")
plt.close(fig)

# --- UMAP by doublet score ---
scores = adata.obs["doublet_score_scrublet"].values.astype(float)
sort_idx = np.argsort(scores)

fig, ax = plt.subplots(figsize=(5, 4))
sc_plot = ax.scatter(
    umap[sort_idx, 0], umap[sort_idx, 1],
    c=scores[sort_idx], cmap="viridis",
    s=0.5, alpha=0.6, linewidths=0, rasterized=True,
)
cbar = fig.colorbar(sc_plot, ax=ax, shrink=0.5, pad=0.02, aspect=20)
cbar.ax.tick_params(labelsize=6)
cbar.set_label("Doublet Score", fontsize=6)
ax.set_title(f"Doublet score (threshold < {SCORE_THRESHOLD})", fontsize=9, pad=6)
ax.set_axis_off()
fig.savefig(os.path.join(FIG_DIR, "08_umap_doublet_score.svg"),
            format="svg", dpi=450, bbox_inches="tight")
plt.close(fig)

# --- UMAP by study ---
study = adata.obs["study"].astype(str).values
study_cats = sorted(set(study))
KELLY_22 = [
    "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
    "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379",
    "#604E97", "#F6A600", "#B3446C", "#882D17", "#8DB600",
    "#654522", "#E25822", "#2B3D26", "#CC79A7", "#56B4E9",
    "#009E73", "#D55E00",
]
study_pal = {s: KELLY_22[i % len(KELLY_22)] for i, s in enumerate(study_cats)}

fig, ax = plt.subplots(figsize=(6, 4))
for cat in study_cats:
    mask = np.array([study[i] == cat for i in idx])
    ax.scatter(
        umap[idx[mask], 0], umap[idx[mask], 1],
        c=study_pal[cat],
        s=0.5, alpha=0.6, linewidths=0, rasterized=True,
    )
handles = [
    Line2D([0], [0], marker="o", color="w",
           markerfacecolor=study_pal[c], markersize=4,
           linewidth=0, label=c)
    for c in study_cats
]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
          ncol=1, fontsize=5, frameon=False, markerscale=1.2,
          handletextpad=0.3, borderpad=0.2)
ax.set_title(f"Study — {len(study_cats)} studies", fontsize=9, pad=6)
ax.set_axis_off()
fig.savefig(os.path.join(FIG_DIR, "08_umap_study.svg"),
            format="svg", dpi=450, bbox_inches="tight")
plt.close(fig)

del adata
gc.collect()

print(f"""
{'='*60}
DONE — Step 08 Re-filter complete
{'='*60}
  Doublet threshold: < {SCORE_THRESHOLD} (was 0.275)
  Total removed:     {total_removed:,} from original {n_before:,}
  Final cells:       {n_before - total_removed:,}
  Output: {OUTPUT_H5AD}
""")
