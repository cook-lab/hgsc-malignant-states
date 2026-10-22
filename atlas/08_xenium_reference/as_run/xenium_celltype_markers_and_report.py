#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/xenium_celltype_markers_and_report.py, the script that ran.
# Kept: the 1,000-cells-per-type downsample written as xenium_celltype_downsampled.h5ad (the Xenium SingleR reference; Figure S14).
# Removed: the one-vs-rest marker CSVs, the panel-coverage analysis and its HTML report.
# The untrimmed script is in the lab project folder.
"""
Xenium Cell Type — Markers + Coverage Report
=============================================
1. Downsample hgsc_atlas_xenium.h5ad to 1000 cells per xenium_celltype
2. Run Wilcoxon one-vs-rest DEG markers
3. Cross-reference with xenium panel genes
4. Produce HTML coverage report
"""

import gc, os, glob, warnings
import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")

BASE       = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD_PATH  = os.path.join(BASE, "hgsc_atlas_xenium.h5ad")
OUT_DIR    = os.path.join(BASE, "output/xenium_celltype")
MARKER_DIR = os.path.join(OUT_DIR, "markers")
OUT_H5AD   = os.path.join(OUT_DIR, "xenium_celltype_downsampled.h5ad")

os.makedirs(MARKER_DIR, exist_ok=True)

MAX_CELLS = 1000

LABEL_ORDER = [
    "Secretory epithelium",
    "Ciliated epithelium",
    "Mesothelial",
    "Fibroblast",
    "Smooth muscle",
    "Pericyte",
    "Endothelial",
    "T cell",
    "NK cell",
    "B cell",
    "Plasma cell",
    "Macrophage",
    "Conventional dendritic cell",
    "Plasmacytoid dendritic cell",
    "Neutrophil",
    "Mast cell",
]

# ── load & downsample ─────────────────────────────────────────────────
print(f"\nLoading {H5AD_PATH}...", flush=True)
adata = sc.read_h5ad(H5AD_PATH, backed="r")
print(f"  Shape: {adata.shape[0]:,} × {adata.shape[1]:,}")

obs = adata.obs[["xenium_celltype"]].copy()
obs["xenium_celltype"] = obs["xenium_celltype"].astype(str)

print(f"\nDownsampling to {MAX_CELLS:,} cells per type...", flush=True)
rng = np.random.default_rng(42)
keep_idx = []
for ct in LABEL_ORDER:
    ct_idx = obs.index[obs["xenium_celltype"] == ct].tolist()
    n_take = min(MAX_CELLS, len(ct_idx))
    if len(ct_idx) > MAX_CELLS:
        ct_idx = rng.choice(ct_idx, MAX_CELLS, replace=False).tolist()
    keep_idx.extend(ct_idx)
    print(f"  {ct:<35} {n_take:>6,} cells")

print(f"  Total: {len(keep_idx):,} cells")

print("Loading subset into memory...", flush=True)
adata_sub = adata[keep_idx].to_memory()
adata.file.close()
del adata
gc.collect()

print(f"  In-memory shape: {adata_sub.shape}")

# ── save downsampled h5ad ──────────────────────────────────────────────
print(f"\nSaving downsampled h5ad → {OUT_H5AD}...", flush=True)
adata_sub.write_h5ad(OUT_H5AD)
fsize = os.path.getsize(OUT_H5AD) / (1024**3)
print(f"  File size: {fsize:.2f} GB")

del adata_sub
gc.collect()
