#!/usr/bin/env python3
# Trimmed copy (2026-10-22) of 202601_final_atlas/scripts/xenium_coverage_report_level1.py, the script that ran
# (2026-02-23; the marker files are from 15:25 that day). Kept: the one-vs-rest Wilcoxon markers per level-1 type
# (at most 5,000 cells per type, seed 42; adjusted p < 0.05 and log fold change > 0), written to
# output/06_celltype_level1/markers/level1_<type>.csv, the source of Table S2. Removed: the Xenium-panel coverage
# analysis and HTML report that followed, and the panel-gene file they read. The untrimmed script, which the docstring
# below describes, is in the lab project folder.

"""
Xenium Panel Coverage Report — celltype_level1
Computes Wilcoxon DEG markers per level-1 cell type (one-vs-rest),
then checks coverage against the Xenium gene panel and produces
an HTML report flagging any cell types at risk for SingleR annotation.

Downsamples to max 5 000 cells per type for tractable DEG computation.
"""

import gc, os, glob, warnings
import pandas as pd
import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")

# ── paths ──────────────────────────────────────────────────────────────
BASE       = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD_PATH  = os.path.join(BASE, "hgsc_atlas_celltype_level1.h5ad")
OUT_DIR    = os.path.join(BASE, "output/06_celltype_level1")
MARKER_DIR = os.path.join(OUT_DIR, "markers")

os.makedirs(MARKER_DIR, exist_ok=True)

MAX_CELLS_PER_TYPE = 5000

# ── load atlas ─────────────────────────────────────────────────────────
print(f"\nLoading {H5AD_PATH}...", flush=True)
adata = sc.read_h5ad(H5AD_PATH, backed="r")
print(f"  Shape: {adata.shape[0]:,} × {adata.shape[1]:,}")

# Get obs with level1 labels
obs = adata.obs[["celltype_level1"]].copy()
obs["celltype_level1"] = obs["celltype_level1"].astype(str)
celltypes = sorted(obs["celltype_level1"].dropna().unique())
celltypes = [c for c in celltypes if c != "nan"]
print(f"  Level-1 cell types: {celltypes}")

# ── downsample ─────────────────────────────────────────────────────────
print(f"\nDownsampling to max {MAX_CELLS_PER_TYPE:,} cells per type...", flush=True)
keep_idx = []
for ct in celltypes:
    ct_idx = obs.index[obs["celltype_level1"] == ct].tolist()
    if len(ct_idx) > MAX_CELLS_PER_TYPE:
        rng = np.random.default_rng(42)
        ct_idx = rng.choice(ct_idx, MAX_CELLS_PER_TYPE, replace=False).tolist()
    keep_idx.extend(ct_idx)
    print(f"  {ct}: {len(ct_idx):,} cells")

print(f"  Total: {len(keep_idx):,} cells")

# Load into memory (subset)
print("Loading subset into memory...", flush=True)
adata_sub = adata[keep_idx].to_memory()
adata.file.close()
del adata
gc.collect()

print(f"  In-memory shape: {adata_sub.shape}")

# ── compute DEG markers ────────────────────────────────────────────────
print("\nComputing Wilcoxon DEG markers (one-vs-rest)...", flush=True)
sc.tl.rank_genes_groups(
    adata_sub,
    groupby="celltype_level1",
    method="wilcoxon",
    key_added="rank_genes_level1",
)
print("  Done.")

# ── extract and save markers per cell type ─────────────────────────────
print("\nExtracting markers...", flush=True)
result = adata_sub.uns["rank_genes_level1"]
groups = result["names"].dtype.names

celltype_markers = {}  # celltype → DataFrame

for ct in groups:
    names = [str(x) for x in result["names"][ct]]
    scores = result["scores"][ct].astype(float)
    logfc = result["logfoldchanges"][ct].astype(float)
    pvals = result["pvals"][ct].astype(float)
    pvals_adj = result["pvals_adj"][ct].astype(float)

    df = pd.DataFrame({
        "names": names,
        "scores": scores,
        "logfoldchanges": logfc,
        "pvals": pvals,
        "pvals_adj": pvals_adj,
    })
    # Filter to significant upregulated markers
    df = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 0)].reset_index(drop=True)

    celltype_markers[ct] = df

    # Save CSV
    safe_name = ct.replace("/", "-").replace(" ", "_")
    csv_path = os.path.join(MARKER_DIR, f"level1_{safe_name}.csv")
    df.to_csv(csv_path, index=False)
    print(f"  {ct}: {len(df)} significant markers → {csv_path}")

del adata_sub
gc.collect()
