#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/sandbox/06b_remap_level1.py, the script that ran.
# Kept: the corrected cluster-to-level-1 map, the cluster 8 smooth muscle / pericyte split and the artifact removal,
# written back to hgsc_atlas_celltype_level1.h5ad (input of level-2 annotation). Removed: the level-1 UMAP plot, summary
# prints and the assignment and summary CSVs, which nothing reads. The untrimmed script is in the lab project folder.
"""
Step 06b — Remap celltype_level1 on existing leiden_r0.3

The first run of 06 already computed leiden_r0.3 on the new atlas and saved it
to hgsc_atlas_celltype_level1.h5ad. But the cluster→label mapping was stale
(from the old atlas numbering). This script:
  1. Loads hgsc_atlas_celltype_level1.h5ad (already has leiden_r0.3)
  2. Applies the CORRECT mapping derived from the new crosstab
  3. Subclusters cluster 8 (SM/Pericyte) at r=0.5
  4. Removes artifact clusters 18/19 (38 cells)
  5. Saves CSVs, UMAP, and overwrites the h5ad

No Leiden recomputation needed — just remapping.
"""

import gc
import os
import warnings

import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

PROJECT_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD_PATH   = os.path.join(PROJECT_DIR, "hgsc_atlas_celltype_level1.h5ad")
OUT_H5AD    = H5AD_PATH  # overwrite in place

# ============================================================================
# CORRECT MAPPING (derived from new atlas crosstab)
# ============================================================================

LEIDEN_TO_LEVEL1 = {
    "0":  "Epithelial",       # 98.2% Epithelial (590,587 cells)
    "1":  "T/NK cell",        # 68%T + 31%NK (472,373 cells)
    "2":  "Macrophage",       # 81.0% Macrophage (379,256 cells)
    "3":  "Fibroblast",       # 99.1% Fibroblast (279,183 cells)
    "4":  "Epithelial",       # 97.7% Epithelial (85,012 cells)
    "5":  "T/NK cell",        # 88%T + 11%NK (83,963 cells)
    "6":  "Mesothelial",      # 54%F + 27%M — domain decision (64,123 cells)
    "7":  "Endothelial",      # 93.6% Endothelial (62,438 cells)
    "8":  "__SUBCLUSTER_SM_PERI__",  # 48%SM + 30%Peri + 21%Fib (61,840 cells)
    "9":  "T/NK cell",        # 77%T + 21%NK (58,987 cells)
    "10": "Fibroblast",       # 92.4% Fibroblast (50,151 cells)
    "11": "B cell",           # 94.9% B_cell (43,834 cells)
    "12": "Plasma cell",      # 97.4% Plasma_cell (34,748 cells)
    "13": "DC",               # 72.7% DC (12,811 cells)
    "14": "Epithelial",       # 99.3% Epithelial (6,127 cells)
    "15": "Mast cell",        # 97.0% Mast (4,911 cells)
    "16": "Neutrophil",       # 86.3% Neutrophil (4,459 cells)
    "17": "Macrophage",       # 100% Macrophage (51 cells — tiny)
    "18": "__ARTIFACT__",     # 33 cells, mixed — remove
    "19": "__ARTIFACT__",     # 5 cells — remove
}

ARTIFACT_CLUSTERS = ["18", "19"]

LEVEL1_ORDER = [
    "Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte",
    "Endothelial", "T/NK cell", "B cell", "Plasma cell", "Macrophage",
    "DC", "Neutrophil", "Mast cell",
]

# ============================================================================
# LOAD
# ============================================================================

print("=" * 60)
print("Step 06b — Remap celltype_level1 (no Leiden recompute)")
print("=" * 60)

print("\n1. Loading h5ad (already has leiden_r0.3)...", flush=True)
adata = sc.read_h5ad(H5AD_PATH)
print(f"   Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"   leiden_r0.3 clusters: {adata.obs['leiden_r0.3'].nunique()}")

# Verify leiden_r0.3 exists
assert "leiden_r0.3" in adata.obs.columns, "leiden_r0.3 not found!"
adata.obs["leiden_r0.3"] = adata.obs["leiden_r0.3"].astype(str)

# ============================================================================
# APPLY CORRECT MAPPING
# ============================================================================

print("\n2. Applying correct cluster → celltype_level1 mapping...", flush=True)
adata.obs["celltype_level1"] = adata.obs["leiden_r0.3"].map(LEIDEN_TO_LEVEL1)

# ============================================================================
# SUBCLUSTER CLUSTER 8 → SM / Pericyte
# ============================================================================

print("\n3. Subclustering cluster 8 to split SM vs Pericyte...", flush=True)

sm_peri_mask = adata.obs["leiden_r0.3"] == "8"
n_sm_peri = sm_peri_mask.sum()
print(f"   {n_sm_peri:,} cells in cluster 8")

sub = adata[sm_peri_mask].copy()

print("   Computing neighbors (n=15, X_scanvi)...", flush=True)
sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_scanvi", random_state=0)

print("   Running Leiden (r=0.5)...", flush=True)
sc.tl.leiden(sub, resolution=0.5, random_state=0)

n_sub = sub.obs["leiden"].nunique()
print(f"   → {n_sub} subclusters")

sub_cross = pd.crosstab(sub.obs["leiden"], sub.obs["celltype_pred"])
sub_cross.index = sub_cross.index.astype(int)
sub_cross = sub_cross.sort_index()
sub_totals = sub_cross.sum(axis=1)

print(f"\n   {'Sub':>5}  {'N cells':>9}  {'SM':>7}  {'Peri':>7}  {'Fibro':>7}  {'Assign':<15}")
print(f"   {'-'*5}  {'-'*9}  {'-'*7}  {'-'*7}  {'-'*7}  {'-'*15}")

sub_mapping = {}
for s in sub_cross.index:
    total = sub_totals[s]
    sm_n = sub_cross.loc[s, "Smooth_Muscle"] if "Smooth_Muscle" in sub_cross.columns else 0
    peri_n = sub_cross.loc[s, "Pericyte"] if "Pericyte" in sub_cross.columns else 0
    fib_n = sub_cross.loc[s, "Fibroblast"] if "Fibroblast" in sub_cross.columns else 0
    sm_pct = sm_n / total * 100
    peri_pct = peri_n / total * 100
    fib_pct = fib_n / total * 100

    if peri_pct > sm_pct:
        sub_mapping[s] = "Pericyte"
    else:
        sub_mapping[s] = "Smooth muscle"

    print(f"   {s:>5}  {int(total):>9,}  {sm_pct:>6.1f}%  {peri_pct:>6.1f}%  {fib_pct:>6.1f}%  {sub_mapping[s]}")

sub_labels = sub.obs["leiden"].astype(int).map(sub_mapping).values
adata.obs.loc[sm_peri_mask, "celltype_level1"] = sub_labels

sm_count = (sub_labels == "Smooth muscle").sum()
peri_count = (sub_labels == "Pericyte").sum()
print(f"\n   ✓ Smooth muscle: {sm_count:,}  Pericyte: {peri_count:,}")

del sub
gc.collect()

# ============================================================================
# REMOVE ARTIFACT CLUSTERS
# ============================================================================

artifact_mask = adata.obs["leiden_r0.3"].isin(ARTIFACT_CLUSTERS)
n_artifact = artifact_mask.sum()
print(f"\n4. Removing artifact clusters {ARTIFACT_CLUSTERS}: {n_artifact:,} cells")
adata = adata[~artifact_mask].copy()
print(f"   Remaining: {adata.shape[0]:,} cells")

# ============================================================================
# SAVE H5AD
# ============================================================================

print(f"\n8. Writing corrected h5ad...", flush=True)

adata.obs["celltype_level1"] = pd.Categorical(
    adata.obs["celltype_level1"].values, categories=LEVEL1_ORDER
)

print(f"   Saving: {OUT_H5AD}")
adata.write_h5ad(OUT_H5AD)
fsize = os.path.getsize(OUT_H5AD) / (1024**3)
print(f"   File size: {fsize:.1f} GB")

del adata
gc.collect()
