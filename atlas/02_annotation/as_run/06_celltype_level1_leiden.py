#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/sandbox/06_celltype_level1_leiden.py, the script that ran.
# Kept: Leiden r = 0.3 and the level-1 labels written to hgsc_atlas_celltype_level1.h5ad (re-mapped by 06b), and the
# Leiden x CellAssign crosstab the level-1 map was read from. Removed: the level-1 UMAP plot, summary prints, and the
# assignment and summary CSVs, which 06b rewrote and nothing reads. The untrimmed script is in the lab project folder.
"""
Step 06 — Leiden r=0.3 Clustering + Cell Type Level 1 Annotation

Computes Leiden clustering at r=0.3 on the X_scanvi neighbor graph, then
assigns celltype_level1 directly from cluster composition of celltype_pred.
Cluster 8 (48% SM + 30% Pericyte + 21% Fibroblast) is subclustered at r=0.5
to separate Smooth muscle from Pericyte. Artifact clusters 18/19 (38 cells)
are removed.

Cluster → Label mapping (derived from new atlas crosstab):
  0, 4, 14               → Epithelial
  1, 5, 9                → T/NK cell
  2, 17                  → Macrophage
  3, 10                  → Fibroblast
  6                      → Mesothelial (domain decision: 54%F + 27%M)
  7                      → Endothelial
  8                      → subclustered into Smooth muscle / Pericyte
  11                     → B cell
  12                     → Plasma cell
  13                     → DC
  15                     → Mast cell
  16                     → Neutrophil
  18, 19                 → Artifact (38 cells total — removed)

Outputs:
  output/06_celltype_level1/leiden_r0.3_crosstab.csv
  output/06_celltype_level1/celltype_level1_assignments.csv
  output/06_celltype_level1/celltype_level1_summary.csv
  figs/06_celltype_level1/06_celltype_level1_umap.svg
  hgsc_atlas_celltype_level1.h5ad  (new file — original untouched)

Style: Cook Lab v1.2, small-point approach (pt_size=0.02, alpha=0.6)
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
H5AD_PATH   = os.path.join(PROJECT_DIR, "hgsc_atlas_scanvi.h5ad")
OUT_H5AD    = os.path.join(PROJECT_DIR, "hgsc_atlas_celltype_level1.h5ad")
OUT_DIR     = os.path.join(PROJECT_DIR, "output", "06_celltype_level1")

os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================
# LEIDEN → CELLTYPE_LEVEL1 MAPPING
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

# Clusters to remove (artifacts — too small for meaningful annotation)
ARTIFACT_CLUSTERS = ["18", "19"]

# Display order for legend (biological grouping)
LEVEL1_ORDER = [
    "Epithelial",
    "Mesothelial",
    "Fibroblast",
    "Smooth muscle",
    "Pericyte",
    "Endothelial",
    "T/NK cell",
    "B cell",
    "Plasma cell",
    "Macrophage",
    "DC",
    "Neutrophil",
    "Mast cell",
]

# ============================================================================
# LOAD ATLAS
# ============================================================================

print("=" * 60)
print("Step 06 — Leiden r=0.3 + Cell Type Level 1")
print("=" * 60)

print("\n1. Loading atlas (full — needed for Leiden + subclustering)...", flush=True)
adata = sc.read_h5ad(H5AD_PATH)
print(f"   Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# ============================================================================
# LEIDEN CLUSTERING r=0.3
# ============================================================================

print("\n2. Computing Leiden r=0.3...", flush=True)
sc.tl.leiden(adata, resolution=0.3, key_added="leiden_r0.3", random_state=0)
adata.obs["leiden_r0.3"] = adata.obs["leiden_r0.3"].astype(str)

n_clusters = adata.obs["leiden_r0.3"].nunique()
print(f"   → {n_clusters} clusters")

# Save crosstab
cross = pd.crosstab(adata.obs["leiden_r0.3"], adata.obs["celltype_pred"])
cross.index = cross.index.astype(int)
cross = cross.sort_index()
cross_path = os.path.join(OUT_DIR, "leiden_r0.3_crosstab.csv")
cross.to_csv(cross_path)
print(f"   Saved crosstab: {cross_path}")

# Print cluster composition
cross_pct = cross.div(cross.sum(axis=1), axis=0) * 100
print(f"\n   Cluster composition (top 2 types per cluster):")
for cl in cross.index:
    top2 = cross_pct.loc[cl].nlargest(2)
    total = cross.loc[cl].sum()
    desc = ", ".join([f"{t} {v:.0f}%" for t, v in top2.items()])
    print(f"   {cl:>3}  (n={total:>9,})  {desc}")

# ============================================================================
# MAP LEIDEN → CELLTYPE_LEVEL1 (initial pass)
# ============================================================================

print("\n3. Mapping leiden clusters → celltype_level1...", flush=True)
adata.obs["celltype_level1"] = adata.obs["leiden_r0.3"].map(LEIDEN_TO_LEVEL1)

# ============================================================================
# SUBCLUSTER CLUSTER 8 → Smooth muscle / Pericyte
# ============================================================================

SM_PERI_CLUSTERS = ["8"]
SUB_RESOLUTION = 0.5

print("\n4. Subclustering cluster 8 to split SM vs Pericyte...", flush=True)

sm_peri_mask = adata.obs["leiden_r0.3"].isin(SM_PERI_CLUSTERS)
n_sm_peri = sm_peri_mask.sum()
print(f"   {n_sm_peri:,} cells in cluster(s) {SM_PERI_CLUSTERS}")

sub = adata[sm_peri_mask].copy()

print("   Computing neighbors (n=15, X_scanvi)...", flush=True)
sc.pp.neighbors(sub, n_neighbors=15, use_rep="X_scanvi", random_state=0)

print(f"   Running Leiden (r={SUB_RESOLUTION})...", flush=True)
sc.tl.leiden(sub, resolution=SUB_RESOLUTION, random_state=0)

n_sub = sub.obs["leiden"].nunique()
print(f"   → {n_sub} subclusters")

# Crosstab to find SM-dominant vs Pericyte-dominant subclusters
sub_cross = pd.crosstab(sub.obs["leiden"], sub.obs["celltype_pred"])
sub_cross.index = sub_cross.index.astype(int)
sub_cross = sub_cross.sort_index()

sub_totals = sub_cross.sum(axis=1)

# For each subcluster, compute SM fraction vs Pericyte fraction
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

    # Assign based on which of SM or Pericyte is dominant
    if sm_pct >= peri_pct:
        sub_mapping[s] = "Smooth muscle"
    else:
        sub_mapping[s] = "Pericyte"

    print(f"   {s:>5}  {int(total):>9,}  {sm_pct:>6.1f}%  {peri_pct:>6.1f}%  {fib_pct:>6.1f}%  {sub_mapping[s]}")

# Write subcluster labels back
sub_labels = sub.obs["leiden"].astype(int).map(sub_mapping).values
adata.obs.loc[sm_peri_mask, "celltype_level1"] = sub_labels

print(f"\n   ✓ Assigned {n_sm_peri:,} cells to Smooth muscle / Pericyte")

del sub
gc.collect()

# ============================================================================
# REMOVE ARTIFACT CLUSTERS
# ============================================================================

artifact_mask = adata.obs["leiden_r0.3"].isin(ARTIFACT_CLUSTERS)
n_artifact = artifact_mask.sum()
print(f"\n5a. Removing artifact clusters {ARTIFACT_CLUSTERS}: {n_artifact:,} cells")
adata = adata[~artifact_mask].copy()
print(f"    Remaining: {adata.shape[0]:,} cells")

# ============================================================================
# HANDLE UNMAPPED
# ============================================================================

n_mapped = adata.obs["celltype_level1"].notna().sum()
n_sentinel = (adata.obs["celltype_level1"] == "__ARTIFACT__").sum()
n_missing = adata.obs["celltype_level1"].isna().sum() + n_sentinel
print(f"\n5b. Completeness check: Mapped {n_mapped - n_sentinel:,}  Missing/artifact: {n_missing:,}")

if n_missing > 0:
    print(f"   WARNING: {n_missing} cells unmapped — these should have been caught above")

# ============================================================================
# WRITE BACK TO H5AD
# ============================================================================

print(f"\n8. Writing celltype_level1 to new h5ad...", flush=True)

adata.obs["celltype_level1"] = pd.Categorical(
    adata.obs["celltype_level1"].values, categories=LEVEL1_ORDER
)

print(f"   Saving: {OUT_H5AD}")
adata.write_h5ad(OUT_H5AD)
fsize = os.path.getsize(OUT_H5AD) / (1024**3)
print(f"   File size: {fsize:.1f} GB")

del adata
gc.collect()
