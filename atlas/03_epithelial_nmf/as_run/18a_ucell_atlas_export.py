#!/usr/bin/env python3
"""
Step 18a — Export atlas secretory epithelial cells for UCell scoring
====================================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Subsets the epithelial h5ad to NMF-defined secretory cells (SecA, Transitioning,
SecB) and exports raw counts in Matrix Market format for UCell scoring in R.

Input:
    - celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
    - output/17_nmf_labels/celltype_nmf_mapping.csv

Output (output/18_ucell_atlas/):
    - atlas_secretory_counts.mtx.gz   (sparse raw counts)
    - atlas_secretory_barcodes.tsv
    - atlas_secretory_genes.tsv
    - atlas_secretory_metadata.csv
    - atlas_gene_list.txt             (full gene list for intersection in R)

Usage:
    python scripts/18a_ucell_atlas_export.py
"""

import os
import warnings

import numpy as np
import pandas as pd
import anndata as ad
import scipy.io
import scipy.sparse
import gzip

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
PROJECT_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas"
EPI_H5AD    = os.path.join(PROJECT_DIR, "celltype_h5ad",
                            "hgsc_atlas_final_epithelial.h5ad")
NMF_MAP     = os.path.join(PROJECT_DIR, "output", "17_nmf_labels",
                            "celltype_nmf_mapping.csv")
OUT_DIR     = os.path.join(PROJECT_DIR, "output", "18_ucell_atlas")
os.makedirs(OUT_DIR, exist_ok=True)

SECRETORY_LABELS = ["SecA epithelium", "Transitioning epithelium", "SecB epithelium"]

# ============================================================================
# 1. LOAD NMF LABELS & FILTER TO SECRETORY
# ============================================================================

print("=" * 60)
print("  Step 18a — Export atlas secretory cells for UCell")
print("=" * 60)

print("\n[1] Loading NMF label mapping...")
nmf_map = pd.read_csv(NMF_MAP, index_col=0)
print(f"    Total cells in mapping: {len(nmf_map):,}")

sec_mask = nmf_map["celltype_nmf"].isin(SECRETORY_LABELS)
sec_barcodes = nmf_map.index[sec_mask]
print(f"    Secretory cells (SecA + Transitioning + SecB): {len(sec_barcodes):,}")

for label in SECRETORY_LABELS:
    n = (nmf_map.loc[sec_mask, "celltype_nmf"] == label).sum()
    print(f"      {label}: {n:,}")

# ============================================================================
# 2. LOAD EPITHELIAL H5AD & SUBSET
# ============================================================================

print("\n[2] Loading epithelial h5ad (backed mode)...")
adata = ad.read_h5ad(EPI_H5AD, backed="r")
print(f"    Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# Match barcodes between NMF map and h5ad
shared_barcodes = sec_barcodes.intersection(adata.obs_names)
print(f"    Matched secretory barcodes in h5ad: {len(shared_barcodes):,}")

if len(shared_barcodes) < len(sec_barcodes):
    diff = len(sec_barcodes) - len(shared_barcodes)
    print(f"    WARNING: {diff} secretory barcodes not found in epithelial h5ad")

# ============================================================================
# 3. EXTRACT RAW COUNTS
# ============================================================================

print("\n[3] Extracting raw counts for secretory cells...")

# Get indices for subsetting
idx = np.where(adata.obs_names.isin(shared_barcodes))[0]
print(f"    Subsetting {len(idx):,} cells...")

# Read in chunks to manage memory for backed h5ad
chunk_size = 50000
chunks = []
for start in range(0, len(idx), chunk_size):
    end = min(start + chunk_size, len(idx))
    chunk_idx = idx[start:end]
    chunk = adata.X[chunk_idx, :].copy()
    if scipy.sparse.issparse(chunk):
        chunks.append(scipy.sparse.csr_matrix(chunk))
    else:
        chunks.append(scipy.sparse.csr_matrix(chunk))
    print(f"    Processed chunk {start:,}–{end:,}")

counts = scipy.sparse.vstack(chunks, format="csr")
barcodes = adata.obs_names[idx].tolist()
genes = adata.var_names.tolist()
print(f"    Final matrix: {counts.shape[0]:,} cells × {counts.shape[1]:,} genes")

# ============================================================================
# 4. EXPORT
# ============================================================================

print("\n[4] Exporting to Matrix Market format...")

# Counts (transpose for R: genes × cells is what ReadMtx expects)
counts_t = counts.T.tocsc()
mtx_path = os.path.join(OUT_DIR, "atlas_secretory_counts.mtx.gz")
with gzip.open(mtx_path, "wb") as f:
    scipy.io.mmwrite(f, counts_t)
print(f"    Saved: atlas_secretory_counts.mtx.gz ({counts_t.shape[0]} genes × {counts_t.shape[1]} cells)")

# Barcodes
bc_path = os.path.join(OUT_DIR, "atlas_secretory_barcodes.tsv")
pd.DataFrame(barcodes).to_csv(bc_path, sep="\t", header=False, index=False)
print(f"    Saved: atlas_secretory_barcodes.tsv ({len(barcodes):,} barcodes)")

# Genes
gene_path = os.path.join(OUT_DIR, "atlas_secretory_genes.tsv")
pd.DataFrame(genes).to_csv(gene_path, sep="\t", header=False, index=False)
print(f"    Saved: atlas_secretory_genes.tsv ({len(genes):,} genes)")

# Full gene list (for intersection computation in R)
gene_list_path = os.path.join(OUT_DIR, "atlas_gene_list.txt")
with open(gene_list_path, "w") as f:
    f.write("\n".join(genes))
print(f"    Saved: atlas_gene_list.txt")

# Metadata
meta = nmf_map.loc[barcodes, ["celltype_level2", "celltype_nmf", "celltype_level1"]].copy()

# Add patient/study from h5ad obs
epi_obs = adata.obs.loc[barcodes, :]
for col in ["patient_id", "study", "sample_id"]:
    if col in epi_obs.columns:
        meta[col] = epi_obs[col].values

meta_path = os.path.join(OUT_DIR, "atlas_secretory_metadata.csv")
meta.to_csv(meta_path)
print(f"    Saved: atlas_secretory_metadata.csv ({len(meta):,} rows)")

# Close backed reader
adata.file.close()

print(f"\n{'='*60}")
print(f"  Step 18a complete!")
print(f"  Output: {OUT_DIR}")
print(f"{'='*60}")
