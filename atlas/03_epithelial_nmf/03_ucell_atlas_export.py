#!/usr/bin/env python3
"""
Export atlas secretory epithelial cells for UCell scoring
=========================================================
HGSC malignant-states atlas backend.

Subsets the epithelial h5ad to NMF-defined secretory cells
(SecA / Intermediate / SecB) and exports raw counts in Matrix Market format
for cross-platform UCell scoring in R (04_ucell_atlas_scoring.R).

INPUTS:
  - <data_root>/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
  - output_root/03_epithelial_nmf/celltype_nmf_mapping.csv  (from 02_prepare_nmf_labels.py)

OUTPUTS (output_root/03_epithelial_nmf/ucell_atlas/):
  - atlas_secretory_counts.mtx.gz   (sparse raw counts, genes x cells)
  - atlas_secretory_barcodes.tsv
  - atlas_secretory_genes.tsv
  - atlas_secretory_metadata.csv
  - atlas_gene_list.txt

MANUSCRIPT PANELS: feeds the atlas-side cross-platform UCell scoring that
  matches the xenium noBCAM signature (Fig 3B organoid comparison, SF11).

RUNTIME TIER: moderate (chunked backed read of secretory subset).

SEEDING: deterministic (no stochastic step).

Usage:
    python 03_ucell_atlas_export.py
"""

import os
import sys
import warnings
import gzip

import numpy as np
import pandas as pd
import anndata as ad
import scipy.io
import scipy.sparse

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
EPI_H5AD = path("data_root", "2026_final_atlas", "celltype_h5ad",
                "hgsc_atlas_final_epithelial.h5ad")
NMF_MAP  = path("output_root", "03_epithelial_nmf", "celltype_nmf_mapping.csv")
OUT_DIR  = path("output_root", "03_epithelial_nmf", "ucell_atlas")
os.makedirs(OUT_DIR, exist_ok=True)

SECRETORY_LABELS = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]

# ── 1. Load NMF labels & filter to secretory ─────────────────
print("=" * 60)
print("  Export atlas secretory cells for UCell")
print("=" * 60)

print("\n[1] Loading NMF label mapping...")
nmf_map = pd.read_csv(NMF_MAP, index_col=0)
print(f"    Total cells in mapping: {len(nmf_map):,}")

sec_mask = nmf_map["celltype_nmf"].isin(SECRETORY_LABELS)
sec_barcodes = nmf_map.index[sec_mask]
print(f"    Secretory cells (SecA + Intermediate + SecB): {len(sec_barcodes):,}")
for label in SECRETORY_LABELS:
    n = (nmf_map.loc[sec_mask, "celltype_nmf"] == label).sum()
    print(f"      {label}: {n:,}")

# ── 2. Load epithelial h5ad & subset ─────────────────────────
print("\n[2] Loading epithelial h5ad (backed mode)...")
adata = ad.read_h5ad(EPI_H5AD, backed="r")
print(f"    Shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

shared_barcodes = sec_barcodes.intersection(adata.obs_names)
print(f"    Matched secretory barcodes in h5ad: {len(shared_barcodes):,}")
if len(shared_barcodes) < len(sec_barcodes):
    print(f"    WARNING: {len(sec_barcodes) - len(shared_barcodes)} barcodes not found")

# ── 3. Extract raw counts ────────────────────────────────────
print("\n[3] Extracting raw counts for secretory cells...")
idx = np.where(adata.obs_names.isin(shared_barcodes))[0]
print(f"    Subsetting {len(idx):,} cells...")

chunk_size = 50000
chunks = []
for start in range(0, len(idx), chunk_size):
    end = min(start + chunk_size, len(idx))
    chunk = adata.X[idx[start:end], :].copy()
    chunks.append(scipy.sparse.csr_matrix(chunk))
    print(f"    Processed chunk {start:,}-{end:,}")

counts = scipy.sparse.vstack(chunks, format="csr")
barcodes = adata.obs_names[idx].tolist()
genes = adata.var_names.tolist()
print(f"    Final matrix: {counts.shape[0]:,} cells x {counts.shape[1]:,} genes")

# ── 4. Export ────────────────────────────────────────────────
print("\n[4] Exporting to Matrix Market format...")

counts_t = counts.T.tocsc()  # R ReadMtx expects genes x cells
mtx_path = os.path.join(OUT_DIR, "atlas_secretory_counts.mtx.gz")
with gzip.open(mtx_path, "wb") as f:
    scipy.io.mmwrite(f, counts_t)
print(f"    Saved: atlas_secretory_counts.mtx.gz "
      f"({counts_t.shape[0]} genes x {counts_t.shape[1]} cells)")

pd.DataFrame(barcodes).to_csv(os.path.join(OUT_DIR, "atlas_secretory_barcodes.tsv"),
                              sep="\t", header=False, index=False)
print(f"    Saved: atlas_secretory_barcodes.tsv ({len(barcodes):,} barcodes)")

pd.DataFrame(genes).to_csv(os.path.join(OUT_DIR, "atlas_secretory_genes.tsv"),
                           sep="\t", header=False, index=False)
print(f"    Saved: atlas_secretory_genes.tsv ({len(genes):,} genes)")

with open(os.path.join(OUT_DIR, "atlas_gene_list.txt"), "w") as f:
    f.write("\n".join(genes))
print("    Saved: atlas_gene_list.txt")

meta = nmf_map.loc[barcodes, ["celltype_level2", "celltype_nmf", "celltype_level1"]].copy()
epi_obs = adata.obs.loc[barcodes, :]
for col in ["patient_id", "study", "sample_id"]:
    if col in epi_obs.columns:
        meta[col] = epi_obs[col].values
meta.to_csv(os.path.join(OUT_DIR, "atlas_secretory_metadata.csv"))
print(f"    Saved: atlas_secretory_metadata.csv ({len(meta):,} rows)")

adata.file.close()

print(f"\n{'='*60}")
print("  Step complete!")
print(f"  Output: {OUT_DIR}")
print(f"{'='*60}")
