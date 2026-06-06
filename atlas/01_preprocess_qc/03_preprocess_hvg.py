#!/usr/bin/env python
"""
Atlas 01 — Step 03: preprocess + HVG selection (seurat_v3)

PURPOSE
    Load the QC'd, Scrublet-filtered concatenated atlas, stash raw counts into
    layers["counts"], select the top 4000 highly variable genes (seurat_v3,
    batch-aware by sample_id), subset to them, then normalize_total(1e4)+log1p so
    X holds log-normalised values and layers["counts"] the raw counts (standard
    scvi-tools convention). Produces the HVG object used by all integration
    methods (steps 05a–05e) and the benchmark (step 06).

INPUTS
    data_root/2026_final_atlas/processed/atlas_concatenated_filtered.h5ad
        (raw counts; written by 02_concat_qc_doublets — the QC/Scrublet step)

OUTPUTS
    output_root/01_preprocess_qc/integration/anndata/preprocessed.h5ad
        (X = log-normalised, layers["counts"] = raw, HVG-subset)

RUNTIME TIER
    heavy (seurat_v3 HVG on full atlas; high-memory CPU).

MANUSCRIPT ROLE
    Pre-integration feature selection; no panel rendered directly.

NOTE
    AUTHORITATIVE official cluster script (01_preprocess.py). Analytical params
    preserved EXACTLY (n_top_genes=4000, flavor="seurat_v3", batch_key="sample_id",
    layer="counts"; mark HVGs subset=False then subset by mask; normalize_total
    target_sum=1e4; log1p). Only the hardcoded cluster paths were centralised.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# =========================
# CONFIG
# =========================
RAW_H5AD = path("data_root", "2026_final_atlas", "processed", "atlas_concatenated_filtered.h5ad")
PREPROCESSED_H5AD = path("output_root", "01_preprocess_qc", "integration", "anndata", "preprocessed.h5ad")
N_TOP_GENES = 4000
BATCH_KEY = "sample_id"
# =========================

import scanpy as sc

print("=" * 60)
print("STEP 03: PREPROCESSING")
print("=" * 60)

print(f"\n📂 Loading raw data: {RAW_H5AD}")
adata = sc.read_h5ad(RAW_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# Preserve raw counts in a layer
print("\n📋 Creating counts layer...")
adata.layers["counts"] = adata.X.copy()

# Select highly variable genes (batch-aware, uses counts layer)
print(f"\n🔍 Identifying {N_TOP_GENES} highly variable genes (batch-aware)...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_TOP_GENES,
    batch_key=BATCH_KEY,
    flavor="seurat_v3",
    layer="counts",
    subset=False  # Don't subset yet, just mark them
)

n_hvg = adata.var["highly_variable"].sum()
print(f"   Found {n_hvg} HVGs")

# Subset to HVGs
print("\n✂️  Subsetting to HVGs...")
adata = adata[:, adata.var["highly_variable"]].copy()
print(f"   Dimensions: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

# Normalize (standard scvi-tools convention: X=normalized, layers["counts"]=raw)
print("\n📊 Normalizing...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Save preprocessed data
print(f"\n💾 Saving preprocessed data: {PREPROCESSED_H5AD}")
adata.write_h5ad(PREPROCESSED_H5AD, compression="gzip")

print("\n✅ Preprocessing complete!")
print(f"   Output: {PREPROCESSED_H5AD}")
print(f"   adata.X = log-normalized counts")
print(f"   adata.layers['counts'] = raw counts")
