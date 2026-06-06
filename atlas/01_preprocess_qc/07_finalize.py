#!/usr/bin/env python
"""
Atlas 01 — Step 07: finalize integrated object

PURPOSE
    Assemble the final integrated atlas: load the QC'd raw counts, attach the
    CellAssign annotations (label / probability / confident), mirror the chosen
    integration embedding into obsm, normalize_total(1e4) + log1p(base=2) X
    (counts preserved in layers["counts"]), then compute neighbours / UMAP /
    Leiden on the embedding. The selected method is scANVI (--method scanvi),
    giving obsm["X_scanvi"] — the object consumed by 08_refilter_umap.

INPUTS
    data_root/2026_final_atlas/processed/atlas_concatenated_filtered.h5ad   (QC raw counts)
    output_root/01_preprocess_qc/integration/cellassign/predictions.csv     (step 04)
    output_root/01_preprocess_qc/integration/embeddings/<method>/embedding.npz (steps 05a–05e)

OUTPUTS
    output_root/01_preprocess_qc/integration/anndata/integrated_<method>.h5ad
        (X = lognorm base2, layers["counts"], obsm["X_<method>"], UMAP, Leiden,
         obs["celltype_pred"/"celltype_probability"/"celltype_confident"])

RUNTIME TIER
    heavy (full-matrix copy + neighbours + UMAP on full atlas; high-memory CPU).

MANUSCRIPT ROLE
    Produces the integrated atlas object. With --method scanvi this yields
    integrated_scanvi.h5ad (obsm["X_scanvi"]), the input to 08_refilter_umap which
    writes the deposited hgsc_atlas_scanvi.h5ad (obj("atlas_scanvi")).

NOTE
    AUTHORITATIVE official cluster script (05_finalize.py). Analytical params
    preserved EXACTLY (argparse --method/--leiden-res/--umap-min-dist;
    normalize_total target_sum=1e4; log1p base=2; UMAP min_dist default 0.3;
    Leiden resolution default 0.2). Only the cluster paths were centralised and
    the original's accidental duplicate sc.pp.neighbors call (idempotent) collapsed
    to one; corrupted emoji in print strings were cleaned (cosmetic only).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# =========================
# CONFIG
# =========================
RAW_H5AD = path("data_root", "2026_final_atlas", "processed", "atlas_concatenated_filtered.h5ad")
CELLASSIGN_PREDICTIONS = path("output_root", "01_preprocess_qc", "integration", "cellassign", "predictions.csv")
EMBEDDINGS_DIR = path("output_root", "01_preprocess_qc", "integration", "embeddings")
OUTPUT_DIR = path("output_root", "01_preprocess_qc", "integration", "anndata")
# Analysis parameters
UMAP_MIN_DIST = 0.3
LEIDEN_RESOLUTION = 0.2
# =========================

import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import argparse
import os

print("=" * 60)
print("STEP 07: FINALIZE INTEGRATED OBJECT")
print("=" * 60)

# Parse command line arguments
parser = argparse.ArgumentParser(description='Assemble final integrated object')
parser.add_argument('--method', type=str, default='scanvi',
                    choices=['scvi', 'scanvi', 'mrvi', 'sysvi', 'harmony'],
                    help='Integration method to use (default: scanvi, the selected method)')
parser.add_argument('--leiden-res', type=float, default=LEIDEN_RESOLUTION,
                    help=f'Leiden clustering resolution (default: {LEIDEN_RESOLUTION})')
parser.add_argument('--umap-min-dist', type=float, default=UMAP_MIN_DIST,
                    help=f'UMAP min_dist parameter (default: {UMAP_MIN_DIST})')

args = parser.parse_args()

METHOD = args.method
LEIDEN_RESOLUTION = args.leiden_res
UMAP_MIN_DIST = args.umap_min_dist
OUTPUT_H5AD = f"{OUTPUT_DIR}/integrated_{METHOD}.h5ad"

print(f"\nConfiguration:")
print(f"   Method: {METHOD}")
print(f"   Leiden resolution: {LEIDEN_RESOLUTION}")
print(f"   UMAP min_dist: {UMAP_MIN_DIST}")

# Load raw data (immutable base)
print(f"\nLoading raw data: {RAW_H5AD}")
adata = sc.read_h5ad(RAW_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")

# Add CellAssign annotations
print(f"\nLoading CellAssign annotations: {CELLASSIGN_PREDICTIONS}")
predictions = pd.read_csv(CELLASSIGN_PREDICTIONS, index_col='cell_id')
adata.obs['celltype_pred'] = predictions['celltype_pred'].reindex(adata.obs_names)
adata.obs['celltype_probability'] = predictions['max_probability'].reindex(adata.obs_names)
adata.obs['celltype_confident'] = predictions['confident'].reindex(adata.obs_names)

print(f"   Cell types: {adata.obs['celltype_pred'].nunique()}")
print(adata.obs['celltype_pred'].value_counts())

# Load selected embedding
embedding_path = f"{EMBEDDINGS_DIR}/{METHOD}/embedding.npz"
print(f"\nLoading {METHOD} embedding: {embedding_path}")

if not os.path.exists(embedding_path):
    raise FileNotFoundError(f"Embedding not found: {embedding_path}\nRun 05x_integrate_{METHOD}.py first!")

embedding_data = np.load(embedding_path, allow_pickle=True)
X_embed = embedding_data['embedding']
print(f"   Shape: {X_embed.shape}")

# Verify cell order matches
obs_names_embed = embedding_data['obs_names']
if not np.array_equal(obs_names_embed, adata.obs_names.to_numpy()):
    print("   Warning: Cell order mismatch, reindexing...")
    # Create temporary DataFrame for proper reindexing
    embed_df = pd.DataFrame(X_embed, index=obs_names_embed)
    X_embed = embed_df.reindex(adata.obs_names).values

adata.obsm[f'X_{METHOD}'] = X_embed

# Normalize X for visualization (preserve counts in layer)
print("\nNormalizing for visualization...")
print("   Creating counts layer...")
adata.layers['counts'] = adata.X.copy()

# Ensure sparse format
if not sp.issparse(adata.X):
    adata.X = sp.csr_matrix(adata.X)

print("   Normalizing to 10,000 counts per cell...")
sc.pp.normalize_total(adata, target_sum=1e4)

print("   Log-transforming (base 2)...")
sc.pp.log1p(adata, base=2)

# Compute neighbors on the embedding
print(f"\nComputing neighborhood graph on {METHOD} embedding...")
sc.pp.neighbors(adata, use_rep=f'X_{METHOD}')

# Compute UMAP
print(f"\nComputing UMAP (min_dist={UMAP_MIN_DIST})...")
sc.tl.umap(adata, min_dist=UMAP_MIN_DIST)

# Compute Leiden clustering
print(f"\nComputing Leiden clustering (resolution={LEIDEN_RESOLUTION})...")
sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, key_added='leiden')

n_clusters = adata.obs['leiden'].nunique()
print(f"   Found {n_clusters} clusters")
print("\n   Cluster sizes:")
print(adata.obs['leiden'].value_counts().sort_index())

# Save final integrated object
print(f"\nSaving integrated object: {OUTPUT_H5AD}")
adata.write_h5ad(OUTPUT_H5AD, compression="gzip")

# Summary
print("\n" + "=" * 60)
print("FINAL OBJECT SUMMARY")
print("=" * 60)
print(f"Created: {OUTPUT_H5AD}")
print(f"\nContents:")
print(f"   - X: Normalized, log-transformed counts ({adata.n_obs:,} x {adata.n_vars:,})")
print(f"   - layers['counts']: Raw integer counts")
print(f"   - obsm['X_{METHOD}']: Integration embedding ({X_embed.shape[1]}D)")
print(f"   - obsm['X_umap']: UMAP coordinates")
print(f"   - obs['celltype_pred']: CellAssign annotations ({adata.obs['celltype_pred'].nunique()} types)")
print(f"   - obs['leiden']: Leiden clusters ({n_clusters} clusters)")
print(f"   - obs['sample_id']: Batch information ({adata.obs['sample_id'].nunique()} samples)")

print("\nPipeline complete!")
