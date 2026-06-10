#!/usr/bin/env python
"""
Atlas 01 — Step 05a: scVI integration (method comparison)

PURPOSE
    Train an scVI model on the HVG atlas (raw counts in layers["counts"], batched
    by sample_id) for batch correction; save the model, the latent embedding
    (NPZ), metadata, and a QC UMAP preview. scVI is also the base model that
    scANVI (05b, the SELECTED method) is built from.

INPUTS
    output_root/01_preprocess_qc/integration/anndata/preprocessed.h5ad   (step 03)
    output_root/01_preprocess_qc/integration/cellassign/predictions.csv  (step 04; optional, for QC UMAP)

OUTPUTS
    output_root/01_preprocess_qc/integration/models/scvi/
    output_root/01_preprocess_qc/integration/embeddings/scvi/embedding.npz (+ metadata.json, umap_preview.png)

RUNTIME TIER
    heavy — GPU. scVI training (max_epochs=400, early stopping). OPTIONAL re-run.

MANUSCRIPT ROLE
    Integration-method comparison (one of 5 benchmarked; SF2/Methods). Base model
    for the selected scANVI integration (05b).

NOTE
    AUTHORITATIVE official cluster script (03a_integrate_scvi.py). Analytical
    params preserved EXACTLY (n_latent=10, n_layers=2, gene_likelihood="nb",
    max_epochs=400, early_stopping_patience=10). Only the cluster paths were
    centralised.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

# =========================
# CONFIG
# =========================
PREPROCESSED_H5AD = path("output_root", "01_preprocess_qc", "integration", "anndata", "preprocessed.h5ad")
CELLASSIGN_PREDICTIONS = path("output_root", "01_preprocess_qc", "integration", "cellassign", "predictions.csv")
OUTPUT_DIR = path("output_root", "01_preprocess_qc", "integration")
MODEL_DIR = f"{OUTPUT_DIR}/models/scvi"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/scvi/embedding.npz"
BATCH_KEY = "sample_id"
N_LATENT = 10
N_LAYERS = 2
# =========================

import scanpy as sc
import numpy as np
import pandas as pd
import scvi
import torch

scvi.settings.seed = SEED  # seed numpy/torch/scvi — best-effort determinism (GPU training
# is not guaranteed bit-reproducible; deposited integrated object is the trust boundary)
import json
import os
from datetime import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Create output directories
os.makedirs(f"{OUTPUT_DIR}/models/scvi", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/embeddings/scvi", exist_ok=True)

# Set precision for faster training
torch.set_float32_matmul_precision("high")

print("=" * 60)
print("STEP 05A: scVI INTEGRATION")
print("=" * 60)

# Load preprocessed data
print(f"\n📂 Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (HVGs): {adata.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# Load cell type predictions if available
has_celltypes = False
if os.path.exists(CELLASSIGN_PREDICTIONS):
    print(f"\n📋 Loading CellAssign predictions: {CELLASSIGN_PREDICTIONS}")
    predictions = pd.read_csv(CELLASSIGN_PREDICTIONS, index_col='cell_id')
    adata.obs['celltype_pred'] = predictions['celltype_pred'].reindex(adata.obs_names)
    has_celltypes = True
    print(f"   Cell types: {adata.obs['celltype_pred'].nunique()}")

# Setup scVI
print("\n⚙️  Setting up scVI model...")
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key=BATCH_KEY
)

vae = scvi.model.SCVI(
    adata,
    n_layers=N_LAYERS,
    n_latent=N_LATENT,
    gene_likelihood="nb"
)

print(f"   Latent dimensions: {N_LATENT}")
print(f"   Layers: {N_LAYERS}")
print(f"   Gene likelihood: negative binomial")

# Train model
print("\n🚀 Training scVI (this may take a while)...")
vae.train(
    max_epochs=400,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    early_stopping_patience=10
)

# Get latent representation
print("\n📊 Extracting latent representation...")
X_scvi = vae.get_latent_representation()

# Save model
print(f"\n💾 Saving model: {MODEL_DIR}")
vae.save(MODEL_DIR, overwrite=True)

# Save embedding as NPZ
print(f"💾 Saving embedding: {EMBEDDING_NPZ}")
np.savez_compressed(
    EMBEDDING_NPZ,
    embedding=X_scvi.astype(np.float32),
    obs_names=adata.obs_names.to_numpy()
)

# Save metadata
metadata = {
    "method": "scVI",
    "date": datetime.now().isoformat(),
    "n_cells": int(adata.n_obs),
    "n_genes": int(adata.n_vars),
    "n_batches": int(adata.obs[BATCH_KEY].nunique()),
    "n_latent": N_LATENT,
    "n_layers": N_LAYERS,
    "gene_likelihood": "nb",
    "scvi_tools_version": scvi.__version__
}

with open(f"{OUTPUT_DIR}/embeddings/scvi/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# Quick UMAP preview for QC
print("\n🗺️  Generating UMAP preview...")
adata.obsm["X_scvi"] = X_scvi
sc.pp.neighbors(adata, use_rep="X_scvi")
sc.tl.umap(adata)

if has_celltypes:
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(adata, color=BATCH_KEY, ax=axes[0], show=False, legend_loc='on data', legend_fontsize=6)
    axes[0].set_title("scVI: UMAP colored by batch")
    sc.pl.umap(adata, color='celltype_pred', ax=axes[1], show=False, legend_loc='right margin')
    axes[1].set_title("scVI: UMAP colored by cell type")
else:
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.umap(adata, color=BATCH_KEY, ax=ax, show=False, legend_loc='on data', legend_fontsize=6)
    ax.set_title("scVI: UMAP colored by batch")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/embeddings/scvi/umap_preview.png", dpi=150)
plt.close()

print("\n✅ scVI integration complete!")
print(f"   Model: {MODEL_DIR}")
print(f"   Embedding: {EMBEDDING_NPZ}")
print(f"   Preview: {OUTPUT_DIR}/embeddings/scvi/umap_preview.png")
