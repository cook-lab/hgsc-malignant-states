#!/usr/bin/env python
"""
Atlas 01 — Step 05a: scVI integration (base model for scANVI)

PURPOSE
    Train an scVI model on the HVG atlas (raw counts in layers["counts"], batched
    by sample_id) for batch correction; save the model and its latent embedding
    (NPZ). scANVI (05b) is built from this scVI model.

INPUTS
    output_root/01_preprocess_qc/integration/anndata/preprocessed.h5ad   (step 03)

OUTPUTS
    output_root/01_preprocess_qc/integration/models/scvi/
    output_root/01_preprocess_qc/integration/embeddings/scvi/embedding.npz

RUNTIME TIER
    heavy — GPU. scVI training (max_epochs=400, early stopping). OPTIONAL re-run.

MANUSCRIPT ROLE
    Methods integration (scVI → scANVI): the base model of the scANVI
    integration (05b). The scVI embedding itself is not read downstream.

NOTE
    AUTHORITATIVE official cluster script (03a_integrate_scvi.py). Analytical
    params preserved EXACTLY (n_latent=10, n_layers=2, gene_likelihood="nb",
    max_epochs=400, early_stopping_patience=10). Only the cluster paths were
    centralised; the cluster script's metadata JSON and QC UMAP preview are not
    included.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

# =========================
# CONFIG
# =========================
PREPROCESSED_H5AD = path("output_root", "01_preprocess_qc", "integration", "anndata", "preprocessed.h5ad")
OUTPUT_DIR = path("output_root", "01_preprocess_qc", "integration")
MODEL_DIR = f"{OUTPUT_DIR}/models/scvi"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/scvi/embedding.npz"
BATCH_KEY = "sample_id"
N_LATENT = 10
N_LAYERS = 2
# =========================

import scanpy as sc
import numpy as np
import scvi
import torch

scvi.settings.seed = SEED  # seed numpy/torch/scvi — best-effort determinism (GPU training
# is not guaranteed bit-reproducible; deposited integrated object is the trust boundary)
import os

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

print("\n✅ scVI integration complete!")
print(f"   Model: {MODEL_DIR}")
print(f"   Embedding: {EMBEDDING_NPZ}")
