#!/usr/bin/env python
"""
Atlas 01 — Step 05b: scANVI integration (the atlas integration of the manuscript)

PURPOSE
    Build scANVI from the trained scVI model (05a) using the CellAssign labels
    (04) as seed labels, fine-tune semi-supervised, and emit the latent
    embedding. This is the integration the manuscript uses (scVI → scANVI). Its
    embedding becomes obsm "X_scanvi" (lowercase) in the object written by
    07_finalize, which as_run/03b_refilter_umap.py reads.

INPUTS
    output_root/01_preprocess_qc/integration/anndata/preprocessed.h5ad   (step 03)
    output_root/01_preprocess_qc/integration/cellassign/predictions.csv  (step 04)
    output_root/01_preprocess_qc/integration/models/scvi/                (step 05a)

OUTPUTS
    output_root/01_preprocess_qc/integration/models/scanvi/
    output_root/01_preprocess_qc/integration/embeddings/scanvi/embedding.npz

RUNTIME TIER
    heavy — GPU. scANVI fine-tuning (max_epochs=800, early stopping). OPTIONAL re-run.

MANUSCRIPT ROLE
    The atlas integration. Its embedding (obsm["X_scanvi"]) underlies the atlas
    UMAP and clustering (Fig 1a; Figures S2, S3) and is the input to
    as_run/03b_refilter_umap.py.

NOTE
    AUTHORITATIVE official cluster script (03b_integrate_scanvi.py). Analytical
    params preserved EXACTLY (from_scvi_model, unlabeled_category="Unknown",
    max_epochs=800, early stop on elbo_validation, plan_kwargs lr=2e-4 /
    weight_decay=1e-4, gradient_clip_val=1.0). Only the cluster paths were
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
CELLASSIGN_PREDICTIONS = path("output_root", "01_preprocess_qc", "integration", "cellassign", "predictions.csv")
SCVI_MODEL_DIR = path("output_root", "01_preprocess_qc", "integration", "models", "scvi")
OUTPUT_DIR = path("output_root", "01_preprocess_qc", "integration")
MODEL_DIR = f"{OUTPUT_DIR}/models/scanvi"
EMBEDDING_NPZ = f"{OUTPUT_DIR}/embeddings/scanvi/embedding.npz"
LABEL_KEY = "celltype_pred"
UNLABELED_CATEGORY = "Unknown"
# =========================

import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import torch

scvi.settings.seed = SEED  # seed numpy/torch/scvi — best-effort determinism (GPU training
# is not guaranteed bit-reproducible; deposited integrated object is the trust boundary)
import os

# Create output directories
os.makedirs(f"{OUTPUT_DIR}/models/scanvi", exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/embeddings/scanvi", exist_ok=True)

# Set precision
torch.set_float32_matmul_precision("high")

print("=" * 60)
print("STEP 05B: scANVI INTEGRATION (SELECTED METHOD)")
print("=" * 60)

# Load preprocessed data
print(f"\n📂 Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes (HVGs): {adata.n_vars:,}")

# Load CellAssign predictions
print(f"\n📋 Loading CellAssign predictions: {CELLASSIGN_PREDICTIONS}")
predictions = pd.read_csv(CELLASSIGN_PREDICTIONS)
predictions = predictions.set_index('cell_id')

# Add labels to adata
labels = predictions['celltype_pred'].reindex(adata.obs_names)

# Handle missing labels
if labels.isna().any():
    print(f"   Warning: {labels.isna().sum()} cells without predictions, labeling as '{UNLABELED_CATEGORY}'")
    labels = labels.fillna(UNLABELED_CATEGORY)

# Convert to categorical
if UNLABELED_CATEGORY not in labels.values:
    labels = labels.astype('category')
else:
    labels = labels.astype('category')
    if UNLABELED_CATEGORY not in labels.cat.categories:
        labels = labels.cat.add_categories([UNLABELED_CATEGORY])

adata.obs[LABEL_KEY] = labels

print(f"\n📊 Label distribution:")
print(adata.obs[LABEL_KEY].value_counts())

# Load scVI model
print(f"\n🔄 Loading scVI model: {SCVI_MODEL_DIR}")
vae = scvi.model.SCVI.load(SCVI_MODEL_DIR, adata=adata)

# Convert to scANVI
print("\n⚙️  Building scANVI from scVI model...")
scanvi = scvi.model.SCANVI.from_scvi_model(
    vae,
    labels_key=LABEL_KEY,
    unlabeled_category=UNLABELED_CATEGORY
)

print(f"   Using labels from: {LABEL_KEY}")
print(f"   Unlabeled category: {UNLABELED_CATEGORY}")

# Train scANVI
print("\n🚀 Training scANVI (fine-tuning with labels)...")
scanvi.train(
    max_epochs=800,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    plan_kwargs={"lr": 2e-4, "weight_decay": 1e-4},
    gradient_clip_val=1.0
)

# Get latent representation
print("\n📊 Extracting latent representation...")
X_scanvi = scanvi.get_latent_representation()

# Save model
print(f"\n💾 Saving model: {MODEL_DIR}")
scanvi.save(MODEL_DIR, overwrite=True)

# Save embedding as NPZ
print(f"💾 Saving embedding: {EMBEDDING_NPZ}")
np.savez_compressed(
    EMBEDDING_NPZ,
    embedding=X_scanvi.astype(np.float32),
    obs_names=adata.obs_names.to_numpy()
)

print("\n✅ scANVI integration complete!")
print(f"   Model: {MODEL_DIR}")
print(f"   Embedding: {EMBEDDING_NPZ}")
