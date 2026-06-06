#!/usr/bin/env python
"""
Atlas 01 — Step 04: scVI training (batch correction)

PURPOSE
    Train an scVI model on the HVG-subset atlas (raw counts in layers["counts"],
    batched by sample_id) for multi-study batch correction, store the latent
    representation in obsm["X_scVI"], save the trained model, and write the
    integrated HVG object.

INPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_preprocess.h5ad   (from step 03)

OUTPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_integrated.h5ad   (obsm["X_scVI"])
    output_root/01_preprocess_qc/integration/scvi_model_hgsc/            (trained model)

RUNTIME TIER
    heavy — GPU. scVI training (max_epochs=800, early stopping). OPTIONAL re-run.

MANUSCRIPT ROLE
    Multi-study integration backbone (scVI batch correction). The scANVI
    refinement (step 06) reloads this model.

NOTE
    AUTHORITATIVE original cluster script (atlas_03_scvi.py). Model hyper-params
    preserved EXACTLY (n_layers=2, n_latent=30, gene_likelihood="nb",
    max_epochs=800, early_stopping on elbo_validation, patience=10). Only the
    hardcoded cluster paths were replaced with central-config roots. The original
    output filename had a typo ("ovca_atlas_integatrated.h5ad"); corrected here to
    "ovca_atlas_integrated.h5ad" so the chain is coherent (step 06 reads it).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

import scanpy as sc, anndata as ad, torch, scvi, gc

IN_PREPROCESS = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_preprocess.h5ad")
OUT           = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_integrated.h5ad")
MODEL_DIR     = path("output_root", "01_preprocess_qc", "integration", "scvi_model_hgsc")

torch.set_float32_matmul_precision("high")

print("Reading raw h5ad file …", flush=True)
adata = sc.read_h5ad(IN_PREPROCESS)

# ---------- scVI ----------
print("Setting up AnnData for scVI …", flush=True)
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",          # raw counts for the chosen HVGs
    batch_key="sample_id",
)

vae = scvi.model.SCVI(
    adata,
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb",
)

print("Training scVI (GPU) …", flush=True)
vae.train(
    max_epochs=800,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    early_stopping_patience=10
)

adata.obsm["X_scVI"] = vae.get_latent_representation()
print("Saving model and updated AnnData …", flush=True)
vae.save(MODEL_DIR, overwrite=True)
adata.write_h5ad(OUT, compression="gzip")

del adata, vae
gc.collect()
print("Done →", OUT)
