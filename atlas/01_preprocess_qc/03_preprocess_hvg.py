#!/usr/bin/env python
"""
Atlas 01 — Step 03: preprocess + HVG selection (seurat_v3)

PURPOSE
    Read the raw atlas, stash raw counts into layers["counts"], select the top
    4000 highly variable genes (seurat_v3 flavour, batched by sample_id) and
    subset to them. Produces the HVG-subset object used to train scVI.

INPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_raw.h5ad   (from step 02)

OUTPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_preprocess.h5ad   (HVG-subset)

RUNTIME TIER
    heavy (seurat_v3 HVG on full atlas; high-memory CPU).

MANUSCRIPT ROLE
    Pre-integration feature selection; no panel rendered directly.

NOTE
    AUTHORITATIVE original cluster script (atlas_01_preprocess.py). Analytical
    logic preserved EXACTLY (n_top_genes=4000, flavor="seurat_v3",
    batch_key="sample_id", layer="counts", subset=True). Only the hardcoded
    cluster paths were replaced with central-config roots.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

import scanpy as sc, anndata as ad, gc, os, tempfile
import numpy as np
import scipy.sparse as sp

IN_RAW          = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_raw.h5ad")
OUT_PREPROCESS  = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_preprocess.h5ad")

adata = sc.read_h5ad(IN_RAW)

# ---------- Normalisation ----------
print("Preprocessing...", flush=True)
adata.layers["counts"] = adata.X.copy()

# ---------- HVGs ----------
print("HVGs...", flush=True)
sc.pp.highly_variable_genes(
    adata, n_top_genes=4000, batch_key="sample_id", flavor="seurat_v3", layer="counts", subset=True
)

# Save preprocessed data
adata.write_h5ad(
    OUT_PREPROCESS,
    compression="gzip"
)
