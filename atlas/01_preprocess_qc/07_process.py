#!/usr/bin/env python
"""
Atlas 01 — Step 07: build final integrated object (UMAP / Leiden)

PURPOSE
    Mirror the scANVI embedding (from the step-06 NPZ) into the full-gene-space
    raw object, preserve counts in layers["counts"], normalize_total(1e4) +
    log1p(base=2) the X, then compute the neighbour graph on X_scANVI, UMAP
    (min_dist=0.3) and Leiden (resolution=0.15). Writes the final atlas object.

INPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_raw.h5ad     (raw + labels; steps 02/05)
    output_root/01_preprocess_qc/integration/X_scANVI_hvg.npz        (embedding; step 06)

OUTPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_final.h5ad
        (full gene space, normalized/log1p X, layers["counts"], obsm["X_scANVI"],
         UMAP, Leiden, obs["celltype_pred"])

RUNTIME TIER
    heavy (high-memory CPU: one full-matrix copy + neighbours + UMAP on full atlas).

MANUSCRIPT ROLE
    Produces the integrated atlas object that downstream annotation/figures build
    on. See README for how this maps to the deposited hgsc_atlas_scanvi.h5ad.

NOTE
    AUTHORITATIVE original cluster script (atlas_05_process.py). Analytical logic
    preserved EXACTLY (normalize_total target_sum=1e4, log1p base=2,
    neighbors(use_rep="X_scANVI"), umap(min_dist=0.3), leiden(resolution=0.15)).
    Only the cluster paths were centralised.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# =========================
# PATHS (central config)
# =========================
IN_RAW  = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_raw.h5ad")
EMB_NPZ = path("output_root", "01_preprocess_qc", "integration", "X_scANVI_hvg.npz")
OUT_RAW = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_final.h5ad")
# =========================

import scanpy as sc
import numpy as np
import scipy.sparse as sp
import pandas as pd

print(f"🔹 Load RAW → {IN_RAW}", flush=True)
adata = sc.read_h5ad(IN_RAW)
print(f"   RAW cells: {adata.n_obs:,}, genes: {adata.n_vars:,}", flush=True)

# ---- Mirror X_scANVI from saved NPZ into RAW ----
print(f"🔹 Load scANVI embedding (NPZ) → {EMB_NPZ}", flush=True)
npz = np.load(EMB_NPZ, allow_pickle=False)
adata.obsm["X_scANVI"] = npz["X_scANVI"].astype(np.float32)

# ---- Create counts layer, then normalize/log-transform X ----
print("🔹 Preserve counts in layers['counts'] and normalize/log1p X", flush=True)
if "counts" not in adata.layers:
    # This is the one big copy; run on a high-mem CPU node.
    adata.layers["counts"] = adata.X.copy()

# Keep sparse & compact where possible
if not sp.isspmatrix_csr(adata.X):
    adata.X = sp.csr_matrix(adata.X)

# Normalize total to 1e4 and log1p base=2 on X
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata, base=2)

# ---- Neighbors/UMAP/Leiden on the scANVI embedding ----
print("🔹 Compute neighbors on X_scANVI", flush=True)
sc.pp.neighbors(adata, use_rep="X_scANVI")

print("🔹 UMAP (min_dist=0.3)", flush=True)
sc.tl.umap(adata, min_dist=0.3)

print("🔹 Leiden (resolution=0.15)", flush=True)
sc.tl.leiden(adata, resolution=0.15, key_added="leiden")

# Optional sanity print
print("Cluster sizes (top 20):")
print(adata.obs["leiden"].value_counts().head(20), flush=True)

print(f"🔹 Write final AnnData → {OUT_RAW}", flush=True)
adata.write_h5ad(OUT_RAW, compression="gzip")
print("✅ Done.", flush=True)
