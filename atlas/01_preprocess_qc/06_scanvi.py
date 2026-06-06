#!/usr/bin/env python
"""
Atlas 01 — Step 06: scANVI label-propagation refinement

PURPOSE
    Load the scVI-integrated HVG object and the CellAssign labels (read from the
    raw object's obs), attach labels by cell ID, build scANVI from the saved scVI
    model, train it to propagate/refine labels, and store the scANVI latent
    representation in obsm["X_scANVI"]. Also writes a small NPZ of the embedding
    for a CPU-side merge into the raw object (step 07).

INPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_integrated.h5ad   (scVI HVG; step 04)
    output_root/01_preprocess_qc/integration/ovca_atlas_raw.h5ad          (CellAssign labels; step 05)
    output_root/01_preprocess_qc/integration/scvi_model_hgsc/             (scVI model; step 04)

OUTPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_integrated.h5ad   (in place: + obsm["X_scANVI"])
    output_root/01_preprocess_qc/integration/X_scANVI_hvg.npz             (embedding for CPU merge)

RUNTIME TIER
    heavy — GPU. scANVI training (max_epochs=800, early stopping).

MANUSCRIPT ROLE
    Final integration embedding (X_scANVI) used to build the atlas neighbour
    graph / UMAP / Leiden in step 07.

NOTE
    AUTHORITATIVE original cluster script (atlas_04b_scanvi.py). Analytical logic
    preserved EXACTLY (unlabeled_category="Unknown", from_scvi_model,
    max_epochs=800, early stopping on elbo_validation, plan_kwargs lr=2e-4 /
    weight_decay=1e-4, gradient_clip_val=1.0, NaN/Inf guard). Only the cluster
    paths were centralised.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# =========================
# PATHS (central config)
# =========================
IN_HVG     = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_integrated.h5ad")
IN_RAWLAB  = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_raw.h5ad")
MODEL_DIR  = path("output_root", "01_preprocess_qc", "integration", "scvi_model_hgsc")
OUT_HVG    = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_integrated.h5ad")
EMB_NPZ    = path("output_root", "01_preprocess_qc", "integration", "X_scANVI_hvg.npz")  # small file for CPU merge
LABEL_COL  = "celltype_pred"
UNLABELED  = "Unknown"
# =========================

import scanpy as sc
import numpy as np
import scvi
import gc
import pandas as pd

print(f"🔹  Read labels only from RAW (backed) → {IN_RAWLAB}", flush=True)
# This does NOT load counts; it just lets us read obs and then closes.
adata_raw_b = sc.read_h5ad(IN_RAWLAB, backed="r")
raw_labels = adata_raw_b.obs[LABEL_COL].copy()
raw_labels.index = adata_raw_b.obs_names.copy()
# Close and drop backed object immediately
try:
    adata_raw_b.file.close()
except Exception:
    pass
del adata_raw_b
gc.collect()

print(f"🔹  Load HVG → {IN_HVG}", flush=True)
adata = sc.read_h5ad(IN_HVG)

print("🔹 Attach labels to HVG by cell ID", flush=True)
# 🔹 Attach labels to HVG by cell ID (minimal & safe)
import pandas as pd  # make sure this is imported once near the top

labs = raw_labels.reindex(adata.obs_names)

# If categorical, add UNLABELED category if needed before fillna
if isinstance(labs.dtype, pd.CategoricalDtype):
    if UNLABELED not in labs.cat.categories:
        labs = labs.cat.add_categories([UNLABELED])
    labs = labs.fillna(UNLABELED)
else:
    # Non-categorical → fill then cast to category (UNLABELED will be included)
    labs = labs.fillna(UNLABELED).astype("category")

labs.name = LABEL_COL
adata.obs[LABEL_COL] = labs

print(labs.value_counts().head(20), flush=True)

print(f"🔹  Load SCVI model → {MODEL_DIR}", flush=True)
vae = scvi.model.SCVI.load(MODEL_DIR, adata=adata)

print(" Checking NaN")
import numpy as np, scipy.sparse as sp

# Use the same field SCVI was trained on
import numpy as np, scipy.sparse as sp

# SCVI was trained on counts; if that layer exists, check it directly.
print("Checking for NA")
X_for_scanvi = adata.layers["counts"] if "counts" in adata.layers else adata.X

if sp.issparse(X_for_scanvi):
    has_bad = np.isnan(X_for_scanvi.data).any() or np.isinf(X_for_scanvi.data).any()
else:
    has_bad = np.isnan(X_for_scanvi).any() or np.isinf(X_for_scanvi).any()

if has_bad:
    raise RuntimeError("Found NaN/Inf in the matrix used by SCANVI (counts/X).")

print("🔹  Build SCANVI from SCVI", flush=True)
scanvi = scvi.model.SCANVI.from_scvi_model(
    vae,
    labels_key=LABEL_COL,
    unlabeled_category=UNLABELED,
)

print("🔹  Train scANVI to propagate labels", flush=True)
# === YOUR ORIGINAL TRAINING CALL (unchanged) ===
scanvi.train(
    max_epochs=800,
    early_stopping=True,
    early_stopping_monitor="elbo_validation",
    plan_kwargs={"lr": 2e-4, "weight_decay": 1e-4},  # stability
    gradient_clip_val=1.0,                            # pass directly (not inside trainer_kwargs)
)
# ===============================================

print("🔹  Get scANVI latent", flush=True)
X_scanvi = scanvi.get_latent_representation()  # (n_cells, latent_dim) float32
adata.obsm["X_scANVI"] = X_scanvi

print(f"🔹  Overwrite HVG with X_scANVI → {OUT_HVG}", flush=True)
adata.write_h5ad(OUT_HVG, compression="gzip")

# Save tiny NPZ for CPU merge (no need to ever load RAW on GPU)
print(f"🔹  Save latent for CPU merge → {EMB_NPZ}", flush=True)
np.savez_compressed(
    EMB_NPZ,
    obs_names=adata.obs_names.to_numpy(),
    X_scANVI=X_scanvi.astype(np.float32)
)

# Cleanup
del vae, scanvi, adata, X_scanvi
gc.collect()
print("✅  SCANVI on GPU finished.", flush=True)
