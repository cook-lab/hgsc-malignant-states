#!/usr/bin/env python
"""
Atlas 01 — Step 04: CellAssign cell-type annotation

PURPOSE
    Run CellAssign on the full-gene-space QC'd raw atlas using the curated v3
    marker matrix (81 genes × 16 cell types). Computes library-size factors on
    the full count matrix, builds a marker-only AnnData, trains CellAssign, and
    writes hard predictions + per-cell probabilities. These predictions seed
    scANVI (05b) and are the cell-type labels carried into the final object
    (07_finalize).

INPUTS
    data_root/2026_final_atlas/processed/atlas_concatenated_filtered.h5ad   (QC raw counts)
    shared/cellassign_markers_v3.csv   (81-gene × 16-type binary marker matrix)

OUTPUTS
    output_root/01_preprocess_qc/integration/cellassign/predictions.csv
    output_root/01_preprocess_qc/integration/cellassign/probabilities.csv

RUNTIME TIER
    heavy — GPU. CellAssign training (~400 epochs, early stopping patience=50).

MANUSCRIPT ROLE
    Seed labels for scANVI (05b). celltype_pred / max_probability / confident are
    carried into the final atlas as celltype_pred / celltype_probability /
    celltype_confident (07_finalize), which Figure S4a plots. 16 types incl.
    Other (81 marker genes).

NOTE
    AUTHORITATIVE official cluster script (02_cellassign.py). Analytical logic
    preserved EXACTLY (size factors on full .X, marker-column subset, signatures
    passed un-binarized, early_stopping_patience=50, CONFIDENCE_THRESHOLD=0.5).
    Only the hardcoded cluster paths were centralised; MARKER_FILE points at the
    repo's shared/cellassign_markers_v3.csv (identical to the staged copy). The
    cluster script's QC plots (training convergence, type distribution,
    confidence histogram, per-sample composition) are not included.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

# =========================
# CONFIG
# =========================
RAW_H5AD = path("data_root", "2026_final_atlas", "processed", "atlas_concatenated_filtered.h5ad")
MARKER_FILE = str(Path(__file__).resolve().parents[2] / "shared" / "cellassign_markers_v3.csv")
OUTPUT_DIR = path("output_root", "01_preprocess_qc", "integration", "cellassign")
PREDICTIONS_CSV = f"{OUTPUT_DIR}/predictions.csv"
PROBABILITIES_CSV = f"{OUTPUT_DIR}/probabilities.csv"
CONFIDENCE_THRESHOLD = 0.5  # For marking "confident" assignments
# =========================

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scvi.external import CellAssign
import scvi
import os

scvi.settings.seed = SEED  # seed numpy/torch/scvi for CellAssign training — best-effort
# (GPU training is not guaranteed bit-reproducible; the deposited integrated object is
# the trust boundary / entry point).

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("STEP 04: CELLASSIGN ANNOTATION")
print("=" * 60)

# Load raw data
print(f"\n📂 Loading raw data: {RAW_H5AD}")
adata = sc.read_h5ad(RAW_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")

# Ensure sparse format
if not sp.issparse(adata.X):
    print("   Converting to sparse CSR format...")
    adata.X = sp.csr_matrix(adata.X)

# Load marker gene definitions
print(f"\n📋 Loading marker genes: {MARKER_FILE}")
marker_gene_mat = pd.read_csv(MARKER_FILE, index_col=0)
print(f"   Marker genes: {len(marker_gene_mat)}")
print(f"   Cell types: {', '.join(marker_gene_mat.columns)}")

# Find intersection with dataset genes
marker_genes = marker_gene_mat.index.intersection(adata.var_names)
if len(marker_genes) == 0:
    raise ValueError("No marker genes found in dataset!")

print(f"   Markers found in dataset: {len(marker_genes)}/{len(marker_gene_mat)}")

missing = marker_gene_mat.index.difference(adata.var_names).tolist()
if missing:
    print(f"   Missing markers: {missing}")

# Subset to marker genes
marker_idx = np.where(adata.var_names.isin(marker_genes))[0]
marker_counts = adata.X[:, marker_idx]

print(f"\n🔬 Computing size factors...")
lib_size = np.ravel(adata.X.sum(axis=1))
size_factors = lib_size / lib_size.mean()
print(f"   Size factor mean: {size_factors.mean():.4f} (should be ~1.0)")

# Create lean AnnData for CellAssign
print("\n🧬 Building marker-only AnnData...")
bdata = ad.AnnData(
    X=marker_counts,
    obs=adata.obs.copy(),
    var=pd.DataFrame(index=adata.var_names[marker_idx])
)
bdata.obs["size_factor"] = size_factors

# Subset marker matrix to bdata genes - pass directly without binarization
signatures = marker_gene_mat.loc[bdata.var_names].copy()

print(f"   Signature matrix shape: {signatures.shape}")
print(f"   Markers per cell type:")
for ct in signatures.columns:
    n_markers = (signatures[ct] == 1).sum()
    print(f"      {ct}: {n_markers}")

# Setup and train CellAssign
print("\n🚀 Training CellAssign...")
CellAssign.setup_anndata(bdata, size_factor_key="size_factor")
ca = CellAssign(bdata, signatures)
ca.train(early_stopping_patience=50)  #400 epochs, early stopping. More patience to avoid stopping at bad local minimum)

print("\n📊 Predicting cell types...")
pred_probs = ca.predict()
pred_probs.index = bdata.obs_names  # Ensure index matches cells

# Get hard predictions
pred_labels = pred_probs.idxmax(axis=1)
max_probs = pred_probs.max(axis=1)
confident = max_probs >= CONFIDENCE_THRESHOLD

# Create predictions DataFrame
predictions_df = pd.DataFrame({
    'cell_id': pred_labels.index,
    'celltype_pred': pred_labels.values,
    'max_probability': max_probs.values,
    'confident': confident.values
})

# Save outputs
print(f"\n💾 Saving predictions: {PREDICTIONS_CSV}")
predictions_df.to_csv(PREDICTIONS_CSV, index=False)

print(f"💾 Saving probabilities: {PROBABILITIES_CSV}")
pred_probs.to_csv(PROBABILITIES_CSV)

# Summary statistics
print("\n📈 Cell Type Distribution:")
print(predictions_df['celltype_pred'].value_counts())
print(f"\nConfident assignments: {confident.sum():,} / {len(confident):,} ({100*confident.sum()/len(confident):.1f}%)")

print("\n✅ CellAssign complete!")
print(f"   Predictions: {PREDICTIONS_CSV}")
print(f"   Probabilities: {PROBABILITIES_CSV}")
