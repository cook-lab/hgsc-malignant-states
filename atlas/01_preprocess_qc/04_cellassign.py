#!/usr/bin/env python
"""
Atlas 01 — Step 04: CellAssign cell-type annotation

PURPOSE
    Run CellAssign on the full-gene-space QC'd raw atlas using the curated v3
    marker matrix (81 genes × 16 cell types). Computes library-size factors on
    the full count matrix, builds a marker-only AnnData, trains CellAssign, and
    writes hard predictions + per-cell probabilities (+ QC plots). These
    predictions seed scANVI (05b) and are the cell-type labels carried into the
    final object (07_finalize).

INPUTS
    data_root/2026_final_atlas/processed/atlas_concatenated_filtered.h5ad   (QC raw counts)
    shared/cellassign_markers_v3.csv   (81-gene × 16-type binary marker matrix)

OUTPUTS
    output_root/01_preprocess_qc/integration/cellassign/predictions.csv
    output_root/01_preprocess_qc/integration/cellassign/probabilities.csv
    output_root/01_preprocess_qc/integration/cellassign/qc_plots/*.png

RUNTIME TIER
    heavy — GPU. CellAssign training (~400 epochs, early stopping patience=50).

MANUSCRIPT ROLE
    Cell-type priors / annotations used by scANVI and the final atlas (Fig 1,
    SF2). 16 types incl. Other (81 marker genes).

NOTE
    AUTHORITATIVE official cluster script (02_cellassign.py). Analytical logic
    preserved EXACTLY (size factors on full .X, marker-column subset, signatures
    passed un-binarized, early_stopping_patience=50, CONFIDENCE_THRESHOLD=0.5).
    Only the hardcoded cluster paths were centralised; MARKER_FILE points at the
    repo's shared/cellassign_markers_v3.csv (identical to the staged copy).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

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
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/qc_plots", exist_ok=True)

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

# Save convergence plot
print("\n📈 Saving convergence plot...")
history = ca.history
fig, ax = plt.subplots(figsize=(8, 5))
for key in ['elbo_validation', 'train_loss']:
    if key in history.keys():
        vals = history[key].values.flatten()
        ax.plot(vals, label=key, alpha=0.7)
ax.set_xlabel('Epoch')
ax.set_ylabel('Loss')
ax.set_title('CellAssign Training Convergence')
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/qc_plots/training_convergence.png", dpi=150)
plt.close()

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

# Generate QC plots
print("\n📊 Generating QC plots...")

# 1. Cell type distribution
fig, ax = plt.subplots(figsize=(10, 6))
ct_counts = predictions_df['celltype_pred'].value_counts()
ct_counts.plot(kind='bar', ax=ax, color='steelblue')
ax.set_xlabel('Cell Type')
ax.set_ylabel('Number of Cells')
ax.set_title('CellAssign: Cell Type Distribution')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/qc_plots/celltype_distribution.png", dpi=150)
plt.close()

# 2. Confidence histogram
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(max_probs, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
ax.axvline(CONFIDENCE_THRESHOLD, color='red', linestyle='--', linewidth=2, label=f'Threshold ({CONFIDENCE_THRESHOLD})')
ax.set_xlabel('Max Probability')
ax.set_ylabel('Number of Cells')
ax.set_title('CellAssign: Prediction Confidence')
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/qc_plots/confidence_histogram.png", dpi=150)
plt.close()

# 3. Per-sample composition
if 'sample_id' in adata.obs.columns:
    sample_comp = pd.crosstab(
        adata.obs['sample_id'],
        predictions_df.set_index('cell_id')['celltype_pred'],
        normalize='index'
    ) * 100

    fig, ax = plt.subplots(figsize=(12, max(6, len(sample_comp) * 0.3)))
    sample_comp.plot(kind='barh', stacked=True, ax=ax, legend=True)
    ax.set_xlabel('Percentage')
    ax.set_ylabel('Sample ID')
    ax.set_title('CellAssign: Cell Type Composition per Sample')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/qc_plots/per_sample_composition.png", dpi=150, bbox_inches='tight')
    plt.close()

print("\n✅ CellAssign complete!")
print(f"   Predictions: {PREDICTIONS_CSV}")
print(f"   Probabilities: {PROBABILITIES_CSV}")
print(f"   QC plots: {OUTPUT_DIR}/qc_plots/")
