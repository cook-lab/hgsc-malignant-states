#!/usr/bin/env python
"""
Atlas 01 — Step 06: benchmark integration methods (scib-metrics)

PURPOSE
    Quantitatively compare the five integration embeddings (scVI, scANVI, MrVI,
    SysVI, Harmony) plus an unintegrated PCA baseline using the scib-metrics
    Benchmarker (bio-conservation + batch-correction metrics). Writes the scaled
    / raw metric tables and the summary plot used to justify selecting scANVI.

INPUTS
    output_root/01_preprocess_qc/integration/anndata/preprocessed.h5ad   (step 03)
    output_root/01_preprocess_qc/integration/cellassign/predictions.csv  (step 04)
    output_root/01_preprocess_qc/integration/embeddings/<method>/embedding.npz (steps 05a–05e)

OUTPUTS
    output_root/01_preprocess_qc/integration/benchmark/metrics_raw.csv, metrics_scaled.csv
    output_root/01_preprocess_qc/integration/benchmark/benchmark_summary.png (+ metadata.json, table)

RUNTIME TIER
    heavy (neighbour graphs per embedding; n_jobs parallel). CPU.

MANUSCRIPT ROLE
    Integration-method benchmark (SF2 / Methods). Basis for selecting scANVI (05b)
    as the manuscript integration.

NOTE
    AUTHORITATIVE official cluster script (04_benchmark.py). Analytical logic
    preserved EXACTLY (auto-detect embeddings, unintegrated PCA n_comps=30
    baseline, BioConservation(isolated_labels=False), default BatchCorrection,
    scIB 40/60 weighting). Only the cluster paths were centralised.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# =========================
# CONFIG
# =========================
PREPROCESSED_H5AD = path("output_root", "01_preprocess_qc", "integration", "anndata", "preprocessed.h5ad")
CELLASSIGN_PREDICTIONS = path("output_root", "01_preprocess_qc", "integration", "cellassign", "predictions.csv")
EMBEDDINGS_DIR = path("output_root", "01_preprocess_qc", "integration", "embeddings")
OUTPUT_DIR = path("output_root", "01_preprocess_qc", "integration", "benchmark")
BATCH_KEY = "sample_id"
LABEL_KEY = "celltype_pred"

# Methods to benchmark. Set to None to auto-detect from available embeddings.
METHODS = None  # or e.g. ["scvi", "harmony"]

N_JOBS = 32
# =========================

import scanpy as sc
import pandas as pd
import numpy as np
import os
import json
from datetime import datetime

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("STEP 06: BENCHMARK INTEGRATION METHODS (scib-metrics)")
print("=" * 60)

# -----------------------------
# Auto-detect methods
# -----------------------------
if METHODS is None:
    available = []
    if os.path.isdir(EMBEDDINGS_DIR):
        for name in os.listdir(EMBEDDINGS_DIR):
            npz_path = os.path.join(EMBEDDINGS_DIR, name, "embedding.npz")
            if os.path.exists(npz_path):
                available.append(name)
    METHODS = sorted(available)
    print(f"\n🔍 Auto-detected methods: {', '.join(METHODS)}")
else:
    print(f"\n📋 Specified methods: {', '.join(METHODS)}")

if len(METHODS) == 0:
    raise ValueError("No embeddings found! Run integration scripts first.")

# -----------------------------
# Load data
# -----------------------------
print(f"\n📂 Loading preprocessed data: {PREPROCESSED_H5AD}")
adata = sc.read_h5ad(PREPROCESSED_H5AD)
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")
print(f"   Batches: {adata.obs[BATCH_KEY].nunique()}")

# Load cell type labels
print(f"\n📋 Loading CellAssign predictions: {CELLASSIGN_PREDICTIONS}")
predictions = pd.read_csv(CELLASSIGN_PREDICTIONS, index_col='cell_id')
adata.obs[LABEL_KEY] = predictions['celltype_pred'].reindex(adata.obs_names)
print(f"   Cell types: {adata.obs[LABEL_KEY].nunique()}")

# -----------------------------
# Compute unintegrated PCA baseline
# -----------------------------
# preprocessed.h5ad is already log-normalized
print("\n🔧 Computing unintegrated PCA baseline...")
sc.tl.pca(adata, n_comps=30)
adata.obsm["X_pca_unintegrated"] = adata.obsm["X_pca"].copy()
print(f"   Shape: {adata.obsm['X_pca_unintegrated'].shape}")

# -----------------------------
# Load pre-computed embeddings
# -----------------------------
print("\n📊 Loading embeddings...")

embedding_keys = ["X_pca_unintegrated"]

for method in METHODS:
    npz_path = f"{EMBEDDINGS_DIR}/{method}/embedding.npz"
    print(f"   Loading {method}...")
    data = np.load(npz_path, allow_pickle=True)
    embedding = data['embedding']

    # Verify cell order matches
    if 'obs_names' in data:
        saved_names = data['obs_names']
        if not np.array_equal(saved_names, adata.obs_names.values):
            print(f"      Reordering to match adata...")
            idx = pd.Index(saved_names).get_indexer(adata.obs_names)
            embedding = embedding[idx]

    key = f"X_{method}"
    adata.obsm[key] = embedding
    embedding_keys.append(key)
    print(f"      Shape: {embedding.shape}")

print(f"\n✅ Methods to benchmark: {[k.replace('X_', '') for k in embedding_keys]}")

# -----------------------------
# Run scib-metrics Benchmarker
# -----------------------------
print("\n" + "=" * 60)
print("RUNNING SCIB-METRICS BENCHMARK")
print("=" * 60)

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

# Use defaults - disable isolated_labels for speed on large datasets
bio_metrics = BioConservation(isolated_labels=False)
batch_metrics = BatchCorrection()

bm = Benchmarker(
    adata,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    embedding_obsm_keys=embedding_keys,
    bio_conservation_metrics=bio_metrics,
    batch_correction_metrics=batch_metrics,
    pre_integrated_embedding_obsm_key="X_pca_unintegrated",
    n_jobs=N_JOBS
)

print("\n⚙️  Preparing (computing neighbors)...")
bm.prepare()

print("\n🚀 Computing metrics...")
bm.benchmark()

# -----------------------------
# Results
# -----------------------------
print("\n" + "=" * 60)
print("RESULTS")
print("=" * 60)

results_raw = bm.get_results(min_max_scale=False, clean_names=True)
results_scaled = bm.get_results(min_max_scale=True, clean_names=True)

print("\nScaled results (0-1):")
print(results_scaled.to_string())

# Save
results_raw.to_csv(f"{OUTPUT_DIR}/metrics_raw.csv")
results_scaled.to_csv(f"{OUTPUT_DIR}/metrics_scaled.csv")

metadata = {
    "date": datetime.now().isoformat(),
    "n_cells": int(adata.n_obs),
    "n_batches": int(adata.obs[BATCH_KEY].nunique()),
    "n_cell_types": int(adata.obs[LABEL_KEY].nunique()),
    "methods": embedding_keys,
    "batch_key": BATCH_KEY,
    "label_key": LABEL_KEY
}
with open(f"{OUTPUT_DIR}/metadata.json", "w") as f:
    json.dump(metadata, f, indent=2)

# -----------------------------
# Plots
# -----------------------------
print("\n📊 Generating plots...")

bm.plot_results_table(min_max_scale=True, show=False, save_dir=OUTPUT_DIR)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 6))
scores = results_scaled.T
scores = scores.drop('Metric Type', axis=1, errors='ignore')
agg_cols = [c for c in scores.columns if c in ['Batch correction', 'Bio conservation', 'Total']]
if agg_cols:
    scores[agg_cols].plot(kind='bar', ax=ax, width=0.8)
    ax.set_ylabel('Score (0-1)')
    ax.set_title('Integration Benchmark Summary\n(scIB weighting: 40% batch, 60% bio)')
    ax.set_xticklabels([k.replace('X_', '') for k in scores.index], rotation=45, ha='right')
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/benchmark_summary.png", dpi=150, bbox_inches='tight')
plt.close()

print("\n✅ Benchmarking complete!")
print(f"   Results: {OUTPUT_DIR}/metrics_scaled.csv")
print(f"   Plot: {OUTPUT_DIR}/benchmark_summary.png")
