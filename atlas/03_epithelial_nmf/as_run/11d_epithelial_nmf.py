#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/11d_epithelial_nmf.py, the script that ran.
# Kept: the k = 10 NMF fit and projection that write 11d_nmf_usage.csv, 11d_nmf_usage_raw.csv, 11d_nmf_loadings.csv and
# 11d_factor_scores.csv (states of every epithelial panel; Tables S5-S6). Removed: the five figure functions and their
# calls, style and palette. The untrimmed script is in the lab project folder.
"""
Step 11d — Epithelial NMF: SecA/SecB Polarization Programs
============================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Non-negative Matrix Factorization on epithelial cells to identify
gene programs underlying the SecA (progenitor-like) → SecB (adaptive/
differentiated) polarization axis.

Strategy:
  1. Downsample 50k cells stratified by subtype for NMF training
  2. Select top 3000 HVGs from training set
  3. Run NMF (k=10) on log-normalized expression
  4. Project all 575k cells onto learned factors
  5. Identify SecA/SecB-aligned factors via gene loading overlap
  6. Generate publication figures

SecA genes (progenitor-like):
  BCAM, MECOM, FBXO21, LGR5, LPAR3, PBX1, SOX17, WT1, RCN2

SecB genes (adaptive/differentiated):
  KRT17, KRT19, KRT7, LCN2, MMP7, PRSS22, SLPI, TACSTD2

Input:   celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
Output:  output/11d_epithelial_nmf/

Usage:
    python scripts/11d_epithelial_nmf.py
    python scripts/11d_epithelial_nmf.py --k 8
    python scripts/11d_epithelial_nmf.py --figures-only
"""

import argparse
import gc
import os
import time
import warnings

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import NMF

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

BASE     = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD     = os.path.join(BASE, "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
OUT_DIR  = os.path.join(BASE, "output", "11d_epithelial_nmf")

# ============================================================================
# GENE SIGNATURES
# ============================================================================

SECA_GENES = ["BCAM", "MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17", "WT1", "RCN2"]
SECB_GENES = ["KRT17", "KRT19", "KRT7", "LCN2", "MMP7", "PRSS22", "SLPI", "TACSTD2"]

# ============================================================================
# NMF PARAMETERS
# ============================================================================

N_TRAIN    = 50_000   # cells for NMF training
N_HVG      = 3000     # highly variable genes
DEFAULT_K  = 10       # number of factors
RNG_SEED   = 42


# ============================================================================
# LOAD & PREPROCESS
# ============================================================================

def load_epithelial():
    """Load epithelial h5ad via backed mode, extract needed data."""
    print("  [LOAD] Opening epithelial h5ad (backed)...", flush=True)
    t0 = time.time()

    bdata = ad.read_h5ad(H5AD, backed="r")
    n_cells = bdata.n_obs
    genes = bdata.var_names.tolist()

    # Extract metadata
    level2 = bdata.obs["celltype_level2"].values.astype(str)
    barcodes = bdata.obs.index.values.copy()
    umap_local = np.array(bdata.obsm["X_umap_local"])

    # Extract counts (sparse) — read in chunks for memory safety
    print(f"    Extracting counts matrix ({n_cells:,} × {len(genes):,})...", flush=True)
    counts = bdata.layers["counts"][:, :]  # reads full sparse matrix

    bdata.file.close()
    del bdata
    gc.collect()

    # Build in-memory AnnData with counts
    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"celltype_level2": level2}, index=barcodes),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["X_umap_local"] = umap_local

    print(f"    {n_cells:,} cells × {len(genes):,} genes loaded in {time.time()-t0:.0f}s")
    return adata


def preprocess_for_nmf(adata):
    """Normalize, log-transform, select HVGs. Returns log-normalized adata."""
    print("  [PREPROCESS] Normalize + HVG selection...", flush=True)
    t0 = time.time()

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # Normalize + log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3",
                                 layer="counts")
    n_hvg = adata.var["highly_variable"].sum()
    print(f"    {n_hvg:,} HVGs selected")

    # Ensure SecA/SecB genes are included
    for gene in SECA_GENES + SECB_GENES:
        if gene in adata.var_names:
            adata.var.loc[gene, "highly_variable"] = True

    n_hvg_final = adata.var["highly_variable"].sum()
    if n_hvg_final > n_hvg:
        print(f"    Added {n_hvg_final - n_hvg} SecA/SecB genes → {n_hvg_final:,} total")

    print(f"    Done in {time.time()-t0:.0f}s")
    return adata


# ============================================================================
# NMF
# ============================================================================

def run_nmf(adata, k=DEFAULT_K):
    """Downsample, train NMF, project all cells."""
    print(f"  [NMF] k={k}, training on {N_TRAIN:,} cells...", flush=True)
    t0 = time.time()

    rng = np.random.RandomState(RNG_SEED)
    level2 = adata.obs["celltype_level2"].values

    # Stratified downsample for training
    train_idx = []
    subtypes = np.unique(level2)
    for st in subtypes:
        st_idx = np.where(level2 == st)[0]
        n_st = len(st_idx)
        # Proportional allocation
        n_sample = max(100, int(N_TRAIN * n_st / len(level2)))
        n_sample = min(n_sample, n_st)
        train_idx.extend(rng.choice(st_idx, size=n_sample, replace=False).tolist())

    train_idx = np.array(sorted(train_idx))
    print(f"    Training set: {len(train_idx):,} cells")

    # Get HVG expression (log-normalized)
    hvg_mask = adata.var["highly_variable"].values
    hvg_genes = adata.var_names[hvg_mask].tolist()

    X_train = adata[train_idx, :][:, hvg_mask].X
    if sparse.issparse(X_train):
        X_train = X_train.toarray()

    # Ensure non-negative (log-normalized should be, but safety check)
    X_train = np.clip(X_train, 0, None)

    print(f"    Training matrix: {X_train.shape[0]:,} × {X_train.shape[1]:,}", flush=True)

    # Fit NMF
    t_fit = time.time()
    model = NMF(
        n_components=k,
        init="nndsvda",
        random_state=RNG_SEED,
        max_iter=500,
    )
    W_train = model.fit_transform(X_train)
    H = model.components_  # (k, n_hvg)
    recon_err = model.reconstruction_err_
    print(f"    NMF fit in {time.time()-t_fit:.0f}s (recon error: {recon_err:.2f})")

    del X_train
    gc.collect()

    # Project all cells
    print("    Projecting all cells...", flush=True)
    t_proj = time.time()
    X_all = adata[:, hvg_mask].X
    if sparse.issparse(X_all):
        X_all = X_all.toarray()
    X_all = np.clip(X_all, 0, None)

    W_all = model.transform(X_all)
    print(f"    Projection done in {time.time()-t_proj:.0f}s")

    del X_all
    gc.collect()

    # Normalize factor usage per cell (row-normalize to sum=1)
    W_norm = W_all / (W_all.sum(axis=1, keepdims=True) + 1e-10)

    # Build results
    factor_names = [f"Factor_{i+1}" for i in range(k)]
    usage_df = pd.DataFrame(W_norm, index=adata.obs.index, columns=factor_names)
    usage_raw = pd.DataFrame(W_all, index=adata.obs.index, columns=factor_names)

    loadings_df = pd.DataFrame(H, index=factor_names, columns=hvg_genes)

    elapsed = time.time() - t0
    print(f"    NMF complete ({elapsed:.0f}s)")

    return usage_df, usage_raw, loadings_df, model


# ============================================================================
# IDENTIFY SECA / SECB FACTORS
# ============================================================================

def identify_polarization_factors(loadings_df):
    """Score each factor for SecA and SecB gene enrichment."""
    print("  [IDENTIFY] Scoring factors for SecA/SecB enrichment...")

    scores = []
    for factor in loadings_df.index:
        loadings = loadings_df.loc[factor]

        # Rank genes by loading (descending)
        ranked = loadings.sort_values(ascending=False)
        top100 = set(ranked.index[:100])
        top200 = set(ranked.index[:200])

        # Count SecA/SecB genes in top loadings
        secA_in_top100 = len(set(SECA_GENES) & top100)
        secA_in_top200 = len(set(SECA_GENES) & top200)
        secB_in_top100 = len(set(SECB_GENES) & top100)
        secB_in_top200 = len(set(SECB_GENES) & top200)

        # Mean loading of signature genes (normalized by factor max)
        factor_max = loadings.max()
        secA_mean = loadings[loadings.index.isin(SECA_GENES)].mean() / factor_max
        secB_mean = loadings[loadings.index.isin(SECB_GENES)].mean() / factor_max

        scores.append({
            "factor": factor,
            "secA_top100": secA_in_top100,
            "secA_top200": secA_in_top200,
            "secB_top100": secB_in_top100,
            "secB_top200": secB_in_top200,
            "secA_loading": secA_mean,
            "secB_loading": secB_mean,
        })

    scores_df = pd.DataFrame(scores).set_index("factor")

    # Best SecA and SecB factors
    best_secA = scores_df["secA_loading"].idxmax()
    best_secB = scores_df["secB_loading"].idxmax()

    print(f"    Best SecA factor: {best_secA} "
          f"(loading={scores_df.loc[best_secA, 'secA_loading']:.3f}, "
          f"{scores_df.loc[best_secA, 'secA_top100']} genes in top100)")
    print(f"    Best SecB factor: {best_secB} "
          f"(loading={scores_df.loc[best_secB, 'secB_loading']:.3f}, "
          f"{scores_df.loc[best_secB, 'secB_top100']} genes in top100)")

    return scores_df, best_secA, best_secB


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="11d: Epithelial NMF — SecA/SecB polarization")
    parser.add_argument("--k", type=int, default=DEFAULT_K,
                        help=f"Number of NMF factors (default: {DEFAULT_K})")
    parser.add_argument("--figures-only", action="store_true",
                        help="Load saved NMF results and regenerate figures")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    t0 = time.time()
    print("=" * 65)
    print("  11d — Epithelial NMF: SecA/SecB Polarization")
    print("=" * 65)

    # Check for saved results
    usage_path = os.path.join(OUT_DIR, "11d_nmf_usage.csv")
    loadings_path = os.path.join(OUT_DIR, "11d_nmf_loadings.csv")

    if args.figures_only and os.path.exists(usage_path) and os.path.exists(loadings_path):
        print("\n  [LOAD] Reading saved NMF results...")
        usage_df = pd.read_csv(usage_path, index_col=0)
        loadings_df = pd.read_csv(loadings_path, index_col=0)
        k = len(loadings_df)

        # Still need UMAP coords and metadata
        print("  [LOAD] Extracting UMAP + metadata from h5ad...")
        with h5py.File(H5AD, "r") as f:
            umap_coords = np.array(f["obsm"]["X_umap_local"])
            level2 = pd.Series(
                [x.decode() if isinstance(x, bytes) else x
                 for x in f["obs"]["celltype_level2"][:]],
                index=[x.decode() if isinstance(x, bytes) else x
                       for x in f["obs"]["_index"][:]]
            )
    else:
        # Full pipeline
        adata = load_epithelial()
        adata = preprocess_for_nmf(adata)
        usage_df, usage_raw, loadings_df, model = run_nmf(adata, k=args.k)
        k = args.k

        umap_coords = adata.obsm["X_umap_local"]
        level2 = adata.obs["celltype_level2"]

        # Save results
        print("\n  [SAVE] Writing NMF results...")
        usage_df.to_csv(usage_path)
        usage_raw.to_csv(os.path.join(OUT_DIR, "11d_nmf_usage_raw.csv"))
        loadings_df.to_csv(loadings_path)
        print(f"    Saved to {OUT_DIR}")

        del adata
        gc.collect()

    # Identify SecA/SecB factors
    scores_df, best_secA, best_secB = identify_polarization_factors(loadings_df)
    scores_df.to_csv(os.path.join(OUT_DIR, "11d_factor_scores.csv"))

    # Summary
    elapsed = time.time() - t0
    print(f"\n{'=' * 65}")
    print(f"  11d COMPLETE ({elapsed:.0f}s)")
    print(f"  Output: {OUT_DIR}")
    print(f"  NMF: k={k}, best SecA={best_secA}, best SecB={best_secB}")
    print(f"  CSVs:")
    print(f"    11d_nmf_usage.csv        — Normalized factor usage (all cells)")
    print(f"    11d_nmf_usage_raw.csv    — Raw factor usage")
    print(f"    11d_nmf_loadings.csv     — Gene loadings (H matrix)")
    print(f"    11d_factor_scores.csv    — SecA/SecB enrichment scores")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
