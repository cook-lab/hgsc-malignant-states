#!/usr/bin/env python3
"""
Epithelial NMF — SecA/SecB polarization programs
================================================
HGSC malignant-states atlas backend. Port of as_run/11d_epithelial_nmf.py, checked
against the deposited 11d run by 01b_verify_against_deposited.py.

Non-negative Matrix Factorization on epithelial cells to learn the gene
programs underlying the SecA (progenitor-like) -> SecB (adaptive/differentiated)
polarization axis. Trains on a stratified ~50,000-cell draw (49,997 cells; 3,000
seurat_v3 HVGs plus the forced SecA/SecB marker lists = 3,007 genes), runs NMF
(k=10) on log-normalized expression, then projects all epithelial cells onto the
learned factors. Factor_2 is the canonical SecB-defining factor used downstream
(as_run/17_prepare_nmf_labels.py and the figures' loader).

INPUTS:
  - <data_root>/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
    (epithelial subset; layers["counts"], obs celltype_level2)
  - SecA/SecB NMF marker lists (9 + 8, below): forced into the HVG set and used to
    label the learned factors, exactly as in the deposited 11d run. They extend the
    7 + 7 scoring signature of shared/signatures.yml with BCAM, WT1 (SecA) and MMP7 (SecB).

OUTPUTS (output_root/03_epithelial_nmf/):
  - 11d_nmf_usage.csv        normalized per-cell factor usage (Factor_2 defines the states)
  - 11d_nmf_loadings.csv     gene loadings (H matrix)
  - 11d_factor_scores.csv    SecA/SecB enrichment scores per factor (Table S5)
  - nmf_factor_mapping.json  dynamically-identified {"SecA":..., "SecB":...} factor
                             map (asserted downstream; reference SecB=Factor_2)

MANUSCRIPT PANELS: the deposited counterparts of these files give the epithelial
  states of every epithelial panel (Factor_2 usage; Figs 1-2, Figures S5-S9), the
  factor UMAPs of Fig 1h-j and Figure S5, Tables S5-S6, Figure S10 and, through
  the gene classes of as_run/11e_nmf_characterization.py, the TCGA-OV scores of Fig 7g-j.

RUNTIME TIER: heavy (loads full epithelial counts; NMF fit + projection).

SEEDING: NMF inits and the training subsample use SEED from config
  (np.random.seed(SEED) + RandomState(SEED)) for determinism.

Usage:
    python 01_epithelial_nmf.py
    python 01_epithelial_nmf.py --k 8
    python 01_epithelial_nmf.py --figures-only   # reuse the saved usage/loadings; rewrite only
                                                 # the factor scores and mapping (no figures)
"""

import argparse
import gc
import json
import os
import sys
import time
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse
from sklearn.decomposition import NMF

# config is 3 levels up: atlas/03_epithelial_nmf/01_epithelial_nmf.py -> repo root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

H5AD    = path("data_root", "2026_final_atlas", "celltype_h5ad",
               "hgsc_atlas_final_epithelial.h5ad")
OUT_DIR = path("output_root", "03_epithelial_nmf")

# ============================================================================
# NMF MARKER LISTS (deposited 11d run: 9 + 8)
# ============================================================================
# The deposited NMF (as_run/11d_epithelial_nmf.py) forced these 17 genes into the 3,000 HVGs
# (-> 3,007 features) and used them to label the factors; this port does the same, so a rerun
# reproduces the deposited model. The 7 + 7 signature in shared/signatures.yml stays the scoring
# signature everywhere else; it must remain a subset of these lists.
SECA_GENES = ["BCAM", "MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17", "WT1", "RCN2"]
SECB_GENES = ["KRT17", "KRT19", "KRT7", "LCN2", "MMP7", "PRSS22", "SLPI", "TACSTD2"]

_SIG_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "shared", "signatures.yml")
with open(os.path.abspath(_SIG_PATH)) as _fh:
    _SIGS = yaml.safe_load(_fh)
assert set(_SIGS["SecA"]) <= set(SECA_GENES) and set(_SIGS["SecB"]) <= set(SECB_GENES), \
    "shared/signatures.yml 7 + 7 must be a subset of the NMF marker lists"

# ============================================================================
# NMF PARAMETERS (validated — do not change)
# ============================================================================

N_TRAIN   = 50_000   # cells for NMF training
N_HVG     = 3000     # highly variable genes
DEFAULT_K = 10       # number of factors


# ============================================================================
# LOAD & PREPROCESS
# ============================================================================

def load_epithelial():
    """Load epithelial h5ad via backed mode, extract needed data."""
    print("  [LOAD] Opening epithelial h5ad (backed)...", flush=True)
    t0 = time.time()

    bdata = ad.read_h5ad(H5AD)  # in-memory: anndata 0.12.2 backed-sparse .X/layers read is broken
    n_cells = bdata.n_obs
    genes = bdata.var_names.tolist()

    level2 = bdata.obs["celltype_level2"].values.astype(str)
    barcodes = bdata.obs.index.values.copy()

    print(f"    Extracting counts matrix ({n_cells:,} x {len(genes):,})...", flush=True)
    counts = bdata.layers["counts"]

    del bdata
    gc.collect()

    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"celltype_level2": level2}, index=barcodes),
        var=pd.DataFrame(index=genes),
    )

    print(f"    {n_cells:,} cells x {len(genes):,} genes loaded in {time.time()-t0:.0f}s")
    return adata


def preprocess_for_nmf(adata):
    """Normalize, log-transform, select HVGs. Returns log-normalized adata."""
    print("  [PREPROCESS] Normalize + HVG selection...", flush=True)
    t0 = time.time()

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

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
        print(f"    Added {n_hvg_final - n_hvg} SecA/SecB genes -> {n_hvg_final:,} total")

    print(f"    Done in {time.time()-t0:.0f}s")
    return adata


# ============================================================================
# NMF
# ============================================================================

def stratified_train_idx(level2):
    """Stratified training subsample (per celltype_level2, >= 100 cells each), seeded from config.

    Same draw as the deposited 11d run (49,997 cells): RandomState(SEED), strata in sorted order.
    """
    np.random.seed(SEED)
    rng = np.random.RandomState(SEED)
    train_idx = []
    subtypes = np.unique(level2)
    for st in subtypes:
        st_idx = np.where(level2 == st)[0]
        n_st = len(st_idx)
        n_sample = max(100, int(N_TRAIN * n_st / len(level2)))
        n_sample = min(n_sample, n_st)
        train_idx.extend(rng.choice(st_idx, size=n_sample, replace=False).tolist())
    return np.array(sorted(train_idx))


def make_nmf(k=DEFAULT_K):
    """The NMF model as fitted in the deposited 11d run."""
    return NMF(n_components=k, init="nndsvda", random_state=SEED, max_iter=500)


def run_nmf(adata, k=DEFAULT_K):
    """Downsample, train NMF, project all cells. Seeded from config."""
    print(f"  [NMF] k={k}, training on {N_TRAIN:,} cells...", flush=True)
    t0 = time.time()

    level2 = adata.obs["celltype_level2"].values

    # Stratified downsample for training
    train_idx = stratified_train_idx(level2)
    print(f"    Training set: {len(train_idx):,} cells")

    hvg_mask = adata.var["highly_variable"].values
    hvg_genes = adata.var_names[hvg_mask].tolist()

    X_train = adata[train_idx, :][:, hvg_mask].X
    if sparse.issparse(X_train):
        X_train = X_train.toarray()
    X_train = np.clip(X_train, 0, None)

    print(f"    Training matrix: {X_train.shape[0]:,} x {X_train.shape[1]:,}", flush=True)

    t_fit = time.time()
    model = make_nmf(k)
    model.fit_transform(X_train)
    H = model.components_  # (k, n_hvg)
    recon_err = model.reconstruction_err_
    print(f"    NMF fit in {time.time()-t_fit:.0f}s (recon error: {recon_err:.2f})")

    del X_train
    gc.collect()

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

    W_norm = W_all / (W_all.sum(axis=1, keepdims=True) + 1e-10)

    factor_names = [f"Factor_{i+1}" for i in range(k)]
    usage_df = pd.DataFrame(W_norm, index=adata.obs.index, columns=factor_names)
    usage_raw = pd.DataFrame(W_all, index=adata.obs.index, columns=factor_names)
    loadings_df = pd.DataFrame(H, index=factor_names, columns=hvg_genes)

    print(f"    NMF complete ({time.time()-t0:.0f}s)")
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
        ranked = loadings.sort_values(ascending=False)
        top100 = set(ranked.index[:100])
        top200 = set(ranked.index[:200])

        factor_max = loadings.max()
        secA_mean = loadings[loadings.index.isin(SECA_GENES)].mean() / factor_max
        secB_mean = loadings[loadings.index.isin(SECB_GENES)].mean() / factor_max

        scores.append({
            "factor": factor,
            "secA_top100": len(set(SECA_GENES) & top100),
            "secA_top200": len(set(SECA_GENES) & top200),
            "secB_top100": len(set(SECB_GENES) & top100),
            "secB_top200": len(set(SECB_GENES) & top200),
            "secA_loading": secA_mean,
            "secB_loading": secB_mean,
        })

    scores_df = pd.DataFrame(scores).set_index("factor")
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
    parser = argparse.ArgumentParser(description="Epithelial NMF — SecA/SecB polarization")
    parser.add_argument("--k", type=int, default=DEFAULT_K,
                        help=f"Number of NMF factors (default: {DEFAULT_K})")
    parser.add_argument("--figures-only", action="store_true",
                        help="Load saved NMF results and regenerate figures")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    t0 = time.time()
    print("=" * 65)
    print("  Epithelial NMF: SecA/SecB Polarization")
    print("=" * 65)

    usage_path = os.path.join(OUT_DIR, "11d_nmf_usage.csv")
    loadings_path = os.path.join(OUT_DIR, "11d_nmf_loadings.csv")

    if args.figures_only and os.path.exists(usage_path) and os.path.exists(loadings_path):
        print("\n  [LOAD] Reading saved NMF results...")
        loadings_df = pd.read_csv(loadings_path, index_col=0)
        k = len(loadings_df)
    else:
        adata = load_epithelial()
        adata = preprocess_for_nmf(adata)
        usage_df, usage_raw, loadings_df, model = run_nmf(adata, k=args.k)
        k = args.k

        print("\n  [SAVE] Writing NMF results...")
        usage_df.to_csv(usage_path)
        loadings_df.to_csv(loadings_path)
        print(f"    Saved to {OUT_DIR}")

        del adata
        gc.collect()

    scores_df, best_secA, best_secB = identify_polarization_factors(loadings_df)
    scores_df.to_csv(os.path.join(OUT_DIR, "11d_factor_scores.csv"))

    # Persist the dynamically-identified SecA/SecB factor mapping so downstream
    # steps can ASSERT that the
    # factor they reference by literal name still matches what was learned here.
    # Guards against a re-run silently reordering factors (NMF factor indices
    # are not stable across fits). Reference values: SecB=Factor_2, SecA=Factor_3.
    factor_mapping = {"SecA": str(best_secA), "SecB": str(best_secB)}
    mapping_path = os.path.join(OUT_DIR, "nmf_factor_mapping.json")
    with open(mapping_path, "w") as _fh:
        json.dump(factor_mapping, _fh, indent=2)
    print(f"    Saved factor mapping: nmf_factor_mapping.json {factor_mapping}")

    print(f"\n{'=' * 65}")
    print(f"  COMPLETE ({time.time() - t0:.0f}s)")
    print(f"  Output: {OUT_DIR}")
    print(f"  NMF: k={k}, best SecA={best_secA}, best SecB={best_secB}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
