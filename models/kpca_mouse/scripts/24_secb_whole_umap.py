#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/24_secb_whole_umap.py, the script that ran. Kept:
# 24_allcell_ucell.csv (read by 27) and the check against the step-22 scores. Removed: the per-compartment SecB
# summary and the whole-UMAP figure. The untrimmed script, which the docstring below describes, is in the lab
# project folder.

"""
Step 24 — SecA/SecB UCell across the ENTIRE UMAP (all compartments)
===================================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Specificity check: is SecB tumour-restricted, or does it bleed into immune /
stromal / NK cells? Scores every cell (not just tumour) with the SAME UCell
mechanism used to build the axis, then plots it over the whole-object UMAP and
quantifies SecB per compartment.

UCell is reimplemented in Python (rank-based, maxRank=1500 — the package default)
and validated against the saved tumour scores (output/22_secb_axis/22_ucell_scores.csv)
before being applied to all cells.

Input:   output/21_working/kpca_working.h5ad
Output:  output/24_whole_umap/
"""
import os, sys
from importlib import import_module
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")
OUT = os.path.join(S.OUTPUT, "24_whole_umap")
SEC = os.path.join(S.OUTPUT, "22_secb_axis")


def ucell(adata, genes, maxrank=1500):
    """UCell score per cell (Mann-Whitney U on within-cell gene ranks)."""
    genes = [g for g in genes if g in adata.var_names]
    n = len(genes)
    gi = [adata.var_names.get_loc(g) for g in genes]
    X = adata.X
    scores = np.empty(adata.n_obs, dtype=np.float32)
    step = 2000
    for start in range(0, adata.n_obs, step):
        end = min(start + step, adata.n_obs)
        blk = X[start:end]
        blk = blk.toarray() if hasattr(blk, "toarray") else np.asarray(blk)
        # rank descending: rank 1 = highest expressed gene in the cell
        order = np.argsort(-blk, axis=1, kind="stable")
        ranks = np.empty_like(order)
        rows = np.arange(blk.shape[0])[:, None]
        ranks[rows, order] = np.arange(1, blk.shape[1] + 1)[None, :]
        sig = ranks[:, gi].astype(np.float64)
        sig[sig > maxrank] = maxrank + 1          # cap
        U = sig.sum(axis=1)
        U_min = n * (n + 1) / 2.0
        scores[start:end] = 1.0 - (U - U_min) / (n * maxrank)
    return np.clip(scores, 0, 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    a = sc.read_h5ad(os.path.join(S.OUTPUT, "21_working", "kpca_working.h5ad"))
    xy = a.obsm["X_umap"]
    secb = pd.read_csv(os.path.join(SEC, "secb_genes.csv")).iloc[:, 0].tolist()
    seca = pd.read_csv(os.path.join(SEC, "seca_genes.csv")).iloc[:, 0].tolist()

    print("scoring UCell on all", a.n_obs, "cells ...")
    a.obs["SecB_UCell"] = ucell(a, secb)
    a.obs["SecA_UCell"] = ucell(a, seca)
    a.obs["polarization"] = a.obs["SecB_UCell"] - a.obs["SecA_UCell"]

    # ---- validate against saved tumour scores ----
    sv = pd.read_csv(os.path.join(SEC, "22_ucell_scores.csv"), index_col=0)
    shared = a.obs_names.intersection(sv.index)
    rho = spearmanr(a.obs.loc[shared, "SecB_UCell"], sv.loc[shared, "SecB_UCell"]).correlation
    rho_a = spearmanr(a.obs.loc[shared, "SecA_UCell"], sv.loc[shared, "SecA_UCell"]).correlation
    print(f"validation vs saved tumour UCell (n={len(shared)}): "
          f"SecB rho={rho:.3f}, SecA rho={rho_a:.3f}")

    a.obs[["SecA_UCell", "SecB_UCell", "polarization", "compartment"]].to_csv(
        os.path.join(OUT, "24_allcell_ucell.csv"))

    print("\nsaved", OUT)


if __name__ == "__main__":
    main()
