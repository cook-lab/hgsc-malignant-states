#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/31_fulltumour_nmf.py, the script that ran. Kept: the
# 31_program_* tables and kpca_fulltumour_programs.h5ad (read by 32, 33 and the Figure S10 script). Removed: the
# renaming of the program figures, which are no longer drawn. The untrimmed script is in the lab project folder.

"""
Step 31 — Unsupervised de-novo consensus NMF on the FULL epithelial subset
==========================================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Sarah's call: run a fresh unsupervised NMF on the ENTIRE epithelial (CNV-defined
tumour) subset — 12,429 cells, secretory + EMT — rather than projecting the old
secretory-only program. Reuses the exact step-04 cNMF machinery (HVG, k=4..12 x
20 cell-subsample replicates, stability-selected k, consensus centroids, per-cell
NNLS usage, cycling/ambient/stress flags) on a fresh tumour-only UMAP.

Input:   output/26_final_annotation/kpca_annotated.h5ad  (is_tumour)
Output:  output/31_fulltumour_nmf/
"""
import os, sys
from importlib import import_module
import numpy as np, pandas as pd, scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")
M04 = import_module("04_mouse_nmf")
OUT = os.path.join(S.OUTPUT, "31_fulltumour_nmf")


def main():
    os.makedirs(OUT, exist_ok=True)
    print("=" * 74)
    print("Step 31 — de-novo consensus NMF on the FULL epithelial subset")
    print("=" * 74)
    A = sc.read_h5ad(os.path.join(S.OUTPUT, "26_final_annotation", "kpca_annotated.h5ad"))
    E = A[A.obs["is_tumour"]].copy()
    print(f"    full epithelial/tumour subset: {E.n_obs:,} cells")

    # fresh tumour-only embedding for the program UMAPs
    print("[0] Fresh PCA / neighbours / UMAP on the tumour subset")
    Ee = E.copy()
    sc.pp.highly_variable_genes(Ee, n_top_genes=S.N_HVG)
    Ee = Ee[:, Ee.var["highly_variable"]].copy()
    sc.pp.scale(Ee, max_value=10)
    sc.tl.pca(Ee, n_comps=30, random_state=S.SEED)
    sc.pp.neighbors(Ee, n_neighbors=15, n_pcs=30, random_state=S.SEED)
    sc.tl.umap(Ee, random_state=S.SEED)
    E.obsm["X_umap_local"] = Ee.obsm["X_umap"]

    ortho = S.load_ortholog_map()
    X, genes = M04.get_hvg_matrix(E)
    per_k = M04.run_nmf_grid(X, S.NMF_KS, S.NMF_SEEDS)
    best_k, sel = M04.select_k(per_k, OUT)
    loadings, summ = M04.consensus_programs(per_k, best_k, genes)
    W, Wn = M04.program_usage(X, loadings)
    usage_n = pd.DataFrame(Wn, index=E.obs_names, columns=loadings.index)
    usage_r = pd.DataFrame(W, index=E.obs_names, columns=loadings.index)
    ann, cyc = M04.annotate_programs(E, loadings, usage_n, ortho, OUT)

    loadings.to_csv(os.path.join(OUT, "31_program_loadings.csv"))
    usage_n.to_csv(os.path.join(OUT, "31_program_usage.csv"))
    usage_r.to_csv(os.path.join(OUT, "31_program_usage_raw.csv"))
    ann.to_csv(os.path.join(OUT, "31_program_annotation.csv"), index=False)
    summ.to_csv(os.path.join(OUT, "31_program_summary.csv"), index=False)

    M04.make_figures(E, loadings, usage_n, ann, summ, OUT)
    E.write_h5ad(os.path.join(OUT, "kpca_fulltumour_programs.h5ad"))
    print(f"\n    selected k={best_k}; wrote {OUT}")
    print("Step 31 complete.")


if __name__ == "__main__":
    main()
