#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/15b_apply_doublet_cut.py, the script that ran. Kept:
# kpca_all_clean.h5ad (read by 18 and step 21) and kpca_epithelial_clean.h5ad. Removed: the per-cluster removal
# table, the doublet-cut figure and the residual double-positive check. The untrimmed script, which the docstring
# below describes, is in the lab project folder.

"""
Step 15b — Apply a manual doublet cut, then re-embed and re-annotate
=====================================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

The scrublet score does NOT discriminate doublets here (max 0.183, compressed by
ambient contamination), so a per-cell scrublet cutoff is useless. The real
doublet signal is CO-EXPRESSION of mutually exclusive lineages, at genuine level,
combined with high complexity (two cells = ~2x the RNA of one).

Manual rule applied per cell:
    doublet  IF  (>= MIN_LINEAGES lineage programs each z > Z_CUT)
             AND (n_genes >= GENES_CUT)
The complexity requirement is what separates a real doublet from an ambient
double-positive: ambient bleeds a contaminating gene in at trace level in a
normal-sized cell, which does not clear both thresholds at once.

Rule is deliberately checked against the scrublet ranking and the Epcam+/Ptprc+
rate before it is applied, and the removed cells are reported per cluster.

Input:   output/01_kpca_qc/kpca_all.h5ad
Output:  output/15_doublet_cleanup/
             kpca_all_clean.h5ad, kpca_epithelial_clean.h5ad

Usage:
    python scripts/15b_apply_doublet_cut.py
    python scripts/15b_apply_doublet_cut.py --z 1.0 --genes 2500 --min-lineages 2
"""

import argparse
import os
import sys
from importlib import import_module

import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")

OUT_DIR = os.path.join(S.OUTPUT, "15_doublet_cleanup")
ALL_H5AD = os.path.join(S.OUTPUT, "01_kpca_qc", "kpca_all.h5ad")

# Lineage anchors: specific, mutually exclusive. Span the doublet types seen
# (epithelial / endothelial / fibroblast / myeloid / lymphoid). A singlet is
# positive for one; a doublet for two.
LINEAGE_ANCHORS = {
    "Epithelial":  ["Epcam", "Krt8", "Krt18", "Cdh1"],
    "Endothelial": ["Pecam1", "Cdh5", "Egfl7", "Emcn"],
    "Fibroblast":  ["Col1a1", "Col1a2", "Dcn", "Pdgfra"],
    "Myeloid":     ["Lyz2", "C1qa", "Csf1r", "Adgre1"],
    "Lymphoid":    ["Cd3e", "Cd3d", "Trac", "Nkg7"],
}


def lineage_z(A):
    """Per-cell z-scored score for each lineage anchor set."""
    Z = {}
    for lin, gl in LINEAGE_ANCHORS.items():
        present = [g for g in gl if g in A.var_names]
        sc.tl.score_genes(A, present, score_name=f"_l_{lin}", random_state=S.SEED)
        v = A.obs[f"_l_{lin}"].to_numpy(dtype=float)
        Z[lin] = (v - v.mean()) / (v.std() + 1e-9)
    return pd.DataFrame(Z, index=A.obs_names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--z", type=float, default=1.0,
                    help="a lineage counts as 'on' when its z-score exceeds this")
    ap.add_argument("--genes", type=int, default=2500,
                    help="high-complexity gene-count floor for a doublet")
    ap.add_argument("--min-lineages", type=int, default=2,
                    help="min simultaneously-on lineages to call a doublet")
    args = ap.parse_args()
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 78)
    print("Step 15b — manual doublet cut")
    print(f"    rule: >= {args.min_lineages} lineages with z > {args.z}  "
          f"AND  n_genes >= {args.genes}")
    print("=" * 78)

    A = sc.read_h5ad(ALL_H5AD)
    A.obs["leiden"] = A.obs["leiden"].astype(str)

    Z = lineage_z(A)
    n_on = (Z > args.z).sum(axis=1)
    A.obs["lineages_on"] = n_on.values
    A.obs["is_doublet"] = ((n_on >= args.min_lineages) &
                           (A.obs["n_genes_by_counts"] >= args.genes)).values

    n_dbl = int(A.obs["is_doublet"].sum())
    print(f"\n    flagged {n_dbl:,} / {A.shape[0]:,} cells "
          f"({n_dbl/A.shape[0]:.1%}) as doublets")
    print(f"    (scDblFinder independently estimated ~8.5%; "
          f"scrublet flag was an implausible 30%)")

    # ---- apply, re-embed, re-annotate -------------------------------------
    print("\n[apply] removing flagged cells and re-embedding")
    C = A[~A.obs["is_doublet"]].copy()
    for k in list(C.obs.columns):
        if k.startswith("_l_"):
            del C.obs[k]
    sc.pp.highly_variable_genes(C, n_top_genes=S.N_HVG)
    sc.pp.pca(C, n_comps=50, svd_solver="arpack", random_state=S.SEED)
    sc.pp.neighbors(C, n_neighbors=15, n_pcs=40, random_state=S.SEED)
    sc.tl.umap(C, random_state=S.SEED)
    sc.tl.leiden(C, resolution=1.0, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False, random_state=S.SEED)
    print(f"    clean object: {C.shape[0]:,} cells, "
          f"{C.obs['leiden'].nunique()} clusters")

    # re-run the ambient-robust compartment annotation from step 01
    imp01 = import_module("01_load_qc_kpca")
    C = imp01.annotate_compartments(C, OUT_DIR)

    print("\n    compartment sizes (clean):")
    print(C.obs["compartment"].value_counts().to_string())

    C.write_h5ad(os.path.join(OUT_DIR, "kpca_all_clean.h5ad"))
    print(f"\n    wrote {os.path.join(OUT_DIR, 'kpca_all_clean.h5ad')}")

    E = C[C.obs["compartment"] == "Epithelial"].copy()
    sc.pp.highly_variable_genes(E, n_top_genes=S.N_HVG)
    sc.pp.pca(E, n_comps=50, svd_solver="arpack", random_state=S.SEED)
    sc.pp.neighbors(E, n_neighbors=15, n_pcs=30, random_state=S.SEED)
    sc.tl.umap(E, random_state=S.SEED)
    E.obsm["X_umap_local"] = E.obsm["X_umap"].copy()
    E.write_h5ad(os.path.join(OUT_DIR, "kpca_epithelial_clean.h5ad"))
    print(f"    epithelium (clean): {E.shape[0]:,} cells -> "
          f"{os.path.join(OUT_DIR, 'kpca_epithelial_clean.h5ad')}")
    print("\nStep 15b complete.")


if __name__ == "__main__":
    main()
