#!/usr/bin/env python3
"""
Check that 01_epithelial_nmf.py reproduces the deposited 11d NMF
================================================================
PURPOSE
    The deposited NMF (data_root/2026_final_atlas/output/11d_epithelial_nmf/, written by
    as_run/11d_epithelial_nmf.py) is what the downstream tables and figures read
    (Factor_2 usage -> SecA/Intermediate/SecB; loadings -> Table S5 and the 177/177 TCGA lists). This
    check runs the repo port's own preprocessing, training draw and NMF fit, and compares them with the
    deposit:
      1. feature set: 3,000 seurat_v3 HVGs + the forced NMF marker lists (9 + 8) vs the 3,007 deposited
         loading columns; for comparison, also the feature set obtained by forcing only the 7 + 7 genes;
      2. training draw: the 49,997 stratified cells;
      3. loadings H (k = 10) vs 11d_nmf_loadings.csv, factor by factor;
      4. usage of the training cells (model.transform, row-normalized as in 01) vs 11d_nmf_usage.csv.
    It does not project all 575,366 cells (the heavy last step of 01); the training cells stand in.

INPUTS
    data_root/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad  (via 01's load_epithelial)
    data_root/2026_final_atlas/output/11d_epithelial_nmf/{11d_nmf_loadings.csv, 11d_nmf_usage.csv}

OUTPUTS
    output_root/03_epithelial_nmf/verify_deposited/verify_summary.json
    output_root/03_epithelial_nmf/verify_deposited/loadings_by_factor.csv
    output_root/03_epithelial_nmf/verify_deposited/feature_set_differences.csv

MANUSCRIPT PANEL(S): none (reproducibility check of the NMF behind Figs 1-2, Figures S5-S10,
    Tables S5-S6, Fig 7g-j)

RUNTIME TIER: heavy (loads the 27 GB epithelial h5ad in memory; NMF fit on 49,997 x 3,007).

USAGE
    DATA_ROOT=<deposit> OUTPUT_ROOT=<dir> python atlas/03_epithelial_nmf/01b_verify_against_deposited.py
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
from config.config import path  # noqa: E402

_spec = importlib.util.spec_from_file_location("nmf01", os.path.join(HERE, "01_epithelial_nmf.py"))
nmf01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nmf01)

DEPOSIT = path("data_root", "2026_final_atlas", "output", "11d_epithelial_nmf")
OUT = path("output_root", "03_epithelial_nmf", "verify_deposited")


def _r(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1])


def main():
    t0 = time.time()
    dep_H = pd.read_csv(os.path.join(DEPOSIT, "11d_nmf_loadings.csv"), index_col=0)
    print(f"Deposited loadings: {dep_H.shape[0]} factors x {dep_H.shape[1]:,} genes", flush=True)

    adata = nmf01.load_epithelial()
    adata = nmf01.preprocess_for_nmf(adata)          # forces the 9 + 8 lists into the HVGs
    var = adata.var
    feat_98 = set(var.index[var["highly_variable"].values])
    base = set(var.index[var["highly_variable_rank"].notna().values])
    assert len(base) == nmf01.N_HVG, f"expected {nmf01.N_HVG} seurat_v3 HVGs, got {len(base)}"
    sig77 = (set(nmf01._SIGS["SecA"]) | set(nmf01._SIGS["SecB"])) & set(var.index)
    feat_77 = base | sig77
    dep_set = set(dep_H.columns)
    diffs = ([("in_9+8_only", g) for g in sorted(feat_98 - dep_set)] +
             [("in_deposit_only_vs_9+8", g) for g in sorted(dep_set - feat_98)] +
             [("in_deposit_only_vs_7+7", g) for g in sorted(dep_set - feat_77)])
    pd.DataFrame(diffs, columns=["comparison", "gene"]).to_csv(
        os.path.join(OUT, "feature_set_differences.csv"), index=False)
    print(f"Features: 9+8 forcing {len(feat_98):,} (equal to deposit: {feat_98 == dep_set}); "
          f"7+7 forcing {len(feat_77):,} (equal: {feat_77 == dep_set})", flush=True)

    level2 = adata.obs["celltype_level2"].values
    train_idx = nmf01.stratified_train_idx(level2)
    hvg_mask = var["highly_variable"].values
    hvg_genes = var.index[hvg_mask].tolist()
    X = adata[train_idx, :][:, hvg_mask].X
    X = np.clip(X.toarray() if sparse.issparse(X) else np.asarray(X), 0, None)
    print(f"Training matrix {X.shape[0]:,} x {X.shape[1]:,}; fitting NMF ...", flush=True)
    t_fit = time.time()
    model = nmf01.make_nmf(nmf01.DEFAULT_K)
    model.fit_transform(X)
    print(f"  fit in {time.time() - t_fit:.0f}s, {model.n_iter_} iterations", flush=True)
    names = [f"Factor_{i + 1}" for i in range(nmf01.DEFAULT_K)]
    H = pd.DataFrame(model.components_, index=names, columns=hvg_genes)

    common = [g for g in dep_H.columns if g in H.columns]
    Hn, Hd = H.loc[names, common], dep_H.loc[names, common]
    rows = []
    for f in names:
        rows.append({"factor": f, "pearson_r": _r(Hn.loc[f], Hd.loc[f]),
                     "max_abs_diff": float((Hn.loc[f] - Hd.loc[f]).abs().max()),
                     "deposit_max": float(Hd.loc[f].max())})
    fac = pd.DataFrame(rows)
    cross = np.corrcoef(Hn.values, Hd.values)[:len(names), len(names):]
    fac["best_match_in_deposit"] = [names[j] for j in cross.argmax(axis=1)]
    rel_err = float(np.linalg.norm(Hn.values - Hd.values) / np.linalg.norm(Hd.values))

    W = model.transform(X)
    W = W / (W.sum(axis=1, keepdims=True) + 1e-10)
    dep_U = pd.read_csv(os.path.join(DEPOSIT, "11d_nmf_usage.csv"), index_col=0)
    Ud = dep_U.reindex(adata.obs_names[train_idx])[names]
    assert Ud.notna().all().all(), "training barcodes missing from the deposited usage"
    fac["usage_pearson_r"] = [_r(W[:, i], Ud.iloc[:, i]) for i in range(len(names))]
    fac["usage_max_abs_diff"] = [float(np.abs(W[:, i] - Ud.iloc[:, i].values).max()) for i in range(len(names))]
    fac.to_csv(os.path.join(OUT, "loadings_by_factor.csv"), index=False)

    summary = {
        "n_cells": int(adata.n_obs), "n_train": int(len(train_idx)),
        "n_features_9_8": len(feat_98), "n_features_7_7": len(feat_77), "n_features_deposit": len(dep_set),
        "features_equal_9_8": feat_98 == dep_set, "features_equal_7_7": feat_77 == dep_set,
        "forced_in_9_8": sorted(feat_98 - base), "forced_in_7_7": sorted(feat_77 - base),
        "loadings_min_r": float(fac.pearson_r.min()), "loadings_rel_frobenius_err": rel_err,
        "factor_order_preserved": bool((fac.best_match_in_deposit == fac.factor).all()),
        "usage_min_r_training_cells": float(fac.usage_pearson_r.min()),
        "usage_max_abs_diff_training_cells": float(fac.usage_max_abs_diff.max()),
        "nmf_n_iter": int(model.n_iter_), "runtime_s": round(time.time() - t0),
    }
    with open(os.path.join(OUT, "verify_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(fac.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
