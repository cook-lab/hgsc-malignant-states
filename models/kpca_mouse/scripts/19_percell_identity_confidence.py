#!/usr/bin/env python3
# Trimmed copy (2026-10-22) of 202608_kpca_secb_mouse/scripts/19_percell_identity_confidence.py, the script that ran
# (2026-09-04). Kept: output/19_identity_confidence/19_percell.csv, whose cnv_score 21a_add_cnv_score.py adds to the
# step-15b object. Removed: the two figures (19_A_percell_identity, 19_B_emt_continuum) and the matplotlib import.
# The untrimmed script, which the docstring below describes, is in the lab project folder.

"""
Step 19 — Per-cell identity confidence: epithelial-tumour vs mesothelial vs
          fibroblast vs pericyte, and the EMT continuum within the tumour
============================================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

CNV (step 18) settled tumour-vs-stroma at cluster level. This resolves what is
left, per cell:
  1. a fast per-cell CNV score for ALL cells (not just the inferCNV subsample),
     so tumour confidence is per-cell — reference = immune, profile correlated
     to the epithelial-anchor consensus, position-ordered and smoothed
  2. per-cell scores for four competing identities using SPECIFIC, minimally
     overlapping markers (mesothelial uses Upk3b/Wt1/Lrrn4 — NOT Msln, which is
     shared with epithelium)
  3. a per-cell call with a confidence = softmax margin, and an entropy, so the
     genuinely ambiguous cells surface rather than being forced
  4. characterisation of the EMT continuum: is the tumour a discrete
     epithelial/mesenchymal split or a gradient, and does CNV hold across it

Input:   output/15_doublet_cleanup/kpca_all_clean.h5ad
         output/18_infercnv/18_cnv_correlation.csv  (cluster CNV anchors)
Output:  output/19_identity_confidence/

Usage:
    python scripts/19_percell_identity_confidence.py
"""

import os
import sys
from importlib import import_module

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")

OUT_DIR = os.path.join(S.OUTPUT, "19_identity_confidence")
CLEAN = os.path.join(S.OUTPUT, "15_doublet_cleanup", "kpca_all_clean.h5ad")
POS = os.path.join(S.OUTPUT, "18_infercnv", "gene_positions.tsv")

# Competing identities with SPECIFIC markers (deliberately non-overlapping;
# Msln excluded from mesothelial because it is shared with epithelium/CAF).
IDENTITY = {
    "Epithelial-tumour": ["Epcam", "Cdh1", "Krt8", "Krt18", "Krt19", "Cldn3",
                          "Cldn4", "Krt7"],
    "Mesothelial":       ["Upk3b", "Wt1", "Lrrn4", "Nkain4", "Gpm6a"],
    "Fibroblast":        ["Pdgfra", "Dcn", "Col1a1", "Lum", "Clec3b", "Pi16",
                          "Col14a1"],
    "Pericyte":          ["Rgs5", "Notch3", "Kcnj8", "Higd1b", "Cspg4"],
}
# clusters at the tumour/stroma boundary worth resolving per cell
FOCUS = ["2", "3", "7", "8", "9", "10", "11", "12", "13", "6", "15"]

# immune reference clusters for CNV baseline (diploid)
IMMUNE_REF = ["4", "5", "16", "17", "0", "1", "14", "18"]
EPI_ANCHOR = ["10", "13"]


def percell_cnv(A):
    """Lightweight per-cell CNV score for all cells (infercnvpy-style).

    reference = immune mean; per-gene residual; position-ordered; rolling mean
    per chromosome; per-cell score = correlation to the epithelial-anchor
    consensus CNV profile. Same principle validated against the full inferCNV
    run in step 18.
    """
    print("[1] Per-cell CNV score (all cells)")
    pos = pd.read_csv(POS, sep="\t", header=None,
                      names=["gene", "chr", "start", "stop"])
    genes = [g for g in pos["gene"] if g in A.var_names]
    pos = pos[pos["gene"].isin(genes)].copy()
    X = A[:, pos["gene"].tolist()].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)   # log-norm
    # residual vs immune reference
    ref = A.obs["leiden"].isin(IMMUNE_REF).values
    resid = X - X[ref].mean(axis=0, keepdims=True)
    resid = np.clip(resid, -3, 3)
    # smooth within each chromosome (rolling mean, window 100)
    sm = np.zeros_like(resid)
    for ch in pos["chr"].unique():
        idx = np.where(pos["chr"].values == ch)[0]
        if len(idx) < 5:
            sm[:, idx] = resid[:, idx]
            continue
        w = min(100, len(idx))
        block = pd.DataFrame(resid[:, idx].T).rolling(
            w, center=True, min_periods=1).mean().to_numpy().T
        sm[:, idx] = block
    # consensus from epithelial anchors, correlation per cell
    anc = A.obs["leiden"].isin(EPI_ANCHOR).values
    consensus = sm[anc].mean(axis=0)
    cons_c = consensus - consensus.mean()
    smc = sm - sm.mean(axis=1, keepdims=True)
    cnv = (smc @ cons_c) / (np.linalg.norm(smc, axis=1) *
                            np.linalg.norm(cons_c) + 1e-12)
    A.obs["cnv_score"] = cnv
    ref_p95 = np.quantile(cnv[ref], 0.95)
    A.obs["is_aneuploid"] = cnv > ref_p95
    print(f"    immune-ref CNV mean={cnv[ref].mean():+.3f}, p95={ref_p95:.3f}")
    print(f"    epithelial-anchor CNV mean={cnv[anc].mean():+.3f}")
    return A, ref_p95


def identity_scores(A):
    print("[2] Per-cell competing-identity scores")
    for name, gl in IDENTITY.items():
        present = [g for g in gl if g in A.var_names]
        sc.tl.score_genes(A, present, score_name=f"id_{name}", random_state=S.SEED)
    cols = [f"id_{k}" for k in IDENTITY]
    Z = A.obs[cols].copy()
    # z-score each identity across cells so they are comparable
    Zz = (Z - Z.mean()) / (Z.std() + 1e-9)
    Zz.columns = list(IDENTITY)
    # softmax over identities -> per-cell probability, call, confidence, entropy
    e = np.exp(Zz.to_numpy() - Zz.to_numpy().max(axis=1, keepdims=True))
    P = e / e.sum(axis=1, keepdims=True)
    Pdf = pd.DataFrame(P, index=A.obs_names, columns=list(IDENTITY))
    call = Pdf.idxmax(axis=1)
    srt = np.sort(P, axis=1)
    margin = srt[:, -1] - srt[:, -2]
    ent = -(P * np.log(P + 1e-12)).sum(axis=1) / np.log(P.shape[1])
    A.obs["id_call"] = call.values
    A.obs["id_confidence"] = margin
    A.obs["id_entropy"] = ent
    for k in IDENTITY:
        A.obs[f"P_{k}"] = Pdf[k].values
    return A, Pdf


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 78)
    print("Step 19 — per-cell identity confidence")
    print("=" * 78)
    A = sc.read_h5ad(CLEAN)
    A.obs["leiden"] = A.obs["leiden"].astype(str)
    order = sorted(A.obs["leiden"].unique(), key=int)

    A, ref_p95 = percell_cnv(A)
    A, Pdf = identity_scores(A)

    # combine CNV + identity into a final per-cell status for the tumour/stroma set
    print("[3] Per-cluster resolution")
    rows = []
    for c in order:
        m = (A.obs["leiden"] == c).values
        sub = A.obs[m]
        vc = sub["id_call"].value_counts(normalize=True)
        rows.append({
            "leiden": c, "n": int(m.sum()),
            "aneuploid%": round(100 * sub["is_aneuploid"].mean()),
            "cnv": round(float(sub["cnv_score"].mean()), 3),
            "top_id": vc.index[0], "top_id%": round(100 * vc.iloc[0]),
            "mean_conf": round(float(sub["id_confidence"].mean()), 2),
            "mean_entropy": round(float(sub["id_entropy"].mean()), 2),
            "ambiguous%": round(100 * (sub["id_confidence"] < 0.15).mean()),
        })
    tab = pd.DataFrame(rows)
    pd.set_option("display.width", 220)
    print(tab.to_string(index=False))
    tab.to_csv(os.path.join(OUT_DIR, "19_cluster_resolution.csv"), index=False)
    A.obs[["leiden", "cnv_score", "is_aneuploid", "id_call", "id_confidence",
           "id_entropy"] + [f"P_{k}" for k in IDENTITY]].to_csv(
        os.path.join(OUT_DIR, "19_percell.csv"))

    n_amb = int((A.obs["id_confidence"] < 0.15).sum())
    print(f"\n    genuinely ambiguous cells (confidence<0.15): {n_amb:,} "
          f"({100*n_amb/A.shape[0]:.1f}%)")
    print(f"    wrote {OUT_DIR}")
    print("\nStep 19 complete.")


if __name__ == "__main__":
    main()
