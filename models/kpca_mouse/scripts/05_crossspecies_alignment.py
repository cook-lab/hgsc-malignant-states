#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/05_crossspecies_alignment.py, the script that ran.
# Kept: align_gene_space and perm_null_pvalues (permutation null, BH FDR), which step 32 imports. Removed: this
# file's own step-05 analysis of the step-04 programs and its figures. The untrimmed script, which the docstring
# below describes, is in the lab project folder.

"""
Step 05 — Cross-species alignment: does mouse recover the human SecB axis?
===========================================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

This is the falsifiable core of the project. Step 03 TRANSFERRED the human
signature onto mouse cells, which produces a score gradient whether or not the
axis exists. Here we ask the opposite question: does a de novo, unbiased mouse
NMF (step 04) RECOVER a program matching human Factor_2 (SecB)?

Tests, in order of how much they can be trusted:

  1. Loading similarity — every mouse program x every human 11d factor, by
     Spearman and cosine over shared orthologs, with a GENE-LABEL PERMUTATION
     NULL (10k) and BH correction. A high correlation with no permutation
     support is not evidence.
  2. Cell-cycle exclusion — the whole battery re-run with cell-cycle orthologs
     removed. KPCA.B is a clonal line and a match that evaporates here is a
     proliferation artefact, not a SecB axis.
  3. Confound audit — mouse programs flagged ambient/cycling/stress in step 04
     are reported alongside, so a match driven by macrophage ambient signal is
     visible rather than silently believed.
  4. Transfer-vs-recovery concordance — does the best-matching mouse program
     mark the SAME cells as the step-03 transferred SecB score?
  5. Secondary: the 11g consensus metaprograms. Reported but weighted lightly —
     MP3 (SecB) recurs in only 12 human patients across 4 studies.

Input:   output/04_mouse_nmf/04_program_loadings.csv, 04_program_annotation.csv
         output/03_transfer/03_percell_scores.csv
         human 11d loadings, 11g MP gene weights
Output:  output/05_alignment/

Usage:
    python scripts/05_crossspecies_alignment.py
    python scripts/05_crossspecies_alignment.py --n-perm 2000
"""

import os
import sys
from importlib import import_module

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")


def align_gene_space(mouse_loadings, human_loadings, ortho):
    """Restrict both loading matrices to shared 1:1-mapped orthologs."""
    mapped, rep = S.map_genes(human_loadings.columns.tolist(), ortho,
                              restrict_to=set(mouse_loadings.columns))
    hum = [g for g in human_loadings.columns if g in mapped]
    mou = [mapped[g] for g in hum]
    # a mouse gene hit by two human genes would be double-counted
    seen, hkeep, mkeep = set(), [], []
    for h, m in zip(hum, mou):
        if m not in seen:
            seen.add(m); hkeep.append(h); mkeep.append(m)
    Hm = human_loadings[hkeep]
    Mm = mouse_loadings[mkeep]
    print(f"    shared orthologs: {len(mkeep):,} "
          f"(human HVG {human_loadings.shape[1]:,}, "
          f"mouse HVG {mouse_loadings.shape[1]:,})")
    return Mm, Hm, mkeep, hkeep


def perm_null_pvalues(Mm, Hm, n_perm, seed=S.SEED):
    """Spearman + cosine for every pair, with a gene-label permutation null.

    Spearman == Pearson on ranks. Under permutation of one vector's gene labels
    only the dot product varies (the norms are permutation-invariant), so the
    null is computed by shuffling ranks and re-taking dot products — exact and
    fast enough for 10k x 80 pairs.
    """
    rng = np.random.default_rng(seed)
    prog, fac = Mm.index.tolist(), Hm.index.tolist()

    # rank-transform once (handles the many tied zero loadings correctly)
    Rm = np.vstack([rankdata(Mm.loc[p].to_numpy()) for p in prog])
    Rh = np.vstack([rankdata(Hm.loc[f].to_numpy()) for f in fac])
    Rm_c = Rm - Rm.mean(axis=1, keepdims=True)
    Rh_c = Rh - Rh.mean(axis=1, keepdims=True)
    Rm_n = np.linalg.norm(Rm_c, axis=1)
    Rh_n = np.linalg.norm(Rh_c, axis=1)
    obs_rho = (Rm_c @ Rh_c.T) / np.outer(Rm_n, Rh_n)

    # cosine on the raw non-negative loadings
    Mv = Mm.to_numpy(); Hv = Hm.to_numpy()
    obs_cos = ((Mv / (np.linalg.norm(Mv, axis=1, keepdims=True) + 1e-12)) @
               (Hv / (np.linalg.norm(Hv, axis=1, keepdims=True) + 1e-12)).T)

    n = Rm.shape[1]
    ge = np.zeros_like(obs_rho, dtype=float)
    for _ in range(n_perm):
        idx = rng.permutation(n)
        null = (Rm_c[:, idx] @ Rh_c.T) / np.outer(Rm_n, Rh_n)
        ge += (np.abs(null) >= np.abs(obs_rho))
    pval = (ge + 1.0) / (n_perm + 1.0)          # two-sided, add-one

    rows = []
    for i, p in enumerate(prog):
        for j, f in enumerate(fac):
            rows.append({"mouse_program": p, "human_factor": f,
                         "spearman": obs_rho[i, j], "cosine": obs_cos[i, j],
                         "p_perm": pval[i, j]})
    df = pd.DataFrame(rows)
    df["p_adj"] = multipletests(df["p_perm"], method="fdr_bh")[1]
    return df
