#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/32_identify_secb_align.py, the script that ran. Kept:
# 32_alignment_stats.csv and 32_verdict.csv (read by the Figure S10 script; the FDR of Figure S10c). Removed: the
# program-identity table and the summary figure. The untrimmed script is in the lab project folder.

"""
Step 32 — Identify the SecB program in the full-tumour NMF and align to human
=============================================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Consumes the fresh unsupervised NMF on the full epithelial subset (step 31) and:

  1. IDENTIFY which de-novo program is SecB, three independent ways:
       (a) SecB marker content — mean loading-rank of the canonical mouse SecB
           genes + conserved core in each program's loadings;
       (b) usage vs the UCell SecB polarization (per-cell Spearman);
       (c) loading correlation with human Factor_2 (SecB) over shared orthologs.
     A program that tops all three is SecB; disagreement is reported honestly.
  2. ALIGN to the human 11d factors — extend loadings to all genes (NNLS, W
     fixed), Spearman/cosine vs every factor with a gene-label permutation null,
     then calibrate the SecB-axis match against a floor (ambient/cycling
     programs) and a ceiling (within-tumour split-half NMF). Same verdict logic
     as step 05, now on the full epithelium.

Input:   output/31_fulltumour_nmf/*   ·   output/26_final_annotation/kpca_annotated.h5ad
Output:  output/32_secb_align/
"""
import os, sys
from importlib import import_module
import numpy as np, pandas as pd, scanpy as sc
from scipy.optimize import nnls
from scipy.stats import spearmanr, rankdata
from sklearn.decomposition import NMF

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")
M05 = import_module("05_crossspecies_alignment")
OUT = os.path.join(S.OUTPUT, "32_secb_align")
D31 = os.path.join(S.OUTPUT, "31_fulltumour_nmf")
SECB_MK = ["Krt7", "Krt19", "Lcn2", "Prss22", "Slpi", "Cldn4", "Cldn3", "Cldn6",
           "Krt80", "Wnt7a", "Sfn", "Krt8", "Krt18"]
SECA_MK = ["Bcam", "Mecom", "Sox17", "Pbx1", "Rcn2", "Lgr5"]


def extend_all_genes(usage_raw, E):
    """Refit loadings over all genes with usage W fixed (NNLS)."""
    W = usage_raw.loc[E.obs_names].to_numpy(float)
    X = E.X; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    X = np.clip(X, 0, None)
    H = np.zeros((W.shape[1], X.shape[1]))
    for j in range(X.shape[1]):
        H[:, j], _ = nnls(W, X[:, j])
    L = pd.DataFrame(H, index=usage_raw.columns, columns=E.var_names)
    return L.div(np.linalg.norm(L.to_numpy(), axis=1, keepdims=True) + 1e-12)


def main():
    os.makedirs(OUT, exist_ok=True)
    print("=" * 74, "\nStep 32 — identify SecB program + align to human\n", "=" * 74)
    ML = pd.read_csv(os.path.join(D31, "31_program_loadings.csv"), index_col=0)
    usage = pd.read_csv(os.path.join(D31, "31_program_usage.csv"), index_col=0)
    usage_raw = pd.read_csv(os.path.join(D31, "31_program_usage_raw.csv"), index_col=0)
    ann = pd.read_csv(os.path.join(D31, "31_program_annotation.csv"))
    E = sc.read_h5ad(os.path.join(D31, "kpca_fulltumour_programs.h5ad"))
    A = sc.read_h5ad(os.path.join(S.OUTPUT, "26_final_annotation", "kpca_annotated.h5ad"))
    pol = A.obs.loc[E.obs_names, "polarization"]
    HL = pd.read_csv(S.H_LOADINGS, index_col=0)
    ortho = S.load_ortholog_map()
    progs = list(ML.index)
    print(f"    {len(progs)} programs: {progs}")

    # ---------- 1. identify SecB ----------
    # (a) marker content: mean percentile of SecB markers in each program's loadings
    rank_pct = ML.rank(axis=1, pct=True)
    secb_mk = [g for g in SECB_MK if g in ML.columns]
    seca_mk = [g for g in SECA_MK if g in ML.columns]
    marker_secb = rank_pct[secb_mk].mean(axis=1)
    marker_seca = rank_pct[seca_mk].mean(axis=1)
    # (b) usage vs UCell polarization
    conc = {p: spearmanr(usage[p], pol).correlation for p in progs}
    # (c) loading vs human Factor_2 over shared orthologs
    ML_all = extend_all_genes(usage_raw, E)
    Mm, Hm, _, _ = M05.align_gene_space(ML_all, HL, ortho)
    rho_f2 = {p: spearmanr(Mm.loc[p], Hm.loc[S.H_SECB_FACTOR]).correlation for p in progs}
    rho_f3 = {p: spearmanr(Mm.loc[p], Hm.loc[S.H_SECA_FACTOR]).correlation for p in progs}

    idf = pd.DataFrame({
        "program": progs,
        "secb_marker_pct": [marker_secb[p] for p in progs],
        "seca_marker_pct": [marker_seca[p] for p in progs],
        "usage_vs_polarization": [conc[p] for p in progs],
        "load_vs_Factor2_SecB": [rho_f2[p] for p in progs],
        "load_vs_Factor3_SecA": [rho_f3[p] for p in progs],
        "top12": [", ".join(ML.loc[p].sort_values(ascending=False).head(12).index) for p in progs],
    }).merge(ann[["program", "flag_cycling", "flag_ambient", "flag_stress", "mean_usage"]], on="program")
    idf["secb_score"] = (idf["usage_vs_polarization"].rank() + idf["secb_marker_pct"].rank()
                         + idf["load_vs_Factor2_SecB"].rank())
    idf = idf.sort_values("secb_score", ascending=False).reset_index(drop=True)
    pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 90)
    print("\nProgram identity (ranked for SecB-ness):")
    print(idf[["program", "secb_marker_pct", "usage_vs_polarization", "load_vs_Factor2_SecB",
               "flag_cycling", "flag_ambient", "top12"]].round(3).to_string(index=False))
    clean = idf[~idf["flag_ambient"] & ~idf["flag_cycling"]]
    secb_prog = (clean.iloc[0]["program"] if len(clean) else idf.iloc[0]["program"])
    seca_prog = idf.sort_values("load_vs_Factor3_SecA", ascending=False).iloc[0]["program"]
    print(f"\n  ==> SecB program = {secb_prog}   SecA-leaning program = {seca_prog}")

    # ---------- 2. full alignment battery + calibration ----------
    Hm = Hm.copy()
    Hm.loc["Polarization_SecB_minus_SecA"] = Hm.loc[S.H_SECB_FACTOR] - Hm.loc[S.H_SECA_FACTOR]
    all_df = M05.perm_null_pvalues(Mm, Hm, S.N_PERM); all_df["comparison"] = "all genes"
    all_df.to_csv(os.path.join(OUT, "32_alignment_stats.csv"), index=False)
    secb_axis = all_df[(all_df.human_factor == "Polarization_SecB_minus_SecA")
                       & (all_df.mouse_program == secb_prog)].iloc[0]

    # floor: ambient/cycling programs vs Factor_2 ; ceiling: within-tumour split-half
    flags = ann.set_index("program")
    junk = [p for p in progs if bool(flags.loc[p, "flag_ambient"]) or bool(flags.loc[p, "flag_cycling"])]
    floor = max([spearmanr(Mm.loc[p], Hm.loc[S.H_SECB_FACTOR]).correlation for p in junk], default=0.0)
    sc.pp.highly_variable_genes(E, n_top_genes=S.N_HVG)
    g = E.var_names[E.var["highly_variable"]].tolist()
    X = E[:, g].X; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X); X = np.clip(X, 0, None)
    rng = np.random.default_rng(S.SEED); idx = rng.permutation(X.shape[0]); h = len(idx)//2
    K = ML.shape[0]; Hs = []
    for part in (idx[:h], idx[h:]):
        mo = NMF(n_components=K, init="random", random_state=0, max_iter=600, tol=1e-4); mo.fit(X[part])
        Hc = mo.components_; Hs.append(Hc/(np.linalg.norm(Hc, axis=1, keepdims=True)+1e-12))
    sim = np.array([[spearmanr(Hs[0][i], Hs[1][j]).correlation for j in range(K)] for i in range(K)])
    ceiling = float(np.median(sim.max(axis=1)))
    rho_secb = rho_f2[secb_prog]
    frac = (rho_secb - floor) / max(ceiling - floor, 1e-9)
    concordance = conc[secb_prog]

    # verdict (same logic as step 05)
    facs = [f for f in Hm.index if f.startswith("Factor_")]
    r = pd.Series({f: spearmanr(Mm.loc[secb_prog], Hm.loc[f]).correlation for f in facs}).sort_values(ascending=False)
    secb_is_best = r.index[0] == S.H_SECB_FACTOR
    gap = r.iloc[0] - r.iloc[1]
    even = 1.0/len(progs); dominance = flags.loc[secb_prog, "mean_usage"]/even
    specific = secb_is_best and gap >= 0.05
    strong = frac >= 0.5
    dominant = dominance >= 1.5
    concordant = concordance > 0.3
    if specific and strong and dominant:
        verdict = "RECOVERED AS DOMINANT"
    elif specific and strong:
        verdict = "RECOVERED, NOT DOMINANT"
    elif concordant:
        verdict = "WEAK / NON-SPECIFIC — gene-level conservation stands; axis not recovered as a dominant program"
    else:
        verdict = "NOT RECOVERED"

    summary = {"secb_program": secb_prog, "seca_program": seca_prog, "k": len(progs),
               "rho_secb_factor2": rho_secb, "rho_secb_axis": secb_axis["spearman"],
               "fdr_axis": secb_axis["p_adj"], "usage_vs_polarization": concordance,
               "floor": floor, "ceiling": ceiling, "frac_floor_to_ceiling": frac,
               "best_human_factor": r.index[0], "gap_to_second": gap,
               "dominance_vs_even": dominance, "verdict": verdict}
    pd.DataFrame([summary]).to_csv(os.path.join(OUT, "32_verdict.csv"), index=False)
    print("\n" + "=" * 74, "\nVERDICT\n", "=" * 74)
    for k, v in summary.items():
        print(f"    {k:22s} {v}")

    print("\nsaved", OUT)


if __name__ == "__main__":
    main()
