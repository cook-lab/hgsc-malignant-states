#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/04_mouse_nmf.py, the script that ran. Kept: the
# consensus-NMF functions that step 31 imports; make_figures keeps only the loop that adds program usage to the
# object step 31 saves. Removed: the plots and this file's own step-04 run (main). The untrimmed script, which the
# docstring below describes, is in the lab project folder.

"""
Step 04 — De novo consensus NMF on KPCA.B epithelium
=====================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Adapts the atlas 11d/11f/11g design to a SINGLE sample. There is no
across-patient consensus available here (n=1), so the consensus is taken across
NMF restarts instead, following the cNMF procedure (Kotliar et al.): for each
k in 4..12, fit 20 replicates on independent 80% cell subsamples, pool the k*20
programs, cluster them into exactly k groups, and score that k by clustering
stability (silhouette) x how evenly replicates populate each cluster.
Replicates differ by cell subsample rather than by random seed because NMF here
converges to the same optimum from any initialisation. k is chosen by that score, and the final
programs are the cluster centroids. A program drawn from few seeds is a fitting
artefact, not a cell state, and is reported as unstable.

Every consensus program is then explicitly flagged for the confounds that
dominate a clonal cell line in an ambient-heavy library:
    cycling  — overlap with mouse S/G2M orthologs
    ambient  — usage correlated with non-epithelial lineage transcripts
               (Col1a1 is in 60% of all cells in this library)
    stress   — immediate-early / heat-shock genes
These flags are what step 05 uses to decide whether a cross-species match to
human Factor_2 is real or a proliferation artefact.

Input:   output/01_kpca_qc/kpca_epithelial.h5ad
Output:  output/04_mouse_nmf/

Usage:
    python scripts/04_mouse_nmf.py
    python scripts/04_mouse_nmf.py --seeds 5 --ks 5,6,7    # quick pass
"""

import os
import sys
import time
from importlib import import_module

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import nnls
from scipy.spatial.distance import squareform
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import cosine_similarity

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")

STRESS_GENES = ["Fos", "Fosb", "Jun", "Junb", "Egr1", "Hspa1a", "Hspa1b",
                "Dnajb1", "Hspb1", "Atf3", "Ier3", "Socs3", "Btg2"]
AMBIENT_PROBES = ["Col1a1", "Col1a2", "Dcn", "Ptprc", "Pecam1", "Lyz2", "Cd3e"]


def get_hvg_matrix(E):
    print("[1] HVG selection")
    sc.pp.highly_variable_genes(E, n_top_genes=S.N_HVG)
    genes = E.var_names[E.var["highly_variable"]].tolist()
    X = E[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    X = np.clip(X, 0, None)          # matches 11d
    print(f"    matrix: {X.shape[0]:,} cells x {X.shape[1]:,} HVGs")
    return X, genes


def run_nmf_grid(X, ks, n_seeds):
    """Fit NMF at every (k, seed). Returns programs grouped by k.

    cNMF procedure (Kotliar et al.): for each k, pool the k*n_seeds programs,
    cluster them into exactly k groups, and score that k by clustering
    stability. This replaces an arbitrary cosine-distance cut — which on a
    single sample collapsed all programs into one cluster, because NMF restarts
    on the same cells are far more similar to each other than programs from
    different patients are in the human cross-patient design.
    """
    print(f"[2] NMF grid: k={ks} x {n_seeds} cell-subsampled replicates")
    per_k = {}
    rng = np.random.default_rng(S.SEED)
    n_cells = X.shape[0]
    n_sub = int(round(0.80 * n_cells))
    t0 = time.time()
    for k in ks:
        progs, errs = [], []
        for s in range(n_seeds):
            # Replicates differ by CELL SUBSAMPLE, not by random_state.
            # NMF on this dataset converges to the same optimum from any
            # initialisation, so seed-only replicates are identical and
            # stability saturates at silhouette = 1.000 for every k, which makes
            # the criterion pick the smallest k regardless of the data.
            # Subsampling 80% of cells asks the question that actually matters:
            # does the program survive a change in cell composition?
            idx = rng.choice(n_cells, size=n_sub, replace=False)
            model = NMF(n_components=k, init="random", random_state=s,
                        max_iter=600, tol=1e-4)
            model.fit(X[idx])
            H = model.components_
            H = H / (np.linalg.norm(H, axis=1, keepdims=True) + 1e-12)
            progs.append(H)
            errs.append(model.reconstruction_err_)
        per_k[k] = {"programs": np.vstack(progs),
                    "seed_of": np.repeat(np.arange(n_seeds), k),
                    "recon_err": float(np.mean(errs))}
        print(f"    k={k:2d}: {per_k[k]['programs'].shape[0]:3d} programs, "
              f"recon_err={per_k[k]['recon_err']:.2f} "
              f"({time.time()-t0:.0f}s elapsed)")
    return per_k


def select_k(per_k, out_dir):
    """Choose k by clustering stability (mean silhouette on cosine distance)."""
    print("[3] Selecting k by consensus stability")
    rows = []
    for k, d in per_k.items():
        P = d["programs"]
        D = np.clip(1.0 - cosine_similarity(P), 0, None)
        np.fill_diagonal(D, 0.0)
        Z = linkage(squareform(D, checks=False), method="average")
        lab = fcluster(Z, t=k, criterion="maxclust")
        # silhouette needs >1 cluster and at least one non-singleton
        sil = (silhouette_score(D, lab, metric="precomputed")
               if len(np.unique(lab)) > 1 else np.nan)
        # how evenly seeds distribute across clusters: a stable solution has
        # each cluster drawing from most seeds
        seed_of = d["seed_of"]
        n_seeds = len(np.unique(seed_of))
        rec = [len(np.unique(seed_of[lab == c])) / n_seeds
               for c in np.unique(lab)]
        rows.append({"k": k, "silhouette": sil,
                     "mean_recurrence": float(np.mean(rec)),
                     "min_recurrence": float(np.min(rec)),
                     "n_clusters": len(np.unique(lab)),
                     "recon_err": d["recon_err"]})
        per_k[k]["labels"] = lab
    sel = pd.DataFrame(rows)
    sel["score"] = sel["silhouette"] * sel["mean_recurrence"]
    # Guard against a degenerate criterion: if stability is saturated (every k
    # looks perfect) the score is uninformative and would pick the smallest k.
    if sel["silhouette"].max() - sel["silhouette"].min() < 0.02:
        print("    WARNING: stability is saturated across k "
              f"(range {sel['silhouette'].min():.4f}-{sel['silhouette'].max():.4f}); "
              "the selection criterion is uninformative here. Inspect "
              "04_k_selection.csv and set --force-k deliberately.")
    print(sel.round(4).to_string(index=False))
    best = int(sel.loc[sel["score"].idxmax(), "k"])
    print(f"    ==> selected k = {best} "
          f"(silhouette x mean recurrence maximised)")
    sel.to_csv(os.path.join(out_dir, "04_k_selection.csv"), index=False)
    return best, sel


def consensus_programs(per_k, k, genes):
    """Cluster centroids at the selected k."""
    print(f"[3a] Building consensus programs at k={k}")
    P, lab = per_k[k]["programs"], per_k[k]["labels"]
    seed_of = per_k[k]["seed_of"]
    n_seeds = len(np.unique(seed_of))

    rows, cents = [], []
    for c in np.unique(lab):
        mem = P[lab == c]
        cen = mem.mean(axis=0)
        cen = cen / (np.linalg.norm(cen) + 1e-12)
        cents.append(cen)
        within = cosine_similarity(mem).mean() if len(mem) > 1 else 1.0
        rows.append({"cluster": int(c), "n_programs": int(len(mem)),
                     "n_seeds": int(len(np.unique(seed_of[lab == c]))),
                     "recurrence": len(np.unique(seed_of[lab == c])) / n_seeds,
                     "stability": float(within)})
    summ = pd.DataFrame(rows).sort_values("recurrence", ascending=False
                                          ).reset_index(drop=True)
    order = summ["cluster"].to_numpy()
    W = np.vstack([cents[np.where(np.unique(lab) == c)[0][0]] for c in order])
    summ["program"] = [f"mP{i+1}" for i in range(len(summ))]
    loadings = pd.DataFrame(W, index=summ["program"], columns=genes)
    print(summ[["program", "n_programs", "n_seeds", "recurrence",
                "stability"]].round(3).to_string(index=False))
    weak = summ[summ["recurrence"] < 0.5]
    if len(weak):
        print(f"    NOTE: {len(weak)} program(s) recur in <50% of seeds "
              f"({', '.join(weak['program'])}) — treat as unstable")
    return loadings, summ


def program_usage(X, loadings):
    """NNLS per cell against the consensus loading matrix."""
    print("[4] Per-cell usage (NNLS)")
    A = loadings.to_numpy(dtype=float).T                 # (genes, programs)
    W = np.zeros((X.shape[0], A.shape[1]))
    for i in range(X.shape[0]):
        W[i], _ = nnls(A, X[i])
    Wn = W / (W.sum(axis=1, keepdims=True) + 1e-10)
    return W, Wn


def annotate_programs(E, loadings, usage_n, ortho, out_dir):
    """Top genes + cycling / ambient / stress flags per program."""
    print("[5] Program annotation and confound flags")

    # --- cycling: mouse orthologs of the scanpy/Regev S and G2M sets ---
    cc_h = (sc.queries.__doc__ or "")   # placeholder to keep import surface small
    s_h = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
           "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "SLBP", "CCNE2",
           "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6",
           "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN",
           "POLA1", "CHAF1B", "MRPL36", "E2F8"]
    g2m_h = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
             "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
             "TACC3", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1",
             "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
             "CDCA3", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2",
             "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA",
             "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2",
             "G2E3", "GAS2L3", "CBX5", "CENPA"]
    s_m, _ = S.map_genes(s_h, ortho, restrict_to=set(E.var_names))
    g2m_m, _ = S.map_genes(g2m_h, ortho, restrict_to=set(E.var_names))
    cycle_genes = set(s_m.values()) | set(g2m_m.values())
    print(f"    cell-cycle orthologs available: {len(cycle_genes)}")
    sc.tl.score_genes_cell_cycle(E, s_genes=list(s_m.values()),
                                 g2m_genes=list(g2m_m.values()))
    print("    phase distribution: " +
          ", ".join(f"{k}={v}" for k, v in E.obs["phase"].value_counts().items()))

    # --- ambient: usage correlated with non-epithelial lineage transcripts ---
    probes = [g for g in AMBIENT_PROBES if g in E.var_names]
    Xa = E[:, probes].X
    Xa = Xa.toarray() if hasattr(Xa, "toarray") else np.asarray(Xa)
    Xa = (Xa - Xa.mean(axis=0)) / (Xa.std(axis=0) + 1e-9)
    ambient_score = pd.Series(Xa.mean(axis=1), index=E.obs_names)
    E.obs["ambient_score"] = ambient_score.values

    stress_m = [g for g in STRESS_GENES if g in E.var_names]

    rows = []
    for p in loadings.index:
        top = loadings.loc[p].sort_values(ascending=False)
        top50 = set(top.head(50).index)
        u = pd.Series(usage_n[p].values, index=E.obs_names)
        rows.append({
            "program": p,
            "top20_genes": ", ".join(top.head(20).index),
            "n_cycle_in_top50": len(top50 & cycle_genes),
            "n_stress_in_top50": len(top50 & set(stress_m)),
            "rho_ambient": float(u.corr(ambient_score, method="spearman")),
            "rho_total_counts": float(u.corr(
                pd.Series(E.obs["total_counts"].values, index=E.obs_names),
                method="spearman")),
            "frac_G2M_top_decile": float(
                (E.obs["phase"].values[u >= u.quantile(0.9)] == "G2M").mean()),
            "mean_usage": float(u.mean()),
        })
    ann = pd.DataFrame(rows)
    ann["flag_cycling"] = ann["n_cycle_in_top50"] >= 5
    ann["flag_stress"] = ann["n_stress_in_top50"] >= 3
    ann["flag_ambient"] = ann["rho_ambient"].abs() >= 0.30
    print("\n" + ann[["program", "n_cycle_in_top50", "n_stress_in_top50",
                      "rho_ambient", "mean_usage", "flag_cycling",
                      "flag_stress", "flag_ambient"]].round(3).to_string(index=False))
    for _, r in ann.iterrows():
        print(f"\n    {r['program']}: {r['top20_genes']}")
    return ann, cycle_genes


def make_figures(E, loadings, usage_n, ann, summ, out_dir):
    for p in loadings.index:
        E.obs[p] = usage_n[p].values
