#!/usr/bin/env python3
"""
Fig 1g prep — per-tumour polarization PC alignment
==================================================
PURPOSE
    For each atlas tumour sample, take the malignant (CopyKAT-aneuploid) epithelial
    cells and ask which principal component of their transcriptional variation the
    SecA/SecB polarization score aligns with (`pol_top_pc`; Fig 1g plots its
    distribution: PC1 in 54% / top-3 in 82% of 250 tumours).

    Method: organoid SecA/SecB scoring sets with scanpy score_genes (polarization =
    SecB - SecA) and score_genes_cell_cycle (Tirosh lists); seurat-flavour top-2000
    HVGs with the NMF, organoid and cell-cycle genes excluded; z-scored target
    matrix; full-SVD PCA with 30 PCs; `pol_top_pc` = the PC with the largest
    |correlation| with the polarization score; >= 150 aneuploid cells per sample.
    partition_pcs() also computes grouped LMG (Shapley) shares of the PC variance
    (polarization, technical, cell cycle, clone; the pc_* columns), which the paper
    does not use; they are computed in the same function as pol_top_pc.

    Data: counts from the final-atlas epithelial subset h5ad (layers["counts"],
    14,470 genes); the sample's cells, study and patient from load_epithelial();
    CopyKAT calls (copykat_prediction.csv, copykat_subclones.csv,
    copykat_subclone_qc.csv) from 19_cnv/per_sample/.

INPUTS
    - data_root/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
        (layers/counts, raw integer counts; rows == load_epithelial() order, asserted)
    - figures/_prep/atlas_final_tables.py :: load_epithelial() (sample_id, study,
        patient_id for the 575,366 epithelial cells)
    - data_root/2026_final_atlas/output/19_cnv/per_sample/<id>/
        copykat_prediction.csv (aneuploid / diploid / not.defined),
        copykat_subclones.csv, copykat_subclone_qc.csv (k, silhouette, monoclonal)

OUTPUTS
    - output_root/_prep_caches/fig1g_polarization_pc_alignment.csv
        one row per sample (sample, study, patient, CopyKAT clone summary, pc_* shares,
        `pol_top_pc`), read by figures/figure1/04_fig1g_polarization_pc_alignment.py
    - printed: samples, PC1 and top-3 shares

MANUSCRIPT PANEL(S): upstream of Fig 1g.

RUNTIME TIER: moderate (~10-20 min; reads ~1/1 of the epithelial counts layer row-wise).

Usage:  PYTHONPATH=. python figures/_prep/fig1g_polarization_pc_alignment.py
        (--samples 10 100  -> debug subset, written to ..._subset.csv)
"""
from __future__ import annotations

import glob
import itertools
import math
import os
import sys
import time
import traceback
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

warnings.simplefilter("ignore")
sc.settings.verbosity = 0

H5AD_EPI = path("data_root", "2026_final_atlas", "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
PS = path("data_root", "2026_final_atlas", "output", "19_cnv", "per_sample")
OUT_CSV = path("output_root", "_prep_caches", "fig1g_polarization_pc_alignment.csv")
MIN_ANEU = 150

# ---- gene sets ------------------------------------------------------------------
# Organoid SecA/SecB signatures: used for the polarization score of this panel.
ORG_SECA = ["SOX17", "MECOM", "FBXO21", "PROM1", "PEG10", "SFRP1", "PAX8"]
ORG_SECB = ["KRT7", "KRT19", "TACSTD2", "CLDN3", "ELF3", "RAB25", "S100A14"]
# NMF discovery signatures (11d) — exclude too, to be safe against circularity.
NMF_SECA = ["BCAM", "MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17", "WT1", "RCN2"]
NMF_SECB = ["KRT17", "KRT19", "KRT7", "LCN2", "MMP7", "PRSS22", "SLPI", "TACSTD2"]
SIG_EXCLUDE = set(ORG_SECA + ORG_SECB + NMF_SECA + NMF_SECB)

S_GENES = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7",
           "DTL", "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
           "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1",
           "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"]
G2M_GENES = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
             "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
             "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK",
             "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
             "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"]
CC_GENES = set(S_GENES + G2M_GENES)


# ---- IO (final atlas) ------------------------------------------------------------
class EpiCounts:
    """Row-wise reader for layers/counts of the epithelial subset h5ad (CSR, int)."""

    def __init__(self, h5ad):
        from anndata.io import read_elem
        self.f = h5py.File(h5ad, "r")
        g = self.f["layers/counts"]
        assert g.attrs["encoding-type"] == "csr_matrix"
        self.data, self.indices = g["data"], g["indices"]
        self.indptr = g["indptr"][:]
        self.shape = tuple(g.attrs["shape"])
        o = self.f["obs"]
        self.obs_names = pd.Index(np.asarray(read_elem(o[o.attrs["_index"]])).astype(str))
        v = self.f["var"]
        self.var_names = pd.Index(np.asarray(read_elem(v[v.attrs["_index"]])).astype(str))

    def rows(self, pos):
        """CSR matrix of the given (sorted) row positions, read in contiguous runs."""
        pos = np.asarray(pos)
        assert np.all(np.diff(pos) > 0)
        breaks = np.where(np.diff(pos) != 1)[0] + 1
        blocks = []
        for run in np.split(pos, breaks):
            a, b = int(run[0]), int(run[-1]) + 1
            lo, hi = int(self.indptr[a]), int(self.indptr[b])
            blocks.append(sp.csr_matrix((self.data[lo:hi], self.indices[lo:hi], self.indptr[a:b + 1] - lo),
                                        shape=(b - a, self.shape[1])))
        return sp.vstack(blocks, format="csr")


def load_sample(sid_dir, key, counts, epi, pos_of):
    """AnnData (cells x genes, raw counts) of all final-atlas epithelial cells of one
    sample (h5ad row order), with obs: pred, subclone (non-reference epithelial cells
    only)."""
    d = f"{PS}/{sid_dir}"
    bcs = epi.index[epi["sample_key"].values == key]
    pos = np.sort(pos_of.reindex(bcs).to_numpy())
    obs_names = counts.obs_names[pos]
    sub = pd.read_csv(f"{d}/copykat_subclones.csv", dtype={"barcode": str})
    pred = pd.read_csv(f"{d}/copykat_prediction.csv", dtype={"cell.names": str})
    A = sc.AnnData(counts.rows(pos), obs=pd.DataFrame(index=obs_names), var=pd.DataFrame(index=counts.var_names))
    A.obs["subclone"] = sub.set_index("barcode")["subclone"].reindex(A.obs_names).values
    A.obs["pred"] = pred.set_index("cell.names")["copykat.pred"].reindex(A.obs_names).values
    assert A.obs["pred"].notna().all(), f"{key}: epithelial cells without a CopyKAT call"
    return A


def malignant_subset(A, which="aneuploid"):
    """which: 'aneuploid' (CopyKAT-malignant, clone defined), the only value used.
    (A already holds only non-reference epithelial cells.)"""
    if which == "aneuploid":
        A = A[A.obs["pred"] == "aneuploid"].copy()
    return A


# ---- LMG / Shapley core (used by partition_pcs for the pc_* shares) ------------
def _lmg_from_cache(cache, groups):
    """Grouped LMG: average incremental R^2 of each group over all orderings. Returns
    {group: lmg_vec(G)} that sum (over groups) to the full-model R^2 per gene."""
    p = len(groups)
    if p == 0:
        return {}
    G = len(next(iter(cache.values())))
    lmg = {g: np.zeros(G) for g in groups}
    fact = [math.factorial(i) for i in range(p + 1)]
    others = list(groups)
    for g in groups:
        rest = [x for x in others if x != g]
        for r in range(len(rest) + 1):
            w = fact[r] * fact[p - r - 1] / fact[p]
            for combo in itertools.combinations(rest, r):
                S = frozenset(combo); Sg = frozenset(combo + (g,))
                lmg[g] += w * (cache[Sg] - cache[S])
    return lmg


def _std(x):
    x = np.asarray(x, float); s = x.std()
    return (x - x.mean()) / (s if s > 0 else 1.0)


def build_design(A, include_cc=True, include_clone=True):
    n = A.n_obs
    groups = {}
    groups["polarization"] = _std(A.obs["polarization"].values.astype(float)).reshape(n, 1)
    groups["technical"] = _std(A.obs["log_counts"].values.astype(float)).reshape(n, 1)
    if include_cc:
        groups["cellcycle"] = np.column_stack([_std(A.obs["S_score"].values.astype(float)),
                                               _std(A.obs["G2M_score"].values.astype(float))])
    if include_clone:
        clones = A.obs["subclone"].astype(str)
        levels = [c for c in clones.unique() if c.startswith("clone_")]
        if len(levels) >= 2:
            D = pd.get_dummies(clones, drop_first=True).values.astype(float)
            D = D - D.mean(0)
            groups["clone"] = D
    return groups


def score_cells(A):
    """Normalize, log1p, score SecA/SecB (organoid), cell cycle (Tirosh), log library size."""
    A.obs["log_counts"] = np.log1p(np.asarray(A.X.sum(1)).ravel())
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    present = set(A.var_names)
    sc.tl.score_genes(A, [g for g in ORG_SECA if g in present], score_name="SecA_score")
    sc.tl.score_genes(A, [g for g in ORG_SECB if g in present], score_name="SecB_score")
    A.obs["polarization"] = A.obs["SecB_score"] - A.obs["SecA_score"]
    sc.tl.score_genes_cell_cycle(A, s_genes=[g for g in S_GENES if g in present],
                                 g2m_genes=[g for g in G2M_GENES if g in present])
    return A


def target_matrix(A, n_hvg=2000, keep_factor_genes=False):
    Asub = A.copy()
    sc.pp.highly_variable_genes(Asub, n_top_genes=min(n_hvg, Asub.n_vars - 1), flavor="seurat")
    keep = Asub.var_names[Asub.var["highly_variable"].values]
    drop = set() if keep_factor_genes else (SIG_EXCLUDE | CC_GENES)
    keep = [g for g in keep if g not in drop]
    X = Asub[:, keep].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    var_raw = X.var(0)
    Xstd = (X - X.mean(0)) / np.where(X.std(0) == 0, 1, X.std(0))
    return Xstd, var_raw, list(keep)


def _cache_grouped(Xstd, groups):
    n, G = Xstd.shape
    keys = list(groups.keys())
    cache = {}
    for r in range(len(keys) + 1):
        for combo in itertools.combinations(keys, r):
            key = frozenset(combo)
            if not combo:
                cache[key] = np.zeros(G); continue
            Xd = np.column_stack([np.ones(n)] + [groups[g] for g in combo])
            beta, *_ = np.linalg.lstsq(Xd, Xstd, rcond=None)
            resid = Xstd - Xd @ beta
            rss = (resid ** 2).sum(0)
            tss = (Xstd ** 2).sum(0)
            cache[key] = 1.0 - rss / np.where(tss == 0, 1, tss)
    return cache


def partition_pcs(A, n_hvg=2000, npc=30, include_cc=True, keep_factor_genes=False):
    from sklearn.decomposition import PCA
    Xstd, var_raw, genes = target_matrix(A, n_hvg=n_hvg, keep_factor_genes=keep_factor_genes)
    npc = int(min(npc, A.n_obs - 1, Xstd.shape[1] - 1))
    pca = PCA(n_components=npc, svd_solver="full").fit(Xstd)
    scores = pca.transform(Xstd)
    lam = pca.explained_variance_
    total_gene_var = Xstd.var(0).sum()
    groups = build_design(A, include_cc=include_cc, include_clone=True)
    if not groups:
        return None
    Sstd = (scores - scores.mean(0)) / np.where(scores.std(0) == 0, 1, scores.std(0))
    cache = _cache_grouped(Sstd, groups)
    lmg = _lmg_from_cache(cache, list(groups.keys()))
    Tstruct = lam.sum()
    out = {f: float((lmg[f] * lam).sum() / Tstruct) for f in lmg}
    full_r2_pc = cache[frozenset(groups.keys())]
    out["residual_within_pcs"] = float(((1 - full_r2_pc) * lam).sum() / Tstruct)
    pol = _std(A.obs["polarization"].values.astype(float))
    cors = np.array([abs(np.corrcoef(Sstd[:, k], pol)[0, 1]) for k in range(npc)])
    out["_pol_top_pc"] = int(np.argmax(cors) + 1)
    out["_pol_corr_pc1"] = float(abs(np.corrcoef(scores[:, 0], pol)[0, 1]))
    out["_pol_max_corr"] = float(cors.max())
    out["_var_top_pcs_frac"] = float(Tstruct / total_gene_var)
    out["_npc"] = npc
    return out


# ---- driver ----------------------------------------------------------------------
def norm(s):
    """canonical integer-string key ("118.0"->"118", "118"->"118")"""
    try:
        return str(int(float(s)))
    except Exception:
        return str(s)


def main(only=None, out_csv=OUT_CSV):
    t0 = time.time()
    epi = load_epithelial(columns=["sample_id_deposited", "study", "patient_id"]).rename(
        columns={"sample_id_deposited": "sample_id"})  # CopyKAT runs are keyed by the deposited sample id (Nath 2021 fix)
    counts = EpiCounts(H5AD_EPI)
    assert counts.obs_names.equals(epi.index), "epithelial h5ad rows != load_epithelial() order"
    pos_of = pd.Series(np.arange(len(counts.obs_names)), index=counts.obs_names)
    # atlas sample_id is "161.0"-style for Olbrecht; key everything by norm()
    epi["sample_key"] = epi["sample_id"].astype(str).map(norm)
    meta = epi.groupby("sample_key", observed=True)[["study", "patient_id"]].first()

    qc = []
    for f in glob.glob(f"{PS}/*/copykat_subclone_qc.csv"):
        if os.path.basename(os.path.dirname(f)).startswith("._"):
            continue
        try:
            qc.append(pd.read_csv(f))
        except Exception:
            pass
    qc = pd.concat(qc, ignore_index=True); qc["key"] = qc["sample_id"].map(norm)
    qcm = qc.set_index("key")[["k", "silhouette", "monoclonal", "n_aneuploid"]].to_dict("index")

    # Eligibility: >= MIN_ANEU CopyKAT-aneuploid cells among the sample's final-atlas
    # epithelial cells (counted here; disagreements with copykat_subclone_qc n_aneuploid are reported).
    epi_sid = epi["sample_key"]
    sids, n_qc_mismatch = [], 0
    for d in sorted(glob.glob(f"{PS}/*/")):
        sid_dir = os.path.basename(d.rstrip("/"))
        if sid_dir.startswith("._"):
            continue
        if not all(os.path.exists(d + x) for x in ["copykat_subclones.csv", "copykat_prediction.csv"]):
            continue
        key = norm(sid_dir)
        if only and key not in only:
            continue
        pred = pd.read_csv(d + "copykat_prediction.csv", dtype={"cell.names": str})
        in_sample = set(epi.index[epi_sid.values == key])
        n_aneu = int(((pred["copykat.pred"] == "aneuploid") & pred["cell.names"].isin(in_sample)).sum())
        if (qcm.get(key, {}).get("n_aneuploid", 0) >= MIN_ANEU) != (n_aneu >= MIN_ANEU):
            n_qc_mismatch += 1
        if n_aneu >= MIN_ANEU:
            sids.append((sid_dir, key))
    print(f"[fig1g] {len(sids)} samples with >= {MIN_ANEU} aneuploid epithelial cells "
          f"({n_qc_mismatch} eligibility disagreements vs copykat_subclone_qc)", flush=True)

    def row_for(sid_dir, key):
        A0 = load_sample(sid_dir, key, counts, epi, pos_of)
        A = score_cells(malignant_subset(A0, which="aneuploid"))
        nclone = sum(1 for c in A.obs["subclone"].astype(str).unique() if c.startswith("clone_"))
        r = {"sample": key, "study": meta.loc[key, "study"], "patient_id": meta.loc[key, "patient_id"],
             "n_aneuploid": A.n_obs, "n_clones": nclone, "k": qcm.get(key, {}).get("k"),
             "silhouette": qcm.get(key, {}).get("silhouette"), "monoclonal": qcm.get(key, {}).get("monoclonal"),
             "df_pol": 1, "df_cc": 2, "df_tech": 1, "df_clone": max(0, nclone - 1)}
        pc = partition_pcs(A, npc=30, include_cc=True)
        for f in pc:
            r[("pc_" + f) if not f.startswith("_") else f.lstrip("_")] = pc[f]
        return r

    rows = []
    for i, (sid_dir, key) in enumerate(sids):
        try:
            rows.append(row_for(sid_dir, key))
            if (i + 1) % 10 == 0 or i == len(sids) - 1:
                pd.DataFrame(rows).to_csv(out_csv, index=False)
                print(f"  [{i + 1}/{len(sids)}] {key} done ({time.time() - t0:.0f}s elapsed)", flush=True)
        except Exception as e:
            print(f"  [{i + 1}/{len(sids)}] {key} FAILED: {e}", flush=True)
            traceback.print_exc()
    out = pd.DataFrame(rows)
    out.to_csv(out_csv, index=False)
    top = out["pol_top_pc"].dropna().astype(int)
    print(f"[fig1g] DONE {len(rows)}/{len(sids)} samples in {time.time() - t0:.0f}s; "
          f"PC1 {int((top == 1).sum())}/{len(top)} ({100 * (top == 1).mean():.1f}%), "
          f"top-3 {100 * (top <= 3).mean():.1f}% -> {out_csv}", flush=True)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*", help="debug: run only these sample ids (writes a _subset CSV)")
    a = ap.parse_args()
    if a.samples:
        main(only=set(a.samples), out_csv=OUT_CSV.replace(".csv", "_subset.csv"))
    else:
        main()
