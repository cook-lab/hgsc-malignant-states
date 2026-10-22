#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/scripts/01_load_qc_kpca.py, the script that ran. Kept:
# kpca_all.h5ad (read by 15b) and kpca_epithelial.h5ad (the gene space for 02). Removed: the QC, annotation and UMAP
# figures, the ambient report, the Cell Ranger cross-check and two summary tables. The untrimmed script, which the
# docstring below describes, is in the lab project folder.

"""
Step 01 — Load KPCA-B, QC, de novo compartment annotation
==========================================================
KPCA.B mouse SecA/SecB subanalysis  ·  Cook Lab  ·  2026

Counts policy (important):
    CellBender ran on the KPCA-B per-sample RAW matrix and OVER-CALLED cells:
    45,585 barcodes vs CellRanger's 28,973, against its own prior of 20,336
    "probable cell barcodes", convergence_indicator 0.505. Its barcode set is a
    strict superset of CellRanger's. So we take CellBender COUNTS but
    CellRanger CELL CALLS — intersect to the 28,973 and assert it.

    Ambient removal was mild (1.6% of counts), so the CellRanger matrix is a
    cheap sensitivity check (--counts cellranger).

Annotation policy:
    cell_types.csv is a CellRanger automated annotation — treated as a PRIOR to
    cross-check, never as truth. We cluster de novo, label by canonical mouse
    markers, and emit the confusion matrix against it.

Input:   /Volumes/CookLab/Sarah/2025_kpca/KPCA-B/
Output:  output/01_kpca_qc/
             kpca_all.h5ad, kpca_epithelial.h5ad
             01_confusion_vs_cellranger.csv
             01_cluster_marker_scores.csv, 01_compartment_summary.csv
             figures

Usage:
    python scripts/01_load_qc_kpca.py
    python scripts/01_load_qc_kpca.py --counts cellranger
    python scripts/01_load_qc_kpca.py --resolution 1.5
"""

import argparse
import os
import sys
from importlib import import_module

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csc_matrix, csr_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")

OUT_DIR = os.path.join(S.OUTPUT, "01_kpca_qc")


# ============================================================================
# LOADING
# ============================================================================

def read_10x_h5_raw(path):
    """Read a CellRanger/CellBender v3-layout h5 without scanpy's name mangling.

    Returns AnnData (cells x genes) with unique var_names and `gene_ids`.
    """
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        data, indices, indptr = m["data"][:], m["indices"][:], m["indptr"][:]
        shape = m["shape"][:]                        # (genes, cells)
        barcodes = np.array([b.decode() for b in m["barcodes"][:]])
        gene_ids = np.array([b.decode() for b in m["features/id"][:]])
        gene_nm = np.array([b.decode() for b in m["features/name"][:]])

    # 10x stores CSC with shape (genes, cells); transpose to cells x genes
    X = csc_matrix((data, indices, indptr), shape=(int(shape[0]), int(shape[1])))
    A = ad.AnnData(X=csr_matrix(X.T.astype(np.float32)))
    A.obs_names = barcodes
    A.var_names = gene_nm
    A.var["gene_ids"] = gene_ids
    A.var_names_make_unique()
    return A


def load_counts(which):
    """CellBender counts on CellRanger cell calls (default), or CellRanger."""
    print(f"[1] Loading counts ({which})")

    cr_bc = read_10x_h5_raw(S.KPCA_CR).obs_names.to_numpy()
    print(f"    CellRanger called cells: {len(cr_bc):,}")
    assert len(cr_bc) == S.N_CELLS_CELLRANGER, (
        f"CellRanger cell count changed: {len(cr_bc)} != {S.N_CELLS_CELLRANGER}")

    if which == "cellranger":
        A = read_10x_h5_raw(S.KPCA_CR)
        A.uns["counts_source"] = "cellranger_filtered"
        print(f"    -> {A.shape[0]:,} cells x {A.shape[1]:,} genes")
        return A

    A = read_10x_h5_raw(S.KPCA_CB)
    print(f"    CellBender barcodes: {A.shape[0]:,} "
          f"(over-calls; its own prior was 20,336)")

    # CellBender's set must be a superset of CellRanger's, else our assumption
    # about which matrix it ran on is wrong and the intersection is unsafe.
    missing = set(cr_bc) - set(A.obs_names)
    assert not missing, (
        f"{len(missing)} CellRanger barcodes absent from CellBender output — "
        "CellBender may not have been run on the KPCA-B per-sample matrix. STOP.")

    A = A[cr_bc].copy()
    assert A.shape[0] == S.N_CELLS_CELLRANGER, "intersection lost cells"
    A.uns["counts_source"] = "cellbender_on_cellranger_calls"
    print(f"    -> intersected to {A.shape[0]:,} cells x {A.shape[1]:,} genes")
    return A


def attach_cellranger_labels(A):
    """Attach the automated annotation as a prior (never as truth)."""
    lab = pd.read_csv(S.KPCA_ANNOT).set_index("barcode")
    for col in ("coarse_cell_type", "fine_cell_type"):
        A.obs[f"cr_{col}"] = lab[col].reindex(A.obs_names).values
    n = A.obs["cr_coarse_cell_type"].notna().sum()
    print(f"    CellRanger labels attached for {n:,}/{A.shape[0]:,} cells")
    return A


# ============================================================================
# QC
# ============================================================================

def run_qc(A, out_dir):
    """QC metrics + data-driven thresholds.

    Flex is probe-based: there is no intronic/mito capture the way there is in
    3' GEX, so mito% is near-zero by construction and is NOT a useful filter
    here. Thresholds come from the observed distribution (MAD-based), not from
    3'-GEX defaults.
    """
    print("[2] QC")
    A.var["mt"] = A.var_names.str.startswith("mt-")
    A.var["ribo"] = A.var_names.str.match(r"^Rp[sl]")
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt", "ribo"], inplace=True,
                               percent_top=None, log1p=False)

    print(f"    counts   median={A.obs['total_counts'].median():,.0f}  "
          f"range=[{A.obs['total_counts'].min():.0f}, {A.obs['total_counts'].max():.0f}]")
    print(f"    genes    median={A.obs['n_genes_by_counts'].median():,.0f}")
    print(f"    pct_mt   median={A.obs['pct_counts_mt'].median():.2f}%  "
          f"(Flex probe assay — expected low, not used as a filter)")
    print(f"    pct_ribo median={A.obs['pct_counts_ribo'].median():.2f}%")

    def mad_bounds(x, n=5):
        med = np.median(x)
        mad = np.median(np.abs(x - med)) or 1e-9
        return med - n * 1.4826 * mad, med + n * 1.4826 * mad

    lg = np.log1p(A.obs["n_genes_by_counts"].values)
    lo_g, hi_g = mad_bounds(lg)
    min_genes = max(200, int(np.expm1(lo_g)))
    max_genes = int(np.expm1(hi_g))
    print(f"    gene bounds (5x MAD on log1p): [{min_genes:,}, {max_genes:,}]")

    A.obs["qc_pass"] = ((A.obs["n_genes_by_counts"] >= min_genes) &
                        (A.obs["n_genes_by_counts"] <= max_genes))
    print(f"    qc_pass: {A.obs['qc_pass'].sum():,}/{A.shape[0]:,} "
          f"({A.obs['qc_pass'].mean():.1%})")

    A.uns["qc_thresholds"] = {"min_genes": min_genes, "max_genes": max_genes}
    return A


def run_doublets(A):
    """Scrublet scores only — deliberately NOT used as a per-cell hard filter.

    On this Flex library Scrublet's automatic threshold calls 41.6% of cells
    doublets, which is not credible (CellRanger reports
    fraction_reads_high_occupancy_gems = 0.0). Probe-based assays break
    Scrublet's simulation assumptions, and hard-filtering at its threshold
    discarded 12k cells and more than half the epithelium.

    So: keep the continuous score, and remove doublets at CLUSTER level after
    clustering (see drop_doublet_clusters) — the standard robust approach.
    """
    print("[3] Doublet scoring (Scrublet — scores kept, not hard-filtered)")
    try:
        sc.pp.scrublet(A, random_state=S.SEED)
        q = A.obs["doublet_score"].quantile([.5, .9, .95, .99]).round(3)
        print(f"    doublet_score quantiles: "
              f"p50={q[.5]} p90={q[.9]} p95={q[.95]} p99={q[.99]}")
        print(f"    Scrublet's own call would flag "
              f"{A.obs['predicted_doublet'].mean():.1%} — not applied per-cell")
    except Exception as e:                                    # noqa: BLE001
        print(f"    Scrublet failed ({e}); scores set to 0")
        A.obs["doublet_score"] = 0.0
        A.obs["predicted_doublet"] = False
    return A


# ============================================================================
# CLUSTER + ANNOTATE
# ============================================================================

def cluster(A, resolution):
    print(f"[4] Normalize / HVG / PCA / Leiden (res={resolution})")
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    A.raw = A
    sc.pp.highly_variable_genes(A, n_top_genes=S.N_HVG)
    print(f"    HVGs: {int(A.var['highly_variable'].sum()):,}")
    sc.pp.pca(A, n_comps=50, svd_solver="arpack", random_state=S.SEED)
    sc.pp.neighbors(A, n_neighbors=15, n_pcs=40, random_state=S.SEED)
    sc.tl.umap(A, random_state=S.SEED)
    sc.tl.leiden(A, resolution=resolution, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=S.SEED)
    print(f"    clusters: {A.obs['leiden'].nunique()}")
    return A


def compartment_scores(A):
    """Per-cell scaled marker-set scores — ambient-robust.

    This library has SEVERE ambient contamination: Col1a1 is detected in 60% of
    ALL cells (29% of T cells, 65% of epithelial cells) because the tumour is
    collagen-rich and CellBender removed only 1.6% of counts
    (convergence_indicator 0.505). Raw or score_genes values therefore give
    every cell a positive "fibroblast" score, and an argmax over marker sets
    calls almost everything fibroblast.

    Fix: z-score each marker gene ACROSS CELLS first, so the ambient floor
    becomes the zero point and only above-average expression scores positively.
    Averaging z-scores also makes different marker sets directly comparable,
    which raw expression levels are not.
    """
    genes = sorted({g for v in S.MOUSE_COMPARTMENT_MARKERS.values() for g in v
                    if g in A.var_names})
    M = A[:, genes].X
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
    M = (M - M.mean(axis=0)) / (M.std(axis=0) + 1e-9)
    Z = pd.DataFrame(M, index=A.obs_names, columns=genes)

    out = {}
    for comp, mks in S.MOUSE_COMPARTMENT_MARKERS.items():
        present = [g for g in mks if g in Z.columns]
        if not present:
            print(f"    WARNING: no markers found for {comp}")
            continue
        out[comp] = Z[present].mean(axis=1)
        if len(present) < len(mks):
            print(f"    {comp}: {len(present)}/{len(mks)} markers "
                  f"(missing {sorted(set(mks) - set(present))})")
    return pd.DataFrame(out)


def annotate_compartments(A, out_dir):
    """Assign compartments per cluster from scaled marker scores."""
    print("[5] De novo compartment annotation (scaled, ambient-robust)")
    sco = compartment_scores(A)
    for c in sco.columns:
        A.obs[f"score_{c}"] = sco[c].values

    lineage = [c for c in sco.columns if c != "Immune_pan"]
    per_cluster = sco.groupby(A.obs["leiden"].values, observed=True)[lineage].mean()
    per_cluster.index.name = "leiden"

    assign = per_cluster.idxmax(axis=1)
    srt = np.sort(per_cluster.values, axis=1)
    best = pd.Series(srt[:, -1], index=per_cluster.index)
    margin = pd.Series(srt[:, -1] - srt[:, -2], index=per_cluster.index)

    # A cluster must be positively enriched for its winner, and clearly ahead
    # of the runner-up, or it is not called.
    ambiguous = (best < 0.10) | (margin < 0.05)
    assign[ambiguous] = "Other"

    A.obs["compartment"] = A.obs["leiden"].map(assign).astype(str)
    A.obs["compartment_margin"] = A.obs["leiden"].map(margin).astype(float)
    A.obs["compartment_best"] = A.obs["leiden"].map(best).astype(float)

    tbl = per_cluster.copy()
    tbl["assigned"] = assign
    tbl["best_score"] = best
    tbl["margin"] = margin
    tbl["n_cells"] = A.obs["leiden"].value_counts().reindex(tbl.index)
    tbl.to_csv(os.path.join(out_dir, "01_cluster_marker_scores.csv"))
    print("\n" + tbl.round(3).to_string())

    print("\n    compartment sizes:")
    print(A.obs["compartment"].value_counts().to_string())
    return A


def flag_mixed_clusters(A, out_dir, thresh=0.25):
    """Doublet/mixed detection by MUTUAL EXCLUSIVITY, not Scrublet.

    Scrublet calls 41.6% of this library doublets, which is not credible
    (CellRanger: fraction_reads_high_occupancy_gems = 0.0). Ambient
    contamination makes real singlets look like mixtures, which is exactly what
    Scrublet is sensitive to. Removing its clusters deleted the entire
    endothelial and mast compartments.

    Instead: a cluster is "mixed" if TWO different lineage compartments are
    both strongly enriched in it. That is interpretable and ambient-resistant,
    since ambient lifts all cells equally and is removed by the z-scoring.
    """
    print(f"[5a] Mixed-cluster check (two lineages both > {thresh})")
    lineage = [c for c in S.MOUSE_COMPARTMENT_MARKERS if c != "Immune_pan"]
    per_cluster = A.obs.groupby("leiden", observed=True)[
        [f"score_{c}" for c in lineage]].mean()
    per_cluster.columns = lineage
    srt = np.sort(per_cluster.values, axis=1)
    second = pd.Series(srt[:, -2], index=per_cluster.index)
    mixed = second[second > thresh].index.tolist()

    A.obs["cluster_mixed"] = A.obs["leiden"].isin(mixed)
    rep = per_cluster.copy()
    rep["second_best"] = second
    rep["mixed"] = rep.index.isin(mixed)
    rep["n_cells"] = A.obs["leiden"].value_counts().reindex(rep.index)
    rep["scrublet_frac"] = A.obs.groupby("leiden", observed=True)[
        "predicted_doublet"].mean()
    # Interpretable corroboration: co-expression of mutually exclusive lineage
    # markers. Ambient gives ~3% Epcam+Ptprc+; a doublet cluster gives ~60%.
    if {"Epcam", "Ptprc"} <= set(A.var_names):
        dp = A[:, ["Epcam", "Ptprc"]].X
        dp = dp.toarray() if hasattr(dp, "toarray") else np.asarray(dp)
        rep["frac_Epcam_Ptprc_dp"] = pd.Series(
            (dp[:, 0] > 0) & (dp[:, 1] > 0), index=A.obs_names
        ).groupby(A.obs["leiden"], observed=True).mean()
    rep.to_csv(os.path.join(out_dir, "01_mixed_cluster_check.csv"))

    if not mixed:
        print("    no mixed clusters detected — none removed")
        return A
    n = int(A.obs["leiden"].isin(mixed).sum())
    print(f"    removing mixed clusters {mixed} ({n:,} cells)")
    print(rep.loc[mixed].round(3).to_string())
    return A[~A.obs["cluster_mixed"]].copy()


def purify_epithelium(E, out_dir, ambient_ref=None):
    """Remove residual doublet sub-populations from the epithelial compartment.

    The global mixed-cluster check works at whole-compartment resolution and can
    leave a doublet sub-population inside the epithelium. Concretely: an
    Epcam-high cluster was 59% Epcam+/Ptprc+ double-positive against a ~3%
    ambient baseline in a clean myeloid cluster — epithelial/immune doublets.

    So sub-cluster the epithelium and drop sub-clusters whose double-positive
    rate is far above the ambient baseline. `ambient_ref` is that baseline,
    measured on a non-epithelial compartment; if absent we use 0.05.
    """
    print("[8a] Epithelial purity check")
    base = 0.05 if ambient_ref is None else float(ambient_ref)
    thresh = max(3.0 * base, 0.20)
    print(f"    ambient Epcam+/Ptprc+ baseline = {base:.1%}; "
          f"dropping sub-clusters > {thresh:.1%}")

    sc.pp.highly_variable_genes(E, n_top_genes=S.N_HVG)
    sc.pp.pca(E, n_comps=30, svd_solver="arpack", random_state=S.SEED)
    sc.pp.neighbors(E, n_neighbors=15, n_pcs=20, random_state=S.SEED)
    sc.tl.leiden(E, resolution=1.0, key_added="_purity_cl", flavor="igraph",
                 n_iterations=2, directed=False, random_state=S.SEED)

    pairs = [("Epcam", "Ptprc"), ("Epcam", "Pecam1"), ("Epcam", "Col1a1")]
    rows = {}
    for a, b in pairs:
        if a in E.var_names and b in E.var_names:
            M = E[:, [a, b]].X
            M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
            rows[f"dp_{a}_{b}"] = pd.Series(
                (M[:, 0] > 0) & (M[:, 1] > 0), index=E.obs_names)
    dp = pd.DataFrame(rows)
    rep = dp.groupby(E.obs["_purity_cl"].values, observed=True).mean()
    rep["n_cells"] = E.obs["_purity_cl"].value_counts().reindex(rep.index)
    rep["drop"] = rep.get("dp_Epcam_Ptprc", pd.Series(0, index=rep.index)) > thresh
    rep.to_csv(os.path.join(out_dir, "01_epithelial_purity.csv"))
    print(rep.round(3).to_string())

    bad = rep.index[rep["drop"]].tolist()
    if bad:
        n = int(E.obs["_purity_cl"].isin(bad).sum())
        print(f"    dropping epithelial sub-clusters {bad} ({n:,} cells)")
        E = E[~E.obs["_purity_cl"].isin(bad)].copy()
    else:
        print("    no impure epithelial sub-clusters — none dropped")
    del E.obs["_purity_cl"]
    print(f"    clean epithelium: {E.shape[0]:,} cells")
    return E


# ============================================================================
# MAIN
# ============================================================================

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", choices=["cellbender", "cellranger"],
                    default="cellbender")
    ap.add_argument("--resolution", type=float, default=1.0)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    sc.settings.verbosity = 1
    print("=" * 74)
    print("Step 01 — KPCA-B load, QC, de novo annotation")
    print("=" * 74)

    A = load_counts(args.counts)
    A = attach_cellranger_labels(A)
    A.obs["sample_id"] = "KPCA-B_B2"

    A = run_qc(A, OUT_DIR)
    A = run_doublets(A)

    keep = A.obs["qc_pass"].astype(bool)
    print(f"\n    retained after QC: {int(keep.sum()):,}/{A.shape[0]:,}")
    A = A[keep].copy()
    sc.pp.filter_genes(A, min_cells=3)
    print(f"    genes with >=3 cells: {A.shape[1]:,}")

    A = cluster(A, args.resolution)
    # Single pass only. Iterating (remove -> re-cluster -> remove) was unstable:
    # each re-clustering changed cluster identity, and successive rounds ate the
    # mast and fibroblast compartments while inflating "Other". Residual
    # epithelial doublets are handled precisely by purify_epithelium() instead.
    A = annotate_compartments(A, OUT_DIR)
    A = flag_mixed_clusters(A, OUT_DIR)
    A = cluster(A, args.resolution)
    A = annotate_compartments(A, OUT_DIR)

    all_path = os.path.join(OUT_DIR, "kpca_all.h5ad")
    A.write_h5ad(all_path)
    print(f"\n    wrote {all_path}")

    E = A[A.obs["compartment"] == "Epithelial"].copy()
    print(f"\n[8] Epithelial subset: {E.shape[0]:,} cells")

    # Ambient double-positive baseline from the largest non-epithelial compartment
    base = None
    if {"Epcam", "Ptprc"} <= set(A.var_names):
        non_epi = A.obs["compartment"].isin(["Myeloid", "T/NK", "Fibroblast"])
        if non_epi.sum() > 200:
            M = A[non_epi, ["Epcam", "Ptprc"]].X
            M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
            base = float(((M[:, 0] > 0) & (M[:, 1] > 0)).mean())
    E = purify_epithelium(E, OUT_DIR, ambient_ref=base)
    if E.shape[0] < 1000:
        print("    WARNING: fewer than 1,000 epithelial cells — NMF will be "
              "underpowered. Check annotation before proceeding.")
    # Re-embed locally: downstream NMF and epitype calls must live in an
    # epithelial-only space, matching the human epithelial-local UMAP.
    sc.pp.highly_variable_genes(E, n_top_genes=S.N_HVG)
    sc.pp.pca(E, n_comps=50, svd_solver="arpack", random_state=S.SEED)
    sc.pp.neighbors(E, n_neighbors=15, n_pcs=30, random_state=S.SEED)
    sc.tl.umap(E, random_state=S.SEED)
    sc.tl.leiden(E, resolution=0.6, key_added="epi_leiden", flavor="igraph",
                 n_iterations=2, directed=False, random_state=S.SEED)
    E.obsm["X_umap_local"] = E.obsm["X_umap"].copy()
    epi_path = os.path.join(OUT_DIR, "kpca_epithelial.h5ad")
    E.write_h5ad(epi_path)
    print(f"    epithelial clusters: {E.obs['epi_leiden'].nunique()}")
    print(f"    wrote {epi_path}")
    print("\nStep 01 complete.")


if __name__ == "__main__":
    main()
