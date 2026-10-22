#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/11f_consensus_nmf.py, the script that ran.
# Kept: per-patient NMF (k 5-10, 30 runs per k), metaprogram clustering, cell scoring and the 11d comparison that feed
# 11f_mp_summary.csv, 11f_mp_top_genes.csv and 11f_patient_program_counts.csv (Table S5). Removed: figures A-K, the HTML
# report and the CSVs nothing reads. The untrimmed script is in the lab project folder.
"""
Step 11f — Consensus NMF: Per-Patient Metaprogram Discovery
=============================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Validates SecA/SecB by running NMF independently per patient, then
finding recurrent metaprograms across patients (Gavish et al. 2023
geneNMF-style approach).

Strategy:
  Phase 1: Per-patient NMF (30 runs × k∈{5..10}), consensus clustering
  Phase 2: Pool patient programs → hierarchical clustering → metaprograms
  Phase 3: Score all 575k cells against metaprograms
  Phase 4: Compare with 11d (gene-weight + cell-usage correlation)

Input:   celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
         output/11d_epithelial_nmf/11d_nmf_loadings.csv
         output/11d_epithelial_nmf/11d_nmf_usage.csv
Output:  output/11f_consensus_nmf/

Usage:
    python scripts/11f_consensus_nmf.py
    python scripts/11f_consensus_nmf.py --min-cells 500 --n-runs 50
    python scripts/11f_consensus_nmf.py --n-jobs 14
    python scripts/11f_consensus_nmf.py --figures-only
"""

import argparse
import gc
import os
import time
import warnings

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import NMF
from sklearn.preprocessing import normalize
# joblib removed — sequential per-patient to avoid macOS fork deadlocks

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

BASE     = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD     = os.path.join(BASE, "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
OUT_DIR  = os.path.join(BASE, "output", "11f_consensus_nmf")
DIR_11D  = os.path.join(BASE, "output", "11d_epithelial_nmf")

# ============================================================================
# GENE SIGNATURES
# ============================================================================

SECA_GENES = ["BCAM", "MECOM", "FBXO21", "LGR5", "LPAR3", "PBX1", "SOX17", "WT1", "RCN2"]
SECB_GENES = ["KRT17", "KRT19", "KRT7", "LCN2", "MMP7", "PRSS22", "SLPI", "TACSTD2"]

# ============================================================================
# DEFAULTS
# ============================================================================

DEFAULT_MIN_CELLS     = 200
DEFAULT_K_RANGE       = [5, 6, 7, 8, 9, 10]
DEFAULT_N_RUNS        = 30
DEFAULT_COS_THRESHOLD = 0.3
DEFAULT_N_JOBS        = 2
MAX_ITER              = 500
RNG_SEED              = 42


# ============================================================================
# PHASE 1: PER-PATIENT NMF
# ============================================================================

def _run_single_patient_nmf(patient_id, X_pat, k_range, n_runs, max_iter):
    """
    Run NMF n_runs times at each k for one patient.
    Returns L2-normalized loading vectors (n_total_factors × n_genes).
    """
    n_cells, n_genes = X_pat.shape
    all_loadings = []

    for k in k_range:
        for run_idx in range(n_runs):
            seed = hash((patient_id, k, run_idx)) % (2**31)
            model = NMF(
                n_components=k,
                init="random",
                random_state=seed,
                max_iter=max_iter,
            )
            try:
                model.fit(X_pat)
                H = model.components_  # (k, n_genes)
                # L2-normalize each factor
                H_norm = normalize(H, norm="l2", axis=1)
                all_loadings.append(H_norm)
            except Exception:
                pass  # rare convergence failure — skip

    if not all_loadings:
        return None

    # Stack all loading vectors: (n_total, n_genes)
    all_loadings = np.vstack(all_loadings)
    return all_loadings


def _consensus_cluster_patient(all_loadings, min_cluster_size=3, cos_sim_cut=0.5):
    """
    Hierarchical clustering of a patient's NMF loading vectors.
    Returns consensus program centroids (L2-normalized).
    """
    n_vectors = all_loadings.shape[0]
    if n_vectors < min_cluster_size:
        return None

    # Cosine distance
    dists = pdist(all_loadings, metric="cosine")
    Z = linkage(dists, method="average")

    # Cut at cosine_distance = 1 - cos_sim_cut
    clusters = fcluster(Z, t=1.0 - cos_sim_cut, criterion="distance")

    # Centroid of each cluster with enough members
    programs = []
    for cid in np.unique(clusters):
        members = all_loadings[clusters == cid]
        if len(members) >= min_cluster_size:
            centroid = members.mean(axis=0)
            # L2-normalize centroid
            centroid = centroid / (np.linalg.norm(centroid) + 1e-10)
            programs.append(centroid)

    if not programs:
        return None

    return np.array(programs)


def run_patient_nmf(patient_id, X_pat, k_range, n_runs, max_iter):
    """Full per-patient pipeline: NMF runs → consensus programs."""
    all_loadings = _run_single_patient_nmf(patient_id, X_pat, k_range, n_runs, max_iter)
    if all_loadings is None:
        return patient_id, None
    programs = _consensus_cluster_patient(all_loadings)
    return patient_id, programs


def phase1_per_patient_nmf(adata, hvg_genes, k_range, n_runs, n_jobs,
                           min_cells, max_iter):
    """Phase 1: Run NMF per patient, return consensus programs.

    Resumable: saves per-patient results to OUT_DIR/phase1_cache/ so
    crashed runs can pick up where they left off.
    """
    print("\n  ╔══════════════════════════════════════════════════════╗")
    print("  ║  PHASE 1: Per-Patient NMF  (resumable)              ║")
    print("  ╚══════════════════════════════════════════════════════╝")
    t0 = time.time()

    cache_dir = os.path.join(OUT_DIR, "phase1_cache")
    os.makedirs(cache_dir, exist_ok=True)

    patient_ids = adata.obs["patient_id"].values
    unique_patients = np.unique(patient_ids)

    # Filter patients by min_cells
    patient_counts = pd.Series(patient_ids).value_counts()
    eligible = patient_counts[patient_counts >= min_cells].index.tolist()
    excluded = patient_counts[patient_counts < min_cells].index.tolist()
    print(f"    {len(eligible)}/{len(unique_patients)} patients with ≥{min_cells} cells")
    if excluded:
        print(f"    Excluded ({len(excluded)}): {', '.join(sorted(excluded)[:5])}...")

    # Check which patients already have cached results
    done_pids = set()
    for pid in eligible:
        cache_file = os.path.join(cache_dir, f"{pid}.npy")
        if os.path.exists(cache_file):
            done_pids.add(pid)

    todo_pids = [p for p in eligible if p not in done_pids]
    print(f"    Already cached: {len(done_pids)}, remaining: {len(todo_pids)}")

    # Get HVG indices
    hvg_idx = np.array([adata.var_names.get_loc(g) for g in hvg_genes])

    # Process patients sequentially (avoids macOS fork deadlocks)
    # Each patient's 180 NMF fits use BLAS threads for speed
    for i, pid in enumerate(todo_pids):
        t_pat = time.time()
        mask = patient_ids == pid
        n_pat_cells = mask.sum()
        X_pat = adata.X[mask][:, hvg_idx]
        if sparse.issparse(X_pat):
            X_pat = X_pat.toarray()
        X_pat = np.clip(X_pat, 0, None).astype(np.float32)

        _, programs = run_patient_nmf(pid, X_pat, k_range, n_runs, max_iter)

        cache_file = os.path.join(cache_dir, f"{pid}.npy")
        if programs is not None:
            np.save(cache_file, programs)
            n_progs = len(programs)
        else:
            np.save(cache_file, np.array([]))
            n_progs = 0

        elapsed_pat = time.time() - t_pat
        print(f"    [{i+1}/{len(todo_pids)}] {pid}: {n_pat_cells} cells, "
              f"{n_progs} programs ({elapsed_pat:.0f}s)", flush=True)

        del X_pat
        gc.collect()

    # Reassemble all results from cache
    patient_programs = {}
    n_total_programs = 0
    n_missing = 0
    for pid in eligible:
        cache_file = os.path.join(cache_dir, f"{pid}.npy")
        if not os.path.exists(cache_file):
            n_missing += 1
            print(f"    WARNING: no cache for patient {pid}, skipping")
            continue
        programs = np.load(cache_file, allow_pickle=True)
        if programs.size > 0:
            patient_programs[pid] = programs
            n_total_programs += len(programs)

    elapsed = time.time() - t0
    print(f"\n    Phase 1 complete: {len(patient_programs)} patients, "
          f"{n_total_programs} consensus programs ({elapsed:.0f}s)")

    # Also record study mapping for each patient
    patient_study = adata.obs.groupby("patient_id")["study"].first().to_dict()

    return patient_programs, patient_study, eligible, excluded


# ============================================================================
# PHASE 2: METAPROGRAM DISCOVERY
# ============================================================================

def phase2_metaprograms(patient_programs, patient_study, hvg_genes,
                        cos_threshold, min_programs=5, min_patients=3):
    """Phase 2: Cluster patient programs into metaprograms."""
    print("\n  ╔══════════════════════════════════════════════════════╗")
    print("  ║  PHASE 2: Metaprogram Discovery                    ║")
    print("  ╚══════════════════════════════════════════════════════╝")
    t0 = time.time()

    # Pool all programs
    all_programs = []
    program_labels = []  # (patient_id, program_idx)
    for pid in sorted(patient_programs.keys()):
        progs = patient_programs[pid]
        for i in range(len(progs)):
            all_programs.append(progs[i])
            program_labels.append((pid, i))

    all_programs = np.array(all_programs)
    n_programs = len(all_programs)
    print(f"    Pooled {n_programs} consensus programs from "
          f"{len(patient_programs)} patients")

    # Hierarchical clustering
    print(f"    Computing pairwise cosine distances...", flush=True)
    dists = pdist(all_programs, metric="cosine")
    cos_sim_matrix = 1.0 - squareform(dists)

    print(f"    Hierarchical clustering (average linkage, threshold={cos_threshold})...")
    Z = linkage(dists, method="average")
    clusters = fcluster(Z, t=1.0 - cos_threshold, criterion="distance")

    # Filter: require ≥min_programs from ≥min_patients
    mp_data = {}
    for cid in np.unique(clusters):
        member_idx = np.where(clusters == cid)[0]
        member_patients = set(program_labels[i][0] for i in member_idx)
        if len(member_idx) >= min_programs and len(member_patients) >= min_patients:
            mp_data[cid] = {
                "member_idx": member_idx,
                "member_patients": member_patients,
                "n_programs": len(member_idx),
                "n_patients": len(member_patients),
            }

    # Rename metaprograms by size (largest first)
    sorted_mps = sorted(mp_data.keys(), key=lambda c: mp_data[c]["n_programs"],
                        reverse=True)
    metaprograms = {}
    mp_rename = {}
    for new_id, old_id in enumerate(sorted_mps, 1):
        mp_name = f"MP{new_id}"
        mp_rename[old_id] = mp_name
        info = mp_data[old_id]

        # Compute centroid (mean of L2-normed vectors)
        member_vecs = all_programs[info["member_idx"]]
        centroid = member_vecs.mean(axis=0)
        centroid = centroid / (np.linalg.norm(centroid) + 1e-10)

        # Per-member stability: cosine sim to centroid
        sims = member_vecs @ centroid
        stability = sims.mean()

        # Study representation
        contributing_studies = set()
        for i in info["member_idx"]:
            pid = program_labels[i][0]
            contributing_studies.add(patient_study.get(pid, "unknown"))

        # Top genes
        top_gene_idx = np.argsort(centroid)[::-1][:50]
        top_genes = [hvg_genes[i] for i in top_gene_idx]
        top_weights = centroid[top_gene_idx]

        # Recurrence fraction
        total_patients = len(patient_programs)
        recurrence = info["n_patients"] / total_patients

        metaprograms[mp_name] = {
            "centroid": centroid,
            "member_idx": info["member_idx"],
            "n_programs": info["n_programs"],
            "n_patients": info["n_patients"],
            "n_studies": len(contributing_studies),
            "contributing_studies": contributing_studies,
            "contributing_patients": info["member_patients"],
            "stability": stability,
            "recurrence": recurrence,
            "top_genes": top_genes,
            "top_weights": top_weights,
            "member_sims": sims,
        }

    # Build cluster assignment vector for all programs (using new names)
    mp_assignments = []
    for i in range(n_programs):
        cid = clusters[i]
        if cid in mp_rename:
            mp_assignments.append(mp_rename[cid])
        else:
            mp_assignments.append("unassigned")
    mp_assignments = np.array(mp_assignments)

    elapsed = time.time() - t0
    print(f"\n    Found {len(metaprograms)} metaprograms ({elapsed:.1f}s)")
    for mp_name in sorted(metaprograms.keys()):
        mp = metaprograms[mp_name]
        print(f"      {mp_name}: {mp['n_programs']} programs, "
              f"{mp['n_patients']} patients ({mp['recurrence']:.0%}), "
              f"{mp['n_studies']} studies, stability={mp['stability']:.3f}")
        print(f"        Top genes: {', '.join(mp['top_genes'][:10])}")

    return (metaprograms, all_programs, program_labels, Z, cos_sim_matrix,
            mp_assignments, clusters)


# ============================================================================
# PHASE 3: CELL SCORING
# ============================================================================

def phase3_cell_scoring(adata, hvg_genes, metaprograms):
    """Phase 3: Score all cells against metaprograms."""
    print("\n  ╔══════════════════════════════════════════════════════╗")
    print("  ║  PHASE 3: Cell Scoring                              ║")
    print("  ╚══════════════════════════════════════════════════════╝")
    t0 = time.time()

    mp_names = sorted(metaprograms.keys())
    n_mp = len(mp_names)
    n_genes = len(hvg_genes)

    # Build weight matrix (n_mp × n_genes)
    W = np.zeros((n_mp, n_genes))
    for i, mp_name in enumerate(mp_names):
        W[i] = metaprograms[mp_name]["centroid"]

    # Get HVG expression for all cells
    hvg_idx = np.array([adata.var_names.get_loc(g) for g in hvg_genes])
    print(f"    Extracting HVG expression ({adata.n_obs:,} × {n_genes})...", flush=True)
    X_hvg = adata.X[:, hvg_idx]
    if sparse.issparse(X_hvg):
        X_hvg = X_hvg.toarray()
    X_hvg = np.clip(X_hvg, 0, None).astype(np.float32)

    # Score: usage = X_hvg @ W^T
    print(f"    Computing usage scores...", flush=True)
    usage = X_hvg @ W.T  # (n_cells × n_mp)

    # Row-normalize
    row_sums = usage.sum(axis=1, keepdims=True)
    usage_norm = usage / (row_sums + 1e-10)

    usage_df = pd.DataFrame(usage_norm, index=adata.obs.index, columns=mp_names)
    usage_raw_df = pd.DataFrame(usage, index=adata.obs.index, columns=mp_names)

    elapsed = time.time() - t0
    print(f"    Scored {adata.n_obs:,} cells × {n_mp} MPs ({elapsed:.0f}s)")

    del X_hvg
    gc.collect()

    return usage_df, usage_raw_df


# ============================================================================
# PHASE 4: COMPARISON WITH 11d
# ============================================================================

def phase4_compare_11d(metaprograms, usage_df, hvg_genes):
    """Phase 4: Compare metaprograms with 11d factors."""
    print("\n  ╔══════════════════════════════════════════════════════╗")
    print("  ║  PHASE 4: Comparison with 11d                      ║")
    print("  ╚══════════════════════════════════════════════════════╝")
    t0 = time.time()

    # Load 11d results
    loadings_11d = pd.read_csv(os.path.join(DIR_11D, "11d_nmf_loadings.csv"), index_col=0)
    usage_11d = pd.read_csv(os.path.join(DIR_11D, "11d_nmf_usage.csv"), index_col=0)

    factor_names = loadings_11d.index.tolist()
    mp_names = sorted(metaprograms.keys())

    # --- Gene-weight Pearson correlation ---
    # Align genes: use intersection
    shared_genes = sorted(set(hvg_genes) & set(loadings_11d.columns))
    print(f"    Shared HVGs for correlation: {len(shared_genes)}")

    loading_corr = pd.DataFrame(index=mp_names, columns=factor_names, dtype=float)
    for mp_name in mp_names:
        mp_weights = pd.Series(metaprograms[mp_name]["centroid"], index=hvg_genes)
        for factor in factor_names:
            f_weights = loadings_11d.loc[factor, shared_genes]
            m_weights = mp_weights[shared_genes]
            loading_corr.loc[mp_name, factor] = np.corrcoef(m_weights, f_weights)[0, 1]

    print("\n    Gene-weight Pearson r (MP × 11d factor):")
    print(loading_corr.to_string())

    # --- Cell-usage Spearman correlation (50k subsample) ---
    shared_cells = sorted(set(usage_df.index) & set(usage_11d.index))
    rng = np.random.RandomState(RNG_SEED)
    n_sample = min(50_000, len(shared_cells))
    sample_cells = rng.choice(shared_cells, size=n_sample, replace=False)

    from scipy.stats import spearmanr
    usage_corr = pd.DataFrame(index=mp_names, columns=factor_names, dtype=float)
    mp_vals = usage_df.loc[sample_cells, mp_names].values
    f_vals = usage_11d.loc[sample_cells, factor_names].values

    for i, mp_name in enumerate(mp_names):
        for j, factor in enumerate(factor_names):
            r, _ = spearmanr(mp_vals[:, i], f_vals[:, j])
            usage_corr.loc[mp_name, factor] = r

    # --- SecA/SecB identification ---
    secA_set = set(SECA_GENES)
    secB_set = set(SECB_GENES)

    mp_sec_scores = {}
    for mp_name in mp_names:
        top50 = set(metaprograms[mp_name]["top_genes"][:50])
        secA_overlap = len(top50 & secA_set)
        secB_overlap = len(top50 & secB_set)

        # Mean loading of signature genes
        mp_w = pd.Series(metaprograms[mp_name]["centroid"], index=hvg_genes)
        secA_loading = mp_w[mp_w.index.isin(SECA_GENES)].mean()
        secB_loading = mp_w[mp_w.index.isin(SECB_GENES)].mean()

        mp_sec_scores[mp_name] = {
            "secA_in_top50": secA_overlap,
            "secB_in_top50": secB_overlap,
            "secA_mean_loading": secA_loading,
            "secB_mean_loading": secB_loading,
        }

    # Identify best SecA/SecB MPs
    best_secA_mp = max(mp_names, key=lambda m: mp_sec_scores[m]["secA_mean_loading"])
    best_secB_mp = max(mp_names, key=lambda m: mp_sec_scores[m]["secB_mean_loading"])

    print(f"\n    Best SecA metaprogram: {best_secA_mp} "
          f"(loading={mp_sec_scores[best_secA_mp]['secA_mean_loading']:.4f}, "
          f"{mp_sec_scores[best_secA_mp]['secA_in_top50']} genes in top50)")
    print(f"    Best SecB metaprogram: {best_secB_mp} "
          f"(loading={mp_sec_scores[best_secB_mp]['secB_mean_loading']:.4f}, "
          f"{mp_sec_scores[best_secB_mp]['secB_in_top50']} genes in top50)")

    elapsed = time.time() - t0
    print(f"    Phase 4 complete ({elapsed:.0f}s)")

    return (loading_corr, usage_corr, mp_sec_scores, best_secA_mp, best_secB_mp,
            loadings_11d, usage_11d)


# ============================================================================
# SAVE CSVs
# ============================================================================

def save_csvs(metaprograms, usage_df, loading_corr, usage_corr,
              mp_sec_scores, hvg_genes, patient_programs, patient_study,
              out_dir):
    """Save all CSV outputs."""
    print("\n  [SAVE] Writing CSVs...")

    mp_names = sorted(metaprograms.keys())

    # 2. MP top genes
    rows = []
    for mp_name in mp_names:
        mp = metaprograms[mp_name]
        for rank, (gene, w) in enumerate(zip(mp["top_genes"], mp["top_weights"]), 1):
            rows.append({"metaprogram": mp_name, "rank": rank, "gene": gene, "weight": w})
    top_df = pd.DataFrame(rows)
    top_df.to_csv(os.path.join(out_dir, "11f_mp_top_genes.csv"), index=False)
    print(f"    11f_mp_top_genes.csv ({len(top_df)} rows)")

    # 3. MP summary
    summary_rows = []
    for mp_name in mp_names:
        mp = metaprograms[mp_name]
        summary_rows.append({
            "metaprogram": mp_name,
            "n_programs": mp["n_programs"],
            "n_patients": mp["n_patients"],
            "n_studies": mp["n_studies"],
            "recurrence": mp["recurrence"],
            "stability": mp["stability"],
            "studies": ", ".join(sorted(mp["contributing_studies"])),
            "top10_genes": ", ".join(mp["top_genes"][:10]),
            "secA_in_top50": mp_sec_scores[mp_name]["secA_in_top50"],
            "secB_in_top50": mp_sec_scores[mp_name]["secB_in_top50"],
            "secA_mean_loading": mp_sec_scores[mp_name]["secA_mean_loading"],
            "secB_mean_loading": mp_sec_scores[mp_name]["secB_mean_loading"],
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(out_dir, "11f_mp_summary.csv"), index=False)
    print(f"    11f_mp_summary.csv ({len(summary_df)} MPs)")

    # 7. Per-patient program counts
    pat_rows = []
    for pid in sorted(patient_programs.keys()):
        pat_rows.append({
            "patient_id": pid,
            "study": patient_study.get(pid, "unknown"),
            "n_programs": len(patient_programs[pid]),
        })
    pat_df = pd.DataFrame(pat_rows)
    pat_df.to_csv(os.path.join(out_dir, "11f_patient_program_counts.csv"), index=False)
    print(f"    11f_patient_program_counts.csv ({len(pat_df)} patients)")


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="11f: Consensus NMF — Per-Patient Metaprogram Discovery")
    parser.add_argument("--min-cells", type=int, default=DEFAULT_MIN_CELLS,
                        help=f"Min epithelial cells per patient (default: {DEFAULT_MIN_CELLS})")
    parser.add_argument("--k-range", type=str, default=",".join(map(str, DEFAULT_K_RANGE)),
                        help="Comma-separated k values (default: 5,6,7,8,9,10)")
    parser.add_argument("--n-runs", type=int, default=DEFAULT_N_RUNS,
                        help=f"NMF runs per k per patient (default: {DEFAULT_N_RUNS})")
    parser.add_argument("--cosine-threshold", type=float, default=DEFAULT_COS_THRESHOLD,
                        help=f"Cosine sim threshold for MP clustering (default: {DEFAULT_COS_THRESHOLD})")
    parser.add_argument("--n-jobs", type=int, default=DEFAULT_N_JOBS,
                        help=f"Parallel workers (default: {DEFAULT_N_JOBS})")
    parser.add_argument("--figures-only", action="store_true",
                        help="Load saved Phase 1 results and regenerate Phase 2+")
    args = parser.parse_args()

    k_range = [int(x) for x in args.k_range.split(",")]

    os.makedirs(OUT_DIR, exist_ok=True)

    t0_total = time.time()
    print("=" * 65)
    print("  11f — Consensus NMF: Per-Patient Metaprogram Discovery")
    print("=" * 65)
    print(f"  Parameters:")
    print(f"    min_cells       = {args.min_cells}")
    print(f"    k_range         = {k_range}")
    print(f"    n_runs          = {args.n_runs}")
    print(f"    cosine_threshold = {args.cosine_threshold}")
    print(f"    n_jobs          = {args.n_jobs}")

    # ---- Load shared HVG list from 11d ----
    loadings_11d_path = os.path.join(DIR_11D, "11d_nmf_loadings.csv")
    if not os.path.exists(loadings_11d_path):
        raise FileNotFoundError(f"11d loadings not found: {loadings_11d_path}")

    hvg_genes = pd.read_csv(loadings_11d_path, index_col=0, nrows=0).columns.tolist()
    print(f"\n  Shared HVG space: {len(hvg_genes)} genes (from 11d)")

    # ---- Check for saved Phase 1 results ----
    npz_path = os.path.join(OUT_DIR, "11f_patient_programs.npz")

    adata = None  # track whether we have expression data loaded

    if args.figures_only and os.path.exists(npz_path):
        # Load saved Phase 1
        print(f"\n  [LOAD] Reading saved Phase 1 from {npz_path}...")
        data = np.load(npz_path, allow_pickle=True)
        patient_programs = data["patient_programs"].item()
        patient_study = data["patient_study"].item()
        eligible = data["eligible"].tolist()
        excluded = data["excluded"].tolist()
        print(f"    {len(patient_programs)} patients, "
              f"{sum(len(v) for v in patient_programs.values())} programs")
    else:
        # ---- Load epithelial data ----
        print("\n  [LOAD] Opening epithelial h5ad (backed)...", flush=True)
        t_load = time.time()

        bdata = ad.read_h5ad(H5AD, backed="r")
        n_cells = bdata.n_obs
        genes = bdata.var_names.tolist()

        # Extract metadata
        level2 = bdata.obs["celltype_level2"].values.astype(str)
        patient_ids = bdata.obs["patient_id"].values.astype(str)
        study_vals = bdata.obs["study"].values.astype(str)
        barcodes = bdata.obs.index.values.copy()
        umap_local = np.array(bdata.obsm["X_umap_local"])

        # Extract counts — enforce sparse float32 to limit memory
        print(f"    Extracting counts ({n_cells:,} × {len(genes):,})...", flush=True)
        counts = bdata.layers["counts"][:, :]
        bdata.file.close()
        del bdata
        gc.collect()

        if sparse.issparse(counts):
            counts = counts.astype(np.float32)
        else:
            print("    Converting dense → sparse CSR float32...")
            counts = sparse.csr_matrix(counts, dtype=np.float32)

        # Build in-memory AnnData
        adata = ad.AnnData(
            X=counts,
            obs=pd.DataFrame({
                "celltype_level2": level2,
                "patient_id": patient_ids,
                "study": study_vals,
            }, index=barcodes),
            var=pd.DataFrame(index=genes),
        )
        adata.obsm["X_umap_local"] = umap_local
        del counts
        gc.collect()

        print(f"    Loaded in {time.time()-t_load:.0f}s")

        # ---- Normalize (match 11d) ----
        print("  [PREPROCESS] normalize_total(1e4) → log1p...", flush=True)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Verify HVGs are present
        missing = [g for g in hvg_genes if g not in adata.var_names]
        if missing:
            print(f"    WARNING: {len(missing)} HVGs not in var_names, dropping them")
            hvg_genes = [g for g in hvg_genes if g in adata.var_names]

        # Subset to HVGs — frees ~75% of expression matrix before NMF
        print(f"  [MEMORY] Subsetting to {len(hvg_genes)} HVGs "
              f"(dropping {len(genes) - len(hvg_genes)} non-HVG genes)...")
        adata = adata[:, hvg_genes].copy()
        gc.collect()
        mem_gb = (adata.X.data.nbytes / 1e9 if sparse.issparse(adata.X)
                  else adata.X.nbytes / 1e9)
        print(f"    Expression matrix after subset: {mem_gb:.1f} GB")

        # ---- Phase 1 ----
        patient_programs, patient_study, eligible, excluded = phase1_per_patient_nmf(
            adata, hvg_genes, k_range, args.n_runs, args.n_jobs,
            args.min_cells, MAX_ITER
        )

        # Save Phase 1 results
        print(f"  [SAVE] Saving Phase 1 to {npz_path}...")
        np.savez_compressed(
            npz_path,
            patient_programs=patient_programs,
            patient_study=patient_study,
            eligible=np.array(eligible),
            excluded=np.array(excluded),
        )

    # ---- Phase 2: Metaprogram discovery ----
    (metaprograms, all_programs, program_labels, Z, cos_sim_matrix,
     mp_assignments, clusters) = phase2_metaprograms(
        patient_programs, patient_study, hvg_genes,
        args.cosine_threshold
    )

    if len(metaprograms) == 0:
        print("\n  ERROR: No metaprograms found. Try adjusting --cosine-threshold.")
        return

    # ---- Load full data for Phase 3 + figures (if not already loaded) ----
    if adata is None:
        print("\n  [LOAD] Re-loading epithelial data for scoring + figures...", flush=True)
        t_load = time.time()

        bdata = ad.read_h5ad(H5AD, backed="r")
        genes = bdata.var_names.tolist()
        level2 = bdata.obs["celltype_level2"].values.astype(str)
        patient_ids = bdata.obs["patient_id"].values.astype(str)
        study_vals = bdata.obs["study"].values.astype(str)
        barcodes = bdata.obs.index.values.copy()
        umap_local = np.array(bdata.obsm["X_umap_local"])
        counts = bdata.layers["counts"][:, :]
        bdata.file.close()
        del bdata
        gc.collect()

        if sparse.issparse(counts):
            counts = counts.astype(np.float32)
        else:
            counts = sparse.csr_matrix(counts, dtype=np.float32)

        adata = ad.AnnData(
            X=counts,
            obs=pd.DataFrame({
                "celltype_level2": level2,
                "patient_id": patient_ids,
                "study": study_vals,
            }, index=barcodes),
            var=pd.DataFrame(index=genes),
        )
        adata.obsm["X_umap_local"] = umap_local
        del counts
        gc.collect()

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Subset to HVGs to limit memory
        adata = adata[:, hvg_genes].copy()
        gc.collect()

        print(f"    Loaded + normalized in {time.time()-t_load:.0f}s")

    # ---- Phase 3: Cell scoring ----
    usage_df, usage_raw_df = phase3_cell_scoring(adata, hvg_genes, metaprograms)

    # ---- Phase 4: Compare with 11d ----
    (loading_corr, usage_corr, mp_sec_scores, best_secA_mp, best_secB_mp,
     loadings_11d, usage_11d) = phase4_compare_11d(metaprograms, usage_df, hvg_genes)

    # ---- Save CSVs ----
    save_csvs(metaprograms, usage_df, loading_corr, usage_corr,
              mp_sec_scores, hvg_genes, patient_programs, patient_study, OUT_DIR)

    elapsed_total = time.time() - t0_total

    # ---- Summary ----
    print(f"\n{'=' * 65}")
    print(f"  11f COMPLETE ({elapsed_total:.0f}s)")
    print(f"  Output: {OUT_DIR}")
    print(f"  Metaprograms: {len(metaprograms)}")
    print(f"  Best SecA: {best_secA_mp} (recurrence={metaprograms[best_secA_mp]['recurrence']:.0%})")
    print(f"  Best SecB: {best_secB_mp} (recurrence={metaprograms[best_secB_mp]['recurrence']:.0%})")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
