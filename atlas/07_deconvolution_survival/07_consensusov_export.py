#!/usr/bin/env python3
"""
Pseudobulk export for consensusOV TCGA subtype mapping
======================================================
HGSC malignant-states atlas backend.

Streams the atlas h5ad once (single sequential CSR pass; the object is large on
a slow volume) and builds per-sample pseudobulk for the primary-untreated
ovarian cohort (anatomic_site == adnexa, treatment_status == pre-treatment;
>=500 total + >=100 epithelial cells/sample): a bulk-tumor pseudobulk and an
epithelial-only pseudobulk, plus NMF epitype proportions
(SecA / Intermediate / SecB / Ciliated). Feeds consensusOV scoring (08).

INPUTS:
  - hgsc_atlas_final.h5ad   [config: atlas_final]   (obs + X via h5py)
  - output_root/03_epithelial_nmf/celltype_nmf_mapping.csv

OUTPUTS (output_root/07_deconvolution_survival/consensusov/):
  - pseudobulk_bulk_counts.tsv.gz, pseudobulk_epi_counts.tsv.gz
  - pseudobulk_metadata.csv, gene_list.tsv

MANUSCRIPT PANELS: upstream of Fig 3H (epitype x TCGA subtype dot matrix).

RUNTIME TIER: heavy (single sequential pass over the full atlas X).

SEEDING: deterministic (streaming aggregation; no stochastic step).

Usage:
    python 07_consensusov_export.py
"""

import os
import sys
import gc
import time
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import h5py
import scipy.sparse as sps

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path  # noqa: E402

warnings.filterwarnings("ignore")

# ── Paths ───────────────────────────────────────────────────────────
ATLAS_H5AD = obj("atlas_final")
NMF_MAP    = path("output_root", "03_epithelial_nmf", "celltype_nmf_mapping.csv")
OUT_DIR    = path("output_root", "07_deconvolution_survival", "consensusov")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Cohort filter ───────────────────────────────────────────────────
ANATOMIC_SITE    = "adnexa"
TREATMENT_STATUS = "pre-treatment"
MIN_CELLS_TOTAL  = 500
MIN_CELLS_EPI    = 100

# NMF epitype labels ("Transitioning" -> "Intermediate")
EPITYPE_LABELS = {
    "SecA epithelium":         "SecA",
    "Intermediate epithelium": "Intermediate",
    "SecB epithelium":         "SecB",
}
CILIATED_LABEL = "Ciliated epithelial cell"

NNZ_PER_CHUNK = 200_000_000


def _read_obs_col(f, key):
    g = f["obs"][key]
    if isinstance(g, h5py.Group):
        codes = g["codes"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
        return np.where(codes >= 0, cats[np.maximum(codes, 0)], None)
    arr = g[:]
    if arr.dtype.kind in ("S", "O"):
        arr = np.array([c.decode() if isinstance(c, bytes) else c for c in arr])
    return arr


def _read_var_names(f):
    arr = f["var"]["_index"][:]
    if arr.dtype.kind in ("S", "O"):
        arr = np.array([c.decode() if isinstance(c, bytes) else c for c in arr])
    return arr


def main():
    t_start = time.time()
    print("=" * 70)
    print("  Pseudobulk export for consensusOV TCGA subtype mapping")
    print("=" * 70, flush=True)

    print("\n[1] Loading NMF label mapping...", flush=True)
    nmf = pd.read_csv(NMF_MAP, index_col=0)
    print(f"    {len(nmf):,} barcodes; columns = {list(nmf.columns)}", flush=True)

    print("\n[2] Reading atlas obs via h5py...", flush=True)
    t0 = time.time()
    f = h5py.File(ATLAS_H5AD, "r")
    obs_cols = ["sample_id", "patient_id", "study", "anatomic_site",
                "treatment_status", "celltype_level1", "celltype_level2"]
    obs = pd.DataFrame({c: _read_obs_col(f, c) for c in obs_cols})
    obs.index = _read_obs_col(f, "_index")
    obs.index.name = "barcode"
    var_names = _read_var_names(f)
    n_cells_total = len(obs)
    n_genes = len(var_names)
    print(f"    obs: {len(obs):,} cells x {len(obs.columns)} cols  ({time.time()-t0:.1f}s)")
    print(f"    var: {n_genes:,} genes", flush=True)

    print("\n[3] Applying cohort filter + sample thresholds...", flush=True)
    cohort_mask = ((obs["anatomic_site"] == ANATOMIC_SITE) &
                   (obs["treatment_status"] == TREATMENT_STATUS))
    cohort = obs.loc[cohort_mask].copy()
    print(f"    Cohort cells (adnexa + pre-treatment): {len(cohort):,}", flush=True)
    cohort = cohort.join(nmf[["celltype_nmf"]].rename(columns={"celltype_nmf": "nmf_label"}),
                         how="left")

    by_s = cohort.groupby("sample_id")
    n_total = by_s.size().rename("n_cells")
    n_epi = by_s.apply(lambda g: int((g["celltype_level1"] == "Epithelial").sum())).rename("n_epi")
    sq = pd.concat([n_total, n_epi], axis=1).fillna(0).astype(int)
    sq["passes"] = (sq["n_cells"] >= MIN_CELLS_TOTAL) & (sq["n_epi"] >= MIN_CELLS_EPI)
    keep_samples = sorted(sq.index[sq["passes"]].tolist())
    print(f"    Samples passing (>= {MIN_CELLS_TOTAL} total / >= {MIN_CELLS_EPI} epi): "
          f"{len(keep_samples)} / {len(sq)}", flush=True)
    cohort = cohort[cohort["sample_id"].isin(keep_samples)]

    print("\n[4] Streaming X via single sequential pass...", flush=True)
    barcode_to_row = pd.Series(np.arange(n_cells_total, dtype=np.int64), index=obs.index)
    cohort_rows = barcode_to_row.loc[cohort.index].to_numpy()
    s2col = {s: i for i, s in enumerate(keep_samples)}
    cohort_sample_col = np.array([s2col[s] for s in cohort["sample_id"].to_numpy()], dtype=np.int32)
    cohort_is_epi = (cohort["celltype_level1"] == "Epithelial").to_numpy()

    sample_for_row = np.full(n_cells_total, -1, dtype=np.int32)
    sample_for_row[cohort_rows] = cohort_sample_col
    epi_for_row = np.zeros(n_cells_total, dtype=bool)
    epi_for_row[cohort_rows[cohort_is_epi]] = True

    print("    Reading indptr...", flush=True)
    indptr_all = f["X"]["indptr"][:]
    nnz_total = int(indptr_all[-1])
    print(f"    indptr: {len(indptr_all):,} entries, nnz_total = {nnz_total:,}", flush=True)

    bulk_pb = np.zeros((n_genes, len(keep_samples)), dtype=np.float64)
    epi_pb = np.zeros((n_genes, len(keep_samples)), dtype=np.float64)
    data_ds = f["X"]["data"]
    indices_ds = f["X"]["indices"]

    chunk_starts = [0]
    cur = 0
    while True:
        target = indptr_all[cur] + NNZ_PER_CHUNK
        next_row = int(np.searchsorted(indptr_all, target, side="right")) - 1
        if next_row <= cur:
            next_row = cur + 1
        next_row = min(next_row, n_cells_total)
        chunk_starts.append(next_row)
        if next_row >= n_cells_total:
            break
        cur = next_row
    chunk_starts = list(dict.fromkeys(chunk_starts))
    n_chunks = len(chunk_starts) - 1
    print(f"    Plan: {n_chunks} sequential chunks", flush=True)

    t_scan = time.time()
    for ci in range(n_chunks):
        r0, r1 = chunk_starts[ci], chunk_starts[ci + 1]
        d0, d1 = int(indptr_all[r0]), int(indptr_all[r1])
        n_chunk_rows = r1 - r0
        t_chunk = time.time()

        data_chunk = data_ds[d0:d1]
        indices_chunk = indices_ds[d0:d1]
        local_indptr = (indptr_all[r0:r1 + 1] - d0).astype(np.int64)
        sp_chunk = sps.csr_matrix((data_chunk, indices_chunk, local_indptr),
                                  shape=(n_chunk_rows, n_genes))

        chunk_sample = sample_for_row[r0:r1]
        cohort_mask_c = chunk_sample >= 0
        onehot_bulk = None
        if cohort_mask_c.any():
            local_idx = np.where(cohort_mask_c)[0]
            local_samples = chunk_sample[local_idx]
            local_is_epi = epi_for_row[r0:r1][local_idx]
            onehot_bulk = sps.csr_matrix(
                (np.ones(len(local_idx), dtype=np.float32), (local_idx, local_samples)),
                shape=(n_chunk_rows, len(keep_samples)))
            bulk_pb += np.asarray((sp_chunk.T @ onehot_bulk).todense(), dtype=np.float64)
            if local_is_epi.any():
                epi_idx = local_idx[local_is_epi]
                epi_samples = local_samples[local_is_epi]
                onehot_epi = sps.csr_matrix(
                    (np.ones(len(epi_idx), dtype=np.float32), (epi_idx, epi_samples)),
                    shape=(n_chunk_rows, len(keep_samples)))
                epi_pb += np.asarray((sp_chunk.T @ onehot_epi).todense(), dtype=np.float64)
        print(f"    chunk {ci+1}/{n_chunks}: rows {r0:,}-{r1:,} "
              f"({time.time()-t_chunk:.1f}s; cum {time.time()-t_scan:.0f}s)", flush=True)
        del data_chunk, indices_chunk, sp_chunk
        if onehot_bulk is not None:
            del onehot_bulk
        gc.collect()

    print("\n[5] Building per-sample metadata + NMF epitype proportions...", flush=True)
    rows_meta = []
    for sid, grp in cohort.groupby("sample_id"):
        if sid not in s2col:
            continue
        epi_grp = grp[grp["celltype_level1"] == "Epithelial"]
        n_total_s, n_epi_s = len(grp), len(epi_grp)
        epi_counts = defaultdict(int)
        for nmf_lab, short in EPITYPE_LABELS.items():
            epi_counts[short] = int((epi_grp["nmf_label"] == nmf_lab).sum())
        epi_counts["Ciliated"] = int((epi_grp["celltype_level2"] == CILIATED_LABEL).sum())
        epi_counts["Other_epi"] = (n_epi_s - sum(epi_counts[k] for k in
                                                  ["SecA", "Intermediate", "SecB", "Ciliated"]))
        denom = max(n_epi_s, 1)
        rows_meta.append({
            "sample_id": sid, "patient_id": grp["patient_id"].iloc[0],
            "study": grp["study"].iloc[0], "n_cells": n_total_s, "n_epi": n_epi_s,
            "n_secA": epi_counts["SecA"], "n_intermediate": epi_counts["Intermediate"],
            "n_secB": epi_counts["SecB"], "n_ciliated": epi_counts["Ciliated"],
            "n_other_epi": epi_counts["Other_epi"],
            "pct_secA": 100.0 * epi_counts["SecA"] / denom,
            "pct_intermediate": 100.0 * epi_counts["Intermediate"] / denom,
            "pct_secB": 100.0 * epi_counts["SecB"] / denom,
            "pct_ciliated": 100.0 * epi_counts["Ciliated"] / denom,
            "pct_other_epi": 100.0 * epi_counts["Other_epi"] / denom,
        })
    meta_df = pd.DataFrame(rows_meta).set_index("sample_id").loc[keep_samples].reset_index()

    print("\n[6] Writing outputs...", flush=True)
    bulk_df = pd.DataFrame(bulk_pb, index=var_names, columns=keep_samples)
    epi_df = pd.DataFrame(epi_pb, index=var_names, columns=keep_samples)
    bulk_df.index.name = epi_df.index.name = "gene"
    bulk_df.to_csv(os.path.join(OUT_DIR, "pseudobulk_bulk_counts.tsv.gz"),
                   sep="\t", compression="gzip", float_format="%.0f")
    epi_df.to_csv(os.path.join(OUT_DIR, "pseudobulk_epi_counts.tsv.gz"),
                  sep="\t", compression="gzip", float_format="%.0f")
    meta_df.to_csv(os.path.join(OUT_DIR, "pseudobulk_metadata.csv"), index=False)
    pd.Series(var_names).to_csv(os.path.join(OUT_DIR, "gene_list.tsv"), index=False, header=False)
    print(f"    Bulk: {bulk_df.shape}  Epi: {epi_df.shape}  Meta: {meta_df.shape}")
    print(f"\n[done] total wall time: {time.time()-t_start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
