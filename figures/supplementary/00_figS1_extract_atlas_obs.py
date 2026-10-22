"""
Figure S1 (data) — per-cell QC metrics of the concatenated atlas, before and after the UMI / gene filter
=======================================================================================================

PURPOSE
    Reads the two concatenated atlas objects chunk-wise with h5py (X is never loaded into RAM) and
    covers ALL cells:
    prefilter : pre-processing and integration/processed/atlas_concat_counts_only_X.h5ad
                (2,731,632 cells, raw concatenated counts, no QC) -> total_counts, n_genes_by_counts
                computed from X over every cell (Figure S1a, b).
    postfilter: pre-processing and integration/processed/atlas_concatenated_filtered.h5ad
                (2,326,532 cells: after the UMI/gene filter, before the 0.25 Scrublet cutoff) obs ->
                doublet_score_scrublet (only the doublet columns are used, by Figure S1c). This object
                still holds the cells that the 0.25 cutoff removed, so panel c shows them.

INPUTS
    the two h5ad files above, in the folder set on the ATLAS line below

OUTPUTS
    $OUTPUT_ROOT/_prep_caches/atlas_obs_prefilter.parquet   (study, sample_id, patient_id, total_counts,
                                                             n_genes_by_counts)
    $OUTPUT_ROOT/_prep_caches/atlas_obs_postfilter.parquet  (study, sample_id, patient_id,
                                                             doublet_score_scrublet, doublet_scrublet)
    Both are read by 01_figS1_atlas_qc_metrics_by_study.py.

MANUSCRIPT PANEL(S)
    Figure S1a-c (data).

RUNTIME TIER
    slow (streams the full CSR X of the 2.7M-cell object). OUTPUT_ROOT must be set in the environment.
"""
import os, time
import anndata as ad
import h5py
import numpy as np
import pandas as pd

ATLAS = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas/pre-processing and integration/processed")
RAW = os.path.join(ATLAS, "atlas_concat_counts_only_X.h5ad")
POST = os.path.join(ATLAS, "atlas_concatenated_filtered.h5ad")
OUT = os.path.join(os.environ["OUTPUT_ROOT"], "_prep_caches")
os.makedirs(OUT, exist_ok=True)


def read_obs_cols(f, cols):
    o = f["obs"]
    idx = ad.io.read_elem(o[o.attrs["_index"]])
    df = pd.DataFrame({c: ad.io.read_elem(o[c]) for c in cols}, index=pd.Index(np.asarray(idx).astype(str)))
    return df


def extract_postfilter():
    t = time.time()
    with h5py.File(POST, "r") as f:
        df = read_obs_cols(f, ["study", "sample_id", "patient_id", "doublet_score_scrublet", "doublet_scrublet"])
    df["doublet_score_scrublet"] = pd.to_numeric(df["doublet_score_scrublet"], errors="coerce")
    df.to_parquet(os.path.join(OUT, "atlas_obs_postfilter.parquet"))
    print(f"postfilter: {len(df):,} cells ({time.time()-t:.0f}s)", flush=True)


def extract_prefilter(chunk=250_000_000):
    t = time.time()
    with h5py.File(RAW, "r") as f:
        df = read_obs_cols(f, ["study", "sample_id", "patient_id"])
        X = f["X"]
        indptr = X["indptr"][:]
        data = X["data"]
        n = len(indptr) - 1
        tot = np.zeros(n, dtype=np.int64)
        ng = np.zeros(n, dtype=np.int64)
        r0 = 0
        while r0 < n:
            s = indptr[r0]
            r1 = int(np.searchsorted(indptr, s + chunk, side="right") - 1)
            r1 = max(r1, r0 + 1)
            r1 = min(r1, n)
            e = indptr[r1]
            d = data[s:e].astype(np.int64)
            cs = np.concatenate([[0], np.cumsum(d)])
            cp = np.concatenate([[0], np.cumsum(d > 0)])
            lo = indptr[r0:r1] - s
            hi = indptr[r0 + 1:r1 + 1] - s
            tot[r0:r1] = cs[hi] - cs[lo]
            ng[r0:r1] = cp[hi] - cp[lo]
            print(f"  rows {r0:,}-{r1:,} / {n:,}  ({time.time()-t:.0f}s)", flush=True)
            r0 = r1
    df["total_counts"] = tot
    df["n_genes_by_counts"] = ng
    df.to_parquet(os.path.join(OUT, "atlas_obs_prefilter.parquet"))
    print(f"prefilter: {len(df):,} cells ({time.time()-t:.0f}s)", flush=True)


if __name__ == "__main__":
    extract_postfilter()
    extract_prefilter()
