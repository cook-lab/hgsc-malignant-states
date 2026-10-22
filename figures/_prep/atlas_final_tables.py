#!/usr/bin/env python3
"""
Single source of atlas data for Figures 1–2 and Supplementary Figures 1–9
=========================================================================
The Fig 1–2 / Figure S2–S9 scripts load atlas cell metadata through this module.

Sources, all keyed by cell barcode:
  * obj("atlas_final")  hgsc_atlas_final.h5ad: obs metadata and obsm["X_umap"]
    (read with h5py, two dozen obs columns only; X is never loaded).
  * 11d_epithelial_nmf/11d_nmf_usage.csv: per-cell usage of the 10 epithelial
    NMF factors (575,366 epithelial cells; deposited model, not re-fitted).
  * Epithelial-only UMAP: the coordinates (UMAP1/UMAP2) of the deposited
    epithelial embedding, stored in fig_secretory_polarization/data/meta.parquet.
    The object that held them (hgsc_atlas_epithelial.h5ad, obsm X_umap_local) is
    no longer on disk; ONLY these two columns are taken from that file. All
    metadata come from the final atlas.

Nath 2021 sample mapping: build() applies figures/_prep/nath2021_fix.py, which re-assigns each Nath
cell to its GEO patient and time point (sample, patient, treatment, TP53, BRCA); the deposited object is not changed.
Caches without the fix flag are rebuilt on first load.

Epithelial state: among non-ciliated epithelial cells, Factor 2
usage < p50 -> SecA, p50–p75 -> Intermediate, >= p75 -> SecB; ciliated cells are
those whose celltype_level2 contains "iliat". Expected 284,344 / 142,172 /
142,172 / 6,678 (asserted).

Usage
    from figures._prep.atlas_final_tables import load_obs, load_epithelial
    obs = load_obs()          # all 1,980,703 cells: metadata + UMAP1/UMAP2
    epi = load_epithelial()   # 575,366 epithelial cells: metadata + epithelial
                              # UMAP1/UMAP2 + Factor_1..Factor_10 + state
Run directly to (re)build the caches:  python figures/_prep/atlas_final_tables.py
Caches: output_root/_prep_caches/atlas_final_{obs,epithelial}.parquet
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path  # noqa: E402
from figures._prep.nath2021_fix import FIX_TAG, apply_nath2021_fix  # noqa: E402

OBS_COLS = [
    "study", "dataset", "patient_id", "sample_id", "sample_num", "anatomic_site",
    "metastatic_site", "treatment_status", "treatment_response", "stage", "age",
    "histological_subtype", "BRCA_status", "HRD_status", "TP53_status",
    "celltype_level1", "celltype_level2", "celltype_pred", "celltype_probability",
    "celltype_confident", "leiden_r0.3", "total_counts", "n_genes_by_counts",
    "doublet_score_scrublet",
]
NMF_USAGE_CSV = path("data_root", "2026_final_atlas", "output", "11d_epithelial_nmf", "11d_nmf_usage.csv")
EPI_UMAP_SRC = path("data_root", "2026_final_atlas", "output", "fig_secretory_polarization", "data", "meta.parquet")
OBS_CACHE = path("output_root", "_prep_caches", "atlas_final_obs.parquet")
EPI_CACHE = path("output_root", "_prep_caches", "atlas_final_epithelial.parquet")
EXPECTED_STATES = {"SecA": 284_344, "Intermediate": 142_172, "SecB": 142_172, "Ciliated": 6_678}
STATE_ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]


def _read_final_obs() -> pd.DataFrame:
    import h5py
    from anndata.io import read_elem
    with h5py.File(obj("atlas_final"), "r") as f:
        o = f["obs"]
        idx = pd.Index(np.asarray(read_elem(o[o.attrs["_index"]])).astype(str))
        df = pd.DataFrame({c: read_elem(o[c]) for c in OBS_COLS if c in o}, index=idx)
        umap = f["obsm/X_umap"][:]
    df["UMAP1"], df["UMAP2"] = umap[:, 0], umap[:, 1]
    return df


def build() -> None:
    obs = _read_final_obs()
    assert obs.index.is_unique and len(obs) == 1_980_703, len(obs)
    obs = apply_nath2021_fix(obs)   # Nath 2021 cells: GEO patient / time point labels
    obs.to_parquet(OBS_CACHE)
    print(f"[atlas_final_tables] obs: {len(obs):,} cells -> {OBS_CACHE}")

    usage = pd.read_csv(NMF_USAGE_CSV, index_col=0)
    usage.index = usage.index.astype(str)
    assert usage.index.isin(obs.index).all(), "NMF usage barcodes missing from the final atlas"
    epi_umap = pd.read_parquet(EPI_UMAP_SRC, columns=["UMAP1", "UMAP2"])
    epi_umap.index = epi_umap.index.astype(str)
    assert set(epi_umap.index) == set(usage.index), "epithelial UMAP and NMF usage cover different cells"

    meta_cols = [c for c in obs.columns if c not in ("UMAP1", "UMAP2")]
    epi = obs.loc[usage.index, meta_cols].copy()
    epi["UMAP1"] = epi_umap.loc[epi.index, "UMAP1"].values
    epi["UMAP2"] = epi_umap.loc[epi.index, "UMAP2"].values
    epi = epi.join(usage)

    ciliated = epi["celltype_level2"].astype(str).str.contains("iliat").values
    f2 = epi["Factor_2"].values
    p50, p75 = np.percentile(f2[~ciliated], [50, 75])
    epi["state"] = pd.Categorical(
        np.where(ciliated, "Ciliated",
                 np.where(f2 >= p75, "SecB", np.where(f2 >= p50, "Intermediate", "SecA"))),
        categories=STATE_ORDER)
    counts = epi["state"].value_counts().to_dict()
    assert counts == EXPECTED_STATES, f"state counts changed: {counts}"
    epi.attrs.update({"factor2_p50": float(p50), "factor2_p75": float(p75), FIX_TAG: True})
    epi.to_parquet(EPI_CACHE)
    print(f"[atlas_final_tables] epithelial: {len(epi):,} cells, p50={p50:.4f}, p75={p75:.4f}, "
          f"states {counts} -> {EPI_CACHE}")


def _load(cache, columns):
    if os.path.exists(cache):
        df = pd.read_parquet(cache, columns=columns)
        if df.attrs.get(FIX_TAG):
            return df
        print(f"[atlas_final_tables] {os.path.basename(cache)} predates the Nath 2021 fix: rebuilding")
    build()
    df = pd.read_parquet(cache, columns=columns)
    assert df.attrs.get(FIX_TAG), "cache written without the Nath 2021 fix"
    return df


def load_obs(columns=None) -> pd.DataFrame:
    return _load(OBS_CACHE, columns)


def load_epithelial(columns=None) -> pd.DataFrame:
    return _load(EPI_CACHE, columns)


if __name__ == "__main__":
    build()
