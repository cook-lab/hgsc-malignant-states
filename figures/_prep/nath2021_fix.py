#!/usr/bin/env python3
"""
Nath 2021 (GSE158722) sample mapping fix for the final atlas
============================================================
PURPOSE
    GEO GSE158722 holds one counts file per PATIENT (P01-P24, no P22), with the time points inside each file:
    P01-P09 are longitudinal (cell ids carry "_Time1_", "_Time2_", "_Time3_"), P10-P24 are single pre- or
    post-treatment samples. The curated sheet (atlas_metadata.xlsx) lists the same design as 40 samples
    (P1-1 ... P9-3, then P10 ... P24), but data preparation gave its first 21 names, then P23 and P24, to the 23
    files in order. So in the deposited atlas a "sample" of Nath 2021 is a whole GEO patient, carrying the name,
    patient number and treatment of another sample (for example GEO P04, all three time points, is "P2-2",
    patient 64), and nine curated "patients" are mixtures of different GEO patients.

    Fix: each Nath cell is re-assigned from its own GEO cell id,
    which the per-study object keeps in obs["barcode"] ("<curated name>:<GEO cell id>"): GEO patient Pk, time
    point TimeN -> curated sample "Pk-N" (k <= 9) or "Pk" (k >= 10), and that sample's curated row gives
    sample number, patient number, treatment status, TP53 and BRCA status. Every other column and every other
    study is untouched. The deposited atlas object is not modified; the correction is applied when the figure
    caches are built (figures/_prep/atlas_final_tables.py) and by the scripts that read the atlas directly.

    Checks (asserted): every Nath cell of the final atlas is found in the per-study object by its index; that
    object's labels equal the atlas labels before the fix; every GEO cell id parses to a curated sample; the
    curated design matches GEO (time points per patient; GEO's pre/post label of P10-P24 equals the curated
    treatment).

INPUTS (read-only)
    data_root/2026_final_atlas/pre-processing and integration/processed/nath_2021.h5ad   obs only (h5py)
    data_root/2026_final_atlas/pre-processing and integration/atlas_metadata.xlsx        curated sample sheet
    data_root/2026_final_atlas/pre-processing and integration/raw/nath_2021/GSE158722.cell_annotations.txt.gz
    obj("atlas_final") obs index + Nath labels (cross-check)

OUTPUTS
    output_root/_prep_caches/nath2021_cell_correction.parquet   one row per Nath atlas cell (index = atlas cell id):
        geo_patient, geo_time, sample_num, sample_id, patient_id, treatment_status, TP53_status, BRCA_status
        (atlas column names: sample_id = sample number, sample_num = sample name, as in the atlas obs)
    output_root/_prep_caches/nath2021_sample_map.csv           old atlas sample -> corrected sample, with cells

USAGE
    from figures._prep.nath2021_fix import apply_nath2021_fix
    obs = apply_nath2021_fix(obs)        # obs indexed by atlas cell id, with a "study" column
    # adds obs["sample_id_deposited"]: the pre-fix sample_id, the key of the deposited CopyKAT runs
    python figures/_prep/nath2021_fix.py # build the correction and print the old -> new map

MANUSCRIPT PANEL(S): every atlas panel and table that uses Nath 2021 sample, patient, treatment, TP53 or BRCA labels
    (Figs 1b-e, 1g, 1k, 2b-g; Figures S3, S9; Table S1; Fig 6h atlas rows).
RUNTIME TIER: fast (obs of the per-study object only).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path  # noqa: E402

STUDY = "nath_2021"
PREP = ("data_root", "2026_final_atlas", "pre-processing and integration")
PER_STUDY_H5AD = path(*PREP, "processed", "nath_2021.h5ad")
CURATED_XLSX = path(*PREP, "atlas_metadata.xlsx")
GEO_ANNOT = path(*PREP, "raw", "nath_2021", "GSE158722.cell_annotations.txt.gz")
CACHE = path("output_root", "_prep_caches", "nath2021_cell_correction.parquet")
MAP_CSV = path("output_root", "_prep_caches", "nath2021_sample_map.csv")
FIX_COLS = ["sample_id", "sample_num", "patient_id", "treatment_status", "TP53_status", "BRCA_status"]
FIX_TAG = "nath2021_fix_v2"   # v2: adds sample_id_deposited


def _h5_obs(fpath, cols):
    import h5py
    with h5py.File(fpath, "r") as f:
        o = f["obs"]

        def col(name):
            g = o[name]
            if isinstance(g, h5py.Group):
                cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
                return cats[g["codes"][:]]
            return np.array([x.decode() if isinstance(x, bytes) else x for x in g[:]])
        return pd.DataFrame({c: col(c) for c in cols}, index=pd.Index(col(o.attrs["_index"])).astype(str))


def _na(x):
    return "NA" if pd.isna(x) else str(x)


def build(atlas_nath: pd.DataFrame) -> pd.DataFrame:
    """Correction table for the Nath cells of `atlas_nath` (atlas obs rows of nath_2021, uncorrected)."""
    per = _h5_obs(PER_STUDY_H5AD, ["barcode", "sample_id", "sample_num", "patient_id", "treatment_status"])
    assert atlas_nath.index.isin(per.index).all(), "Nath atlas cells missing from the per-study object"
    per = per.loc[atlas_nath.index]
    # the per-study object stores name in sample_id and number in sample_num; the atlas the other way round
    for a, b in [("sample_num", "sample_id"), ("sample_id", "sample_num"), ("patient_id", "patient_id"),
                 ("treatment_status", "treatment_status")]:
        assert (atlas_nath[a].astype(str).values == per[b].astype(str).values).all(), f"{a} differs"
    geo = per["barcode"].str.split(":", n=1).str[1]
    pat = geo.str.extract(r"^P(\d+)_")[0].astype(int)
    tim = geo.str.extract(r"_Time(\d)_")[0]
    assert tim[pat <= 9].notna().all() and tim[pat >= 10].isna().all(), "time points not where expected"
    name = np.where(pat <= 9, "P" + pat.astype(str) + "-" + tim.fillna(""), "P" + pat.astype(str))

    cur = pd.read_excel(CURATED_XLSX)
    cur = cur[cur["study"] == STUDY].set_index("sample_id")
    assert set(name) <= set(cur.index), f"unmapped GEO samples: {set(name) - set(cur.index)}"
    # curated design == GEO design (time points per longitudinal patient; P10-P24 pre/post)
    ann = pd.read_csv(GEO_ANNOT, sep="\t", usecols=["Sample.name", "time"]).drop_duplicates()
    for k in range(10, 25):
        g = ann[ann["Sample.name"] == f"P{k}"]
        if len(g):
            geo_pre = g["time"].iloc[0] == "Pre-treatment"
            assert geo_pre == (cur.at[f"P{k}", "treatment_status"] == "pre-treatment"), f"P{k} treatment"
    longit = cur.index[cur.index.str.contains("-")]
    n_cur = pd.Series(longit.str.split("-").str[0]).value_counts()
    n_geo = pd.Series(np.unique(name[pat <= 9])).str.split("-").str[0].value_counts()
    assert all(n_geo[p] <= n_cur[p] for p in n_geo.index), "more GEO time points than curated samples"

    c = cur.loc[name]
    fix = pd.DataFrame({
        "geo_patient": "P" + pat.astype(str).str.zfill(2).values,
        "geo_time": np.where(tim.isna(), "", "Time" + tim.fillna("")),
        "sample_num": name, "sample_id": c["sample_num"].astype(int).astype(str).values,
        "patient_id": c["patient_id"].astype(int).astype(str).values,
        "treatment_status": [_na(v) for v in c["treatment_status"]],
        "TP53_status": [_na(v) for v in c["TP53_status"]], "BRCA_status": [_na(v) for v in c["BRCA_status"]],
    }, index=atlas_nath.index)
    return fix


def correction(atlas_nath: pd.DataFrame | None = None) -> pd.DataFrame:
    """Cached correction table (built from the uncorrected atlas Nath rows when absent)."""
    if os.path.exists(CACHE):
        return pd.read_parquet(CACHE)
    if atlas_nath is None:
        raise FileNotFoundError(f"{CACHE} missing: run python figures/_prep/nath2021_fix.py")
    fix = build(atlas_nath)
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    fix.to_parquet(CACHE)
    m = (pd.DataFrame({"old_sample_name": atlas_nath["sample_num"].astype(str).values,
                       "old_patient": atlas_nath["patient_id"].astype(str).values,
                       "old_treatment": atlas_nath["treatment_status"].astype(str).values})
         .assign(**{k: fix[k].values for k in ["geo_patient", "geo_time", "sample_num", "sample_id", "patient_id",
                                                 "treatment_status"]})
         .groupby(["old_sample_name", "old_patient", "old_treatment", "geo_patient", "geo_time", "sample_num",
                   "sample_id", "patient_id", "treatment_status"], dropna=False).size().rename("cells").reset_index())
    m.to_csv(MAP_CSV, index=False)
    return fix


def apply_nath2021_fix(obs: pd.DataFrame) -> pd.DataFrame:
    """Return `obs` with the six Nath 2021 label columns corrected (other rows and columns untouched).
    Adds sample_id_deposited (the atlas sample_id before the fix, all studies) for joins to deposited per-sample
    runs such as CopyKAT. Idempotent: obs.attrs[FIX_TAG] marks a corrected frame."""
    if obs.attrs.get(FIX_TAG):
        return obs
    obs = obs.copy()
    is_n = obs["study"].astype(str).values == STUDY
    if "sample_id" in obs.columns and "sample_id_deposited" not in obs.columns:
        # the deposited per-sample runs (CopyKAT, keyed by atlas sample_id) keep their key: for Nath 2021 each run is
        # one GEO patient with its time points pooled; CNV scripts join runs to cells by this column
        obs["sample_id_deposited"] = obs["sample_id"].astype(str).astype("category")
    fix = correction(obs.loc[is_n] if is_n.any() else None)
    rows = obs.index[is_n]
    assert rows.isin(fix.index).all(), "Nath cells without a correction"
    for c in [c for c in FIX_COLS if c in obs.columns]:
        new = fix.loc[rows, c].values
        if isinstance(obs[c].dtype, pd.CategoricalDtype):
            cats = list(obs[c].cat.categories) + sorted(set(new) - set(obs[c].cat.categories))
            col = obs[c].cat.set_categories(cats)
            col[is_n] = new
            obs[c] = col.cat.remove_unused_categories()
        else:
            obs.loc[rows, c] = new
    obs.attrs[FIX_TAG] = True
    return obs


if __name__ == "__main__":
    from figures._prep.atlas_final_tables import _read_final_obs
    o = _read_final_obs()
    n = o[o["study"].astype(str) == STUDY]
    if os.path.exists(CACHE):
        os.remove(CACHE)
    fix = correction(n)
    print(f"Nath cells: {len(n):,}; samples {n['sample_num'].nunique()} -> {fix['sample_num'].nunique()}, "
          f"patients {n['patient_id'].nunique()} -> {fix['patient_id'].nunique()}")
    print(pd.read_csv(MAP_CSV).to_string(index=False))
    o2 = apply_nath2021_fix(o)
    print(f"atlas: samples {o['sample_id'].nunique()} -> {o2['sample_id'].nunique()}, "
          f"patients {o['patient_id'].nunique()} -> {o2['patient_id'].nunique()}")
    print("treatment (Nath cells):", o2.loc[o2.study == STUDY, "treatment_status"].value_counts().to_dict())
