#!/usr/bin/env python3
"""
Fig 2f prep — consensusOV subtype calls from the final atlas
============================================================
PURPOSE
    Builds the whole-sample pseudobulks and per-sample state composition behind
    Fig 2f (figures/figure2/03_fig2f_subtype_state_boxes.py), from the final atlas
    only, and runs one consensusOV TCGA-subtype call (seed 42, fig2f_consensusov.R).
    The Fig 2f subtypes are the 25-seed consensus calls of fig2f_consensusov_multiseed.R,
    run on this script's pseudobulk_counts.tsv.gz; 03_fig2f_... maps them onto the
    joined table in place of the single-seed column.

      * Samples: pre-treatment adnexal (anatomic_site == "adnexa",
        treatment_status == "pre-treatment") with >= 500 cells and >= 100
        epithelial cells. Expected 96 samples / 77 patients (asserted).
      * Whole-sample pseudobulk: raw counts (obj atlas_final, layers/counts)
        summed over ALL cells of each sample. Only the rows of the selected
        samples are read (contiguous runs, chunked); X is never loaded.
      * fig2f_consensusov.R: edgeR TMM -> log2(CPM + 1) -> Entrez IDs
        (org.Hs.eg.db) -> consensusOV::get.subtypes(method = "consensusOV").
      * Epithelial composition per sample from load_epithelial()['state']
        (% of all epithelial cells of the sample: SecA / Intermediate / SecB /
        Ciliated).

INPUTS
    obj("atlas_final")  (layers/counts rows of the 96 samples)
    figures/_prep/atlas_final_tables.py  (load_obs, load_epithelial)

OUTPUTS (output_root/figure2/fig2f_consensusov/)
    pseudobulk_counts.tsv.gz        genes x samples raw summed counts
    pseudobulk_metadata.csv         per-sample n's + epithelial composition
    consensusov_calls.csv           single-seed calls, written by fig2f_consensusov.R
                                    (with gene_id_mapping.csv, the symbol -> Entrez map)
    fig2f_per_sample_joined.csv     composition + single-seed subtype (input of
                                    figures/figure2/03_fig2f_subtype_state_boxes.py)

    --cohort all_tumours (every non-healthy sample; written to fig2f_consensusov_all/) is a
    sensitivity variant that is not in the paper; the published run uses the default.

MANUSCRIPT PANEL(S): Fig 2f (input).
RUNTIME TIER: moderate (reads ~0.6 M rows of the counts layer; ~10-20 min on SMB).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial, load_obs  # noqa: E402

import argparse  # noqa: E402
_ap = argparse.ArgumentParser()
_ap.add_argument("--cohort", choices=["primary_pretreatment", "all_tumours"], default="primary_pretreatment",
                 help="primary_pretreatment = adnexal pre-treatment samples (Fig 2f, 96 samples); "
                      "all_tumours = every non-healthy sample (primary, metastasis, ascites; any "
                      "treatment status) — sensitivity variant, written to fig2f_consensusov_all/")
COHORT = _ap.parse_args().cohort
OUT_DIR = path("output_root", "figure2",
               "fig2f_consensusov" if COHORT == "primary_pretreatment" else "fig2f_consensusov_all")
os.makedirs(OUT_DIR, exist_ok=True)
R_SCRIPT = Path(__file__).with_suffix(".R")
MIN_CELLS, MIN_EPI = 500, 100
EXPECTED = (96, 77)
CHUNK_ROWS = 20_000


# ---------- 1. Samples ----------
obs = load_obs(columns=["sample_id", "patient_id", "study", "anatomic_site", "metastatic_site",
                        "treatment_status"])
epi = load_epithelial(columns=["sample_id", "state"])
obs["sample_id"] = obs["sample_id"].astype(str)
epi["sample_id"] = epi["sample_id"].astype(str)
if COHORT == "primary_pretreatment":
    cand = obs[(obs.anatomic_site == "adnexa") & (obs.treatment_status == "pre-treatment")]
else:
    cand = obs[obs.metastatic_site.astype(str).isin(["primary", "metastasis", "ascites"])]
n_cells = cand.groupby("sample_id").size()
n_epi = epi.groupby("sample_id").size().reindex(n_cells.index).fillna(0).astype(int)
samples = sorted(n_cells[(n_cells >= MIN_CELLS) & (n_epi >= MIN_EPI)].index)
n_pat = cand.loc[cand.sample_id.isin(samples), "patient_id"].nunique()
print(f"[1] samples: {len(samples)}, patients: {n_pat}")
if COHORT == "primary_pretreatment":
    assert (len(samples), n_pat) == EXPECTED, (len(samples), n_pat)

# ---------- 2. Epithelial composition (loader) ----------
comp = pd.crosstab(epi.loc[epi.sample_id.isin(samples), "sample_id"],
                   epi.loc[epi.sample_id.isin(samples), "state"])
comp = comp.reindex(index=samples, columns=["SecA", "Intermediate", "SecB", "Ciliated"], fill_value=0)
first = cand[cand.sample_id.isin(samples)].groupby("sample_id")[
    ["patient_id", "study", "metastatic_site", "treatment_status", "anatomic_site"]].first()
meta = pd.DataFrame({
    "sample_id": samples,
    "patient_id": first.loc[samples, "patient_id"].astype(str).values,
    "study": first.loc[samples, "study"].astype(str).values,
    "site": first.loc[samples, "metastatic_site"].astype(str).values,
    "anatomic_site": first.loc[samples, "anatomic_site"].astype(str).values,
    "treatment_status": first.loc[samples, "treatment_status"].astype(str).values,
    "n_cells": n_cells.loc[samples].values,
    "n_epi": comp.sum(axis=1).values,
    "n_secA": comp["SecA"].values, "n_trans": comp["Intermediate"].values,
    "n_secB": comp["SecB"].values, "n_ciliated": comp["Ciliated"].values,
})
for k in ("secA", "trans", "secB", "ciliated"):
    meta[f"pct_{k}"] = 100 * meta[f"n_{k}"] / meta["n_epi"]
meta.to_csv(os.path.join(OUT_DIR, "pseudobulk_metadata.csv"), index=False)

# ---------- 3. Whole-sample pseudobulk from layers/counts ----------
counts_out = os.path.join(OUT_DIR, "pseudobulk_counts.tsv.gz")
if not os.path.exists(counts_out):
    import h5py
    from anndata.io import read_elem
    code_of = {s: i for i, s in enumerate(samples)}
    n_outside = int((obs.sample_id.isin(samples) & ~obs.index.isin(cand.index)).sum())
    assert n_outside == 0, f"{n_outside} cells of the selected samples are not adnexa/pre-treatment"
    row_code = obs["sample_id"].map(code_of).fillna(-1).astype(np.int64).values
    pos = np.where(row_code >= 0)[0]
    assert len(pos) == int(n_cells.loc[samples].sum())
    with h5py.File(obj("atlas_final"), "r") as f:
        genes = np.asarray(read_elem(f["var"][f["var"].attrs["_index"]])).astype(str)
        idx_chk = np.asarray(read_elem(f["obs"][f["obs"].attrs["_index"]])).astype(str)
        assert (idx_chk == obs.index.values).all(), "loader obs order != atlas_final obs order"
        G, S = len(genes), len(samples)
        cnt = f["layers/counts"]
        indptr = cnt["indptr"][:]
        acc = np.zeros(S * G, dtype=np.float64)
        breaks = np.where(np.diff(pos) != 1)[0]
        starts = np.r_[pos[0], pos[breaks + 1]]
        ends = np.r_[pos[breaks], pos[-1]] + 1
        print(f"[3] reading {len(pos):,} rows in {len(starts)} contiguous runs ...", flush=True)
        for r0, r1 in zip(starts, ends):
            for c0 in range(r0, r1, CHUNK_ROWS):
                c1 = min(c0 + CHUNK_ROWS, r1)
                a, b = indptr[c0], indptr[c1]
                data = cnt["data"][a:b].astype(np.float64)
                ind = cnt["indices"][a:b].astype(np.int64)
                rc = np.repeat(row_code[c0:c1], np.diff(indptr[c0:c1 + 1]))
                acc += np.bincount(rc * G + ind, weights=data, minlength=S * G)
            print(f"    rows {r0:,}-{r1:,} done", flush=True)
    pb = pd.DataFrame(np.rint(acc.reshape(S, G)).astype(np.int64).T, index=genes, columns=samples)
    pb.index.name = "gene"
    pb.to_csv(counts_out, sep="\t", compression="gzip")
    print(f"    pseudobulk: {pb.shape[0]} genes x {pb.shape[1]} samples -> {counts_out}")

# ---------- 4. consensusOV (R) ----------
calls_csv = os.path.join(OUT_DIR, "consensusov_calls.csv")
subprocess.run(["Rscript", str(R_SCRIPT), counts_out, calls_csv], check=True)
calls = pd.read_csv(calls_csv, dtype={"sample_id": str})
calls["sample_id"] = calls["sample_id"].astype(str)

# ---------- 5. Joined table (input of figures/figure2/03_fig2f_subtype_state_boxes.py) ----------
joined = meta.merge(calls[["sample_id", "subtype", "margin_top1_top2"]], on="sample_id")
joined = joined.rename(columns={"subtype": "subtype_bulk", "margin_top1_top2": "margin_bulk"})
assert len(joined) == len(samples)
joined.to_csv(os.path.join(OUT_DIR, "fig2f_per_sample_joined.csv"), index=False)

if COHORT != "primary_pretreatment":
    print(f"[5] subtype counts: {joined.subtype_bulk.value_counts().to_dict()}")
    sys.exit(0)
print(f"[5] subtype counts new: {joined.subtype_bulk.value_counts().to_dict()}")
