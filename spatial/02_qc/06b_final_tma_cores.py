#!/usr/bin/env python3
"""
Final TMA core set: the one set of TMA cores behind every manuscript panel
==========================================================================
PURPOSE
    Author decision (2026-10-19): every TMA analysis uses the same cores, chosen the same way in TMA A and TMA B.
      1. Tumour cores (on-core cells of a TMA patient; the FTE cores are a separate arm).
      2. Core QC: 05b_core_qc.R run once per array (TMA A; TMA B blocks G and H) with the same method on the
         final SingleR labels. Cores with >= 2 of its 5 flags are excluded; pass and review cores are kept.
      3. Patients with primary HGSC (central review; unreviewed histology kept) who were treatment-naive at
         sampling (no neoadjuvant chemotherapy). Same rules and sources as Table S8 and the Fig 7 survival cohort.
      4. Patients with >= 500 cells on their retained cores.
    A core is final when it passes 1-3 and belongs to a patient who passes 4. Panel-specific minimums (for example
    >= 500 secretory cells per core in Fig 4d) are applied on top of this set by each panel. TMA B pools blocks G
    and H: a patient lies in one block only (asserted).

INPUTS
    output_root/07_core_qc/core_qc_summary_singler_{tma,COEUR_TMA_G,COEUR_TMA_H}.csv   (05b_core_qc.R, 2026-10-19)
    output_root/singler_recompute/{tma,COEUR_TMA_G,COEUR_TMA_H}_singler_cells.csv.gz      (labels per cell)
    data_root/202609_final_hgsc_spatial/data/clinical data/chum_tma_clinical_data_updated.xlsx                     (TMA A)
    data_root/202609_final_hgsc_spatial/data/clinical data/coeur_clinical_tmagh.xlsx, sheet "clinical data"         (TMA B)

OUTPUTS (raw identifiers: lab-internal, never committed)
    output_root/07_core_qc/final_tma_cores.csv      one row per profiled tumour core: QC flags, clinical flags, `final`
    output_root/07_core_qc/final_tma_patients.csv   one row per TMA patient: cores and cells at each step, `final`,
                                                    `n_final_cores` (the Fig 7 survival cohort needs >= 2)
    output_root/07_core_qc/final_tma_counts.csv     counts by TMA at each step (printed; quoted in the manuscript)

MANUSCRIPT PANEL(S): the TMA arm of Figs 4a, 4d, 4g, 4i, 5b, 5j, 6g, 7a-f, Figures S12-S15; Tables S8, S11.
RUNTIME TIER: fast (seconds).

USAGE
    OUTPUT_ROOT=<outputs> python spatial/02_qc/06b_final_tma_cores.py
"""
import os
import sys

import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO)
from config.config import path  # noqa: E402

CLIN = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial", "data", "clinical data")
QC_DIR = path("output_root", "07_core_qc")
CELLS = path("output_root", "singler_recompute")
ARRAYS = [("tma", "A", "A"), ("COEUR_TMA_G", "B", "G"), ("COEUR_TMA_H", "B", "H")]   # (file stem, TMA, block)
MIN_PATIENT_CELLS = 500
SEC = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]


def num(x):
    return pd.to_numeric(x, errors="coerce")


# ------------------------------------------------------------------ 1-2. cores, core QC (per array, SingleR labels)
cores = []
for stem, tma, block in ARRAYS:
    qc = pd.read_csv(os.path.join(QC_DIR, f"core_qc_summary_singler_{stem}.csv"),
                     dtype={"core_id_raw": str, "patient_id_raw": str})
    d = pd.read_csv(os.path.join(CELLS, f"{stem}_singler_cells.csv.gz"),
                    usecols=["core_id", "patient_id", "sample_type", "cell_label_singler"],
                    dtype={"core_id": str, "patient_id": str})
    d = d[(d.sample_type == "tumour") & d.core_id.notna() & (d.core_id != "Off core")]
    n_sec = d[d.cell_label_singler.isin(SEC)].groupby("core_id").size()
    flags = [c for c in qc.columns if c.startswith("flag_")]
    c = qc[["core_id_raw", "patient_id_raw", "n_cells", "n_flags", "recommendation"] + flags].copy()
    c.insert(0, "block", block); c.insert(0, "tma", tma)
    c["n_secretory"] = c.core_id_raw.map(n_sec).fillna(0).astype(int)
    # the QC table and the cell table describe the same cores and cells
    assert set(c.core_id_raw) == set(d.core_id.unique()), stem
    assert (c.set_index("core_id_raw").n_cells == d.groupby("core_id").size().reindex(c.core_id_raw).values).all(), stem
    # a core that no patient was assigned to (one in block H) is QC'd with its array but enters no analysis
    no_pat = c.patient_id_raw.isna()
    if no_pat.any():
        print(f"{stem}: {int(no_pat.sum())} core(s) without a patient ({int(c.n_cells[no_pat].sum())} cells) left out")
    cores.append(c[~no_pat])
cores = pd.concat(cores, ignore_index=True)
cores["patient_key"] = cores.tma + ":" + cores.patient_id_raw
blocks = cores.groupby("patient_key").block.nunique()
assert (blocks == 1).all(), "a TMA B patient spans blocks G and H"
cores["core_qc_pass"] = cores.recommendation != "exclude"

# ------------------------------------------------------------------ 3. clinical flags (as Table S8 / Fig 7)
chum = pd.read_excel(os.path.join(CLIN, "chum_tma_clinical_data_updated.xlsx"))
a = pd.DataFrame({"patient_key": "A:" + chum.patient_id.astype(str),
                  "treatment": chum.treatment_status.map({"naïve": "Treatment-naive",
                                                          "post-chemotherapy": "Post-chemotherapy (NACT)"}),
                  "histology": pd.Series("HGSC", index=chum.index).where(chum.stage.notna())})
co = pd.read_excel(os.path.join(CLIN, "coeur_clinical_tmagh.xlsx"), sheet_name="clinical data")
co = co[co.Bloc.isin(["G", "H"])]
HIST_CO = {1: "HGSC", 2: "Endometrioid", 3: "Low-grade serous", 4: "Clear cell", 5: "Mucinous", 6: "Mixed"}
b = pd.DataFrame({"patient_key": "B:" + co["COEUR-TFRI"].astype(str),
                  "treatment": num(co.Post_Chemo_Code).map({0: "Treatment-naive", 1: "Post-chemotherapy (NACT)"}),
                  "histology": num(co.Histopathologyrevision).map(HIST_CO)})
clin = pd.concat([a, b], ignore_index=True)
assert not clin.patient_key.duplicated().any()
cores = cores.merge(clin, on="patient_key", how="left")
cores["nact"] = cores.treatment.eq("Post-chemotherapy (NACT)")
cores["non_hgsc"] = cores.histology.notna() & cores.histology.ne("HGSC")
cores["clinical_pass"] = ~cores.nact & ~cores.non_hgsc

# ------------------------------------------------------------------ 4. patients with >= 500 cells on retained cores
kept = cores[cores.core_qc_pass & cores.clinical_pass]
cells_kept = kept.groupby("patient_key").n_cells.sum()
cores["patient_cells_retained"] = cores.patient_key.map(cells_kept).fillna(0).astype(int)
cores["final"] = cores.core_qc_pass & cores.clinical_pass & (cores.patient_cells_retained >= MIN_PATIENT_CELLS)
cores["reason_excluded"] = np.select(
    [~cores.core_qc_pass, cores.nact, cores.non_hgsc, cores.patient_cells_retained < MIN_PATIENT_CELLS],
    ["core QC (>= 2 flags)", "patient: neoadjuvant chemotherapy", "patient: not HGSC on central review",
     "patient: < 500 cells on retained cores"], default="")
cores.loc[cores.final, "reason_excluded"] = ""

pat = cores.groupby("patient_key").agg(
    tma=("tma", "first"), block=("block", "first"), patient_id_raw=("patient_id_raw", "first"),
    treatment=("treatment", "first"), histology=("histology", "first"),
    n_cores_profiled=("core_id_raw", "size"), n_cores_qc_pass=("core_qc_pass", "sum"),
    n_final_cores=("final", "sum"), cells_profiled=("n_cells", "sum"),
    cells_final=("n_cells", lambda s: int(s[cores.loc[s.index, "final"]].sum())),
    secretory_final=("n_secretory", lambda s: int(s[cores.loc[s.index, "final"]].sum())),
    nact=("nact", "first"), non_hgsc=("non_hgsc", "first")).reset_index()
pat["final"] = pat.n_final_cores > 0
pat["two_final_cores"] = pat.n_final_cores >= 2

cores.to_csv(os.path.join(QC_DIR, "final_tma_cores.csv"), index=False)
pat.to_csv(os.path.join(QC_DIR, "final_tma_patients.csv"), index=False)

# ------------------------------------------------------------------ counts by step
rows = []
for name, sel in [("profiled tumour cores", cores.index == cores.index),
                  ("after core QC", cores.core_qc_pass),
                  ("after core QC, HGSC, treatment-naive", cores.core_qc_pass & cores.clinical_pass),
                  ("final (patients with >= 500 cells on retained cores)", cores.final)]:
    s = cores[sel]
    for t in ["A", "B", "A+B"]:
        u = s if t == "A+B" else s[s.tma == t]
        rows.append({"step": name, "TMA": t, "patients": u.patient_key.nunique(), "cores": len(u),
                     "cells": int(u.n_cells.sum()), "secretory_cells": int(u.n_secretory.sum())})
for t in ["A", "B", "A+B"]:
    u = pat if t == "A+B" else pat[pat.tma == t]
    rows.append({"step": "final patients with >= 2 final cores", "TMA": t,
                 "patients": int(u.two_final_cores.sum()), "cores": int(u.n_final_cores[u.two_final_cores].sum()),
                 "cells": int(u.cells_final[u.two_final_cores].sum()),
                 "secretory_cells": int(u.secretory_final[u.two_final_cores].sum())})
counts = pd.DataFrame(rows)
counts.to_csv(os.path.join(QC_DIR, "final_tma_counts.csv"), index=False)
pd.set_option("display.width", 200)
print(counts.to_string(index=False))
print("\nCores excluded, by first reason:")
print(cores[~cores.final].groupby(["tma", "reason_excluded"]).size().to_string())
print("\nPatients dropped by clinical rules (profiled):",
      pat.groupby("tma")[["nact", "non_hgsc"]].sum().to_dict())
