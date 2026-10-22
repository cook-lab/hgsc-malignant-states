#!/usr/bin/env python3
"""
Coexistence of SecA, Intermediate and SecB cells within patients and within CopyKAT clones (Discussion)
====================================================================================================
PURPOSE
    Backs the Discussion sentence "93.5% of patients with ≥ 100 epithelial cells contain SecA, Intermediate, and SecB
    cells, all existing within individual clones" (patient row: 144 / 154 patients). The clone row (CopyKAT clones that
    contain all three states) supports the second half of the sentence; its count is not quoted in the text.

    * Patients: per patient, all epithelial cells of the final atlas (corrected labels, load_epithelial()). Eligible =
      ≥ 100 epithelial cells (SecA, Intermediate, SecB and ciliated). A patient "contains" a state with ≥ 5 cells of it.
    * Clones: the CopyKAT subclones of Figure S7 (≥ 20 aneuploid epithelial cells; within_clone_coexistence.csv written by
      figures/supplementary/13_figS7_cnv_clone_validation.py). A clone contains a state with ≥ 1 cell of it.

INPUTS   load_epithelial() (output_root/_prep_caches/atlas_final_epithelial.parquet; Nath 2021 fix applied)
         output_root/_prep_caches/sf7_cnv/within_clone_coexistence.csv
OUTPUTS  output_root/_prep_caches/text_state_coexistence.csv (one row per count, with numerator and denominator)
MANUSCRIPT PANEL(S): none (Discussion text).
RUNTIME TIER: fast (seconds).
USAGE (repo root): OUTPUT_ROOT=<atlas figure outputs> python figures/_prep/text_state_coexistence.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

MIN_EPI, MIN_STATE = 100, 5
STATES = ["SecA", "Intermediate", "SecB"]
COEX_CSV = path("output_root", "_prep_caches", "sf7_cnv", "within_clone_coexistence.csv")
OUT = path("output_root", "_prep_caches", "text_state_coexistence.csv")

epi = load_epithelial(columns=["patient_id", "state"])
tab = epi.groupby(["patient_id", epi["state"].astype(str)], observed=True).size().unstack(fill_value=0)
for s in STATES:
    if s not in tab:
        tab[s] = 0
eligible = tab[tab.sum(axis=1) >= MIN_EPI]
n_pat = int((eligible[STATES] >= MIN_STATE).all(axis=1).sum())

clones = pd.read_csv(COEX_CSV)
n_clone = int((clones[["n_secA", "n_transitioning", "n_secB"]] >= 1).all(axis=1).sum())

res = pd.DataFrame([
    {"unit": "patient", "rule": f">= {MIN_EPI} epithelial cells; >= {MIN_STATE} cells of each of SecA, Intermediate, SecB",
     "n_with_all_three": n_pat, "n_eligible": len(eligible), "percent": round(100 * n_pat / len(eligible), 1)},
    {"unit": "CopyKAT clone (Figure S7)", "rule": ">= 20 aneuploid epithelial cells; >= 1 cell of each state",
     "n_with_all_three": n_clone, "n_eligible": len(clones), "percent": round(100 * n_clone / len(clones), 1)},
])
os.makedirs(os.path.dirname(OUT), exist_ok=True)
res.to_csv(OUT, index=False)
print(res.to_string(index=False))
