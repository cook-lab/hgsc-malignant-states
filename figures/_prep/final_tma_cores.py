"""The final TMA core set for figure and table scripts.

PURPOSE
    Every TMA analysis uses the cores of spatial/02_qc/06b_final_tma_cores.py: tumour
    cores that pass per-array core QC on SingleR labels, from patients with primary HGSC who were treatment-naive at
    sampling and have >= 500 cells on those cores. This module is the one place figure scripts read that set from.

INPUTS
    $FINAL_TMA_CORES, default $OUTPUT_ROOT/07_core_qc/final_tma_cores.csv (raw ids; lab-internal)

USAGE
    from figures._prep.final_tma_cores import final_core_ids, final_cores, final_patients
    keep = final_core_ids("COEUR_TMA_G")             # set of raw core ids (str) in that array file
    d = d[d.core_id.astype(str).isin(keep)]          # cells of final cores only

    Array names follow the cell tables: "tma" (TMA A), "COEUR_TMA_G" and "COEUR_TMA_H" (TMA B blocks G and H);
    the "sfe_" prefix is accepted.
"""
import os

import pandas as pd

CORES_CSV = os.environ.get("FINAL_TMA_CORES", os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "07_core_qc/final_tma_cores.csv"))
PATIENTS_CSV = os.path.join(os.path.dirname(CORES_CSV), "final_tma_patients.csv")
ARRAY_BLOCK = {"tma": ("A", "A"), "COEUR_TMA_G": ("B", "G"), "COEUR_TMA_H": ("B", "H")}


def final_cores():
    """Per-core table (every profiled tumour core, with `final` and the reason for any exclusion)."""
    return pd.read_csv(CORES_CSV, dtype={"core_id_raw": str, "patient_id_raw": str})


def final_patients():
    """Per-patient table (`final`, `n_final_cores`, `two_final_cores`, cells at each step)."""
    return pd.read_csv(PATIENTS_CSV, dtype={"patient_id_raw": str})


def final_core_ids(array):
    """Raw core ids of the final cores in one array file ("tma", "COEUR_TMA_G", "COEUR_TMA_H")."""
    tma, block = ARRAY_BLOCK[array.removeprefix("sfe_")]
    c = final_cores()
    return set(c.loc[c.final & (c.tma == tma) & (c.block == block), "core_id_raw"])
