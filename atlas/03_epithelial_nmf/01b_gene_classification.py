#!/usr/bin/env python3
"""
01b_gene_classification.py — per-gene SecA/SecB program classification (11e)
============================================================================
PURPOSE
    Classify every HVG into the SecA-specific / SecB-specific / Shared /
    Background program from the NMF gene loadings, and record each gene's SecA
    (Factor_3) and SecB (Factor_2) loading. A gene is "SecA-specific" if its
    SecA loading is in the top decile AND its SecB loading is not; "SecB-specific"
    is the mirror; "Shared" is top-decile on both; everything else is "Background".
    (Faithful port of the gene-classification step of the source
    11e_nmf_characterization.py — the per-gene loading split, not the figures.)

INPUTS
    - output_root/03_epithelial_nmf/11d_nmf_loadings.csv   (factors x genes; the
      H matrix written by 01_epithelial_nmf.py)

OUTPUTS
    - output_root/03_epithelial_nmf/11e_gene_classification.csv
      (columns: gene, secB_loading, secA_loading, class)

MANUSCRIPT PANEL(S): upstream of Fig 7E/F/G — the signature-survival cache
    (atlas/07_deconvolution_survival/05_signature_survival.py + 06_validate_survival.py)
    reads this to build the SecA/SecB survival signatures.

RUNTIME TIER: fast

NOTE: SecA = Factor_3, SecB = Factor_2 (config polarization.factor = Factor_2;
    the seeded NMF is deterministic, so the factor->program mapping is stable —
    see docs/REPRODUCIBILITY.md). Migrated to close the missing-producer gap
    (audit A4/H4): 05/06 previously read this table from a bare data_root literal
    into the untouched source tree, with no in-repo generator.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

SECB_FACTOR = "Factor_2"   # SecB-defining factor (config polarization.factor)
SECA_FACTOR = "Factor_3"   # SecA-defining factor

LOAD_CSV = path("output_root", "03_epithelial_nmf", "11d_nmf_loadings.csv")
OUT_CSV  = path("output_root", "03_epithelial_nmf", "11e_gene_classification.csv")


def main():
    loadings = pd.read_csv(LOAD_CSV, index_col=0)   # factors x genes (H matrix)
    hvg_genes = list(loadings.columns)

    f2_load = loadings.loc[SECB_FACTOR]   # SecB loading per gene
    f3_load = loadings.loc[SECA_FACTOR]   # SecA loading per gene

    thr_f2 = np.percentile(f2_load.values, 90)
    thr_f3 = np.percentile(f3_load.values, 90)

    is_secA   = (f3_load > thr_f3) & (f2_load <= thr_f2)
    is_secB   = (f2_load > thr_f2) & (f3_load <= thr_f3)
    is_shared = (f2_load > thr_f2) & (f3_load > thr_f3)

    gene_class = pd.DataFrame({
        "gene": hvg_genes,
        "secB_loading": f2_load.values,
        "secA_loading": f3_load.values,
        "class": np.where(is_secA, "SecA-specific",
                 np.where(is_secB, "SecB-specific",
                 np.where(is_shared, "Shared", "Background"))),
    }).sort_values("secA_loading", ascending=False)

    gene_class.to_csv(OUT_CSV, index=False)
    counts = gene_class["class"].value_counts()
    print(f"Wrote {OUT_CSV}")
    print(counts.to_string())


if __name__ == "__main__":
    main()
