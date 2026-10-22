#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/11e_nmf_characterization.py, the script that ran.
# Kept: the 90th-percentile split of the Factor 3 (SecA) and Factor 2 (SecB) loadings, which writes 11e_gene_classification.csv
# (Tables S5-S6, TCGA scores, Figure S10). Removed: figures A-L, the enrichment, PROGENy, TF, cell-cycle and clinical
# analyses with their CSVs, and the HTML report. The untrimmed script is in the lab project folder.
"""
Step 11e — NMF Factor Characterization: SecA vs SecB Programs
=============================================================
Comprehensive biological characterization of the two epithelial
transcriptional programs identified by NMF in Step 11d.

Reads pre-computed NMF CSVs + h5ad metadata. Expression loaded in
backed mode only for PROGENy (~150 genes) and cell-cycle (~90 genes).

Produces 12 figures (A–L) as SVG+PDF, supporting CSVs, and an HTML report.

Usage
-----
    python scripts/11e_nmf_characterization.py
    python scripts/11e_nmf_characterization.py --skip-enrichment --skip-tf
"""

import argparse
import os
import sys
import time
import warnings

import h5py
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ========================================================================
# PATHS
# ========================================================================
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR   = os.path.join(BASE, "output", "11e_nmf_characterization")
USAGE_CSV = os.path.join(BASE, "output", "11d_epithelial_nmf", "11d_nmf_usage.csv")
LOAD_CSV  = os.path.join(BASE, "output", "11d_epithelial_nmf", "11d_nmf_loadings.csv")

os.makedirs(OUT_DIR, exist_ok=True)

# ========================================================================
# CONSTANTS
# ========================================================================
SECA_FACTOR = "Factor_3"          # progenitor-like
SECB_FACTOR = "Factor_2"          # adaptive / differentiated


# ========================================================================
# HELPERS
# ========================================================================
def save_csv(df, name, out_dir=OUT_DIR):
    """Save DataFrame as CSV."""
    path = os.path.join(out_dir, f"{name}.csv")
    df.to_csv(path)
    print(f"    Saved: {name}.csv")


# ========================================================================
# CLI
# ========================================================================
def parse_args():
    p = argparse.ArgumentParser(description="11e — NMF factor characterization")
    p.add_argument("--skip-enrichment", action="store_true",
                   help="Skip GO/Hallmark enrichment (no gseapy needed)")
    p.add_argument("--skip-tf", action="store_true",
                   help="Skip TF activity analysis (no decoupler needed)")
    return p.parse_args()


# ========================================================================
# MAIN
# ========================================================================
def main():
    args = parse_args()
    t0 = time.time()
    print("=" * 70)
    print("Step 11e — NMF Factor Characterization: SecA vs SecB Programs")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 0. LOAD DATA
    # ------------------------------------------------------------------
    print("\n[0] Loading data ...")
    loadings = pd.read_csv(LOAD_CSV, index_col=0)
    usage    = pd.read_csv(USAGE_CSV, index_col=0)
    hvg_genes = list(loadings.columns)
    n_genes   = len(hvg_genes)
    n_cells   = len(usage)
    print(f"    Loadings: {loadings.shape[0]} factors × {n_genes} genes")
    print(f"    Usage: {n_cells} cells × {loadings.shape[0]} factors")

    # ==================================================================
    # PART 1: TRANSCRIPTIONAL PROGRAM DEFINITION
    # ==================================================================
    print("\n" + "=" * 70)
    print("PART 1: Transcriptional Program Definition")
    print("=" * 70)

    # --- A: Loading scatter (F3 vs F2) --------------------------------
    print("\n[A] Loading scatter: Factor 3 (SecA) vs Factor 2 (SecB) ...")
    f2_load = loadings.loc[SECB_FACTOR]
    f3_load = loadings.loc[SECA_FACTOR]

    thr_f2 = np.percentile(f2_load.values, 90)
    thr_f3 = np.percentile(f3_load.values, 90)

    # Classify genes into quadrants
    is_secA = (f3_load > thr_f3) & (f2_load <= thr_f2)
    is_secB = (f2_load > thr_f2) & (f3_load <= thr_f3)
    is_shared = (f2_load > thr_f2) & (f3_load > thr_f3)

    # Save gene classification CSV
    gene_class = pd.DataFrame({
        "gene": hvg_genes,
        "secB_loading": f2_load.values,
        "secA_loading": f3_load.values,
        "class": np.where(is_secA, "SecA-specific",
                 np.where(is_secB, "SecB-specific",
                 np.where(is_shared, "Shared", "Background"))),
    }).sort_values("secA_loading", ascending=False)
    save_csv(gene_class, "11e_gene_classification")

    # ==================================================================
    # DONE
    # ==================================================================
    elapsed = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"Step 11e complete in {elapsed / 60:.1f} min")
    print(f"Output: {OUT_DIR}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
