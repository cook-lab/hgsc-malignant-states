#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/22d_signature_survival.py, the script that ran.
# Kept: the 177 + 177 gene mean z-scores, the BayesPrism epithelial fraction, the clinical join and the polarization groups,
# written to 22d_signature_scores.csv (read by Fig 7g-j). Removed: the KM and Cox analyses with their figures, tables and HTML
# report, and the score plots (Fig 7g-j fits its own models). The untrimmed script is in the lab project folder.
"""
Step 22d — Signature-based survival with epithelial fraction adjustment
========================================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Instead of trying to deconvolve SecA vs SecB as separate cell types
(which fails because their expression profiles are too similar),
this script:
  1. Uses BayesPrism epithelial fraction (sum of all epithelial subtypes)
     to quantify TME contamination per sample
  2. Scores each TCGA-OV sample with SecA and SecB gene signatures
     (177 genes each, from NMF Step 11e) using mean z-score
  3. Tests survival (OS + PFS, full + 5-year) using:
     - Signature scores alone (KM + Cox)
     - Signature scores adjusted for epithelial purity (Cox)
     - SecA-SecB polarization score (difference)

Input:
    - output/22_tcga_deconvolution/CIBERSORTx_Results.txt  (BayesPrism fractions)
    - output/11e_nmf_characterization/11e_gene_classification.csv
    - data/cibersort_data_prev/tcga_ecotyper.txt
    - data/cibersort_data_prev/tcga_hla_clinical.csv

Output (output/22_tcga_deconvolution/):
    - 22d_signature_survival_report.html
    - 22d_km_*.svg/pdf
    - 22d_cox_forest_*.svg/pdf
    - 22d_signature_scores.csv
    - 22d_cox_results.csv
    - 22d_km_results.csv

Usage:
    python scripts/22d_signature_survival.py
"""

import os, sys, warnings, base64

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
PROJECT   = "/Volumes/CookLab/Sarah/2026_final_atlas"
OUT_DIR   = os.path.join(PROJECT, "output", "22_tcga_deconvolution")
FRACTIONS = os.path.join(OUT_DIR, "CIBERSORTx_Results.txt")
GENE_CLS  = os.path.join(PROJECT, "output", "11e_nmf_characterization",
                          "11e_gene_classification.csv")
TCGA_EXPR = os.path.join(PROJECT, "data", "cibersort_data_prev",
                          "tcga_ecotyper.txt")
CLINICAL  = os.path.join(PROJECT, "data", "cibersort_data_prev",
                          "tcga_hla_clinical.csv")

# ============================================================================
# HELPERS
# ============================================================================

def truncate_id(bc):
    return "-".join(str(bc).replace(".","-").split("-")[:3])


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 65)
    print("  Step 22d — Signature survival with epithelial adjustment")
    print("=" * 65)

    # ================================================================
    # 1. LOAD GENE SIGNATURES
    # ================================================================
    print("\n[1] Loading gene signatures...")
    gc = pd.read_csv(GENE_CLS)
    secA_genes = gc[gc["class"] == "SecA-specific"]["gene"].tolist()
    secB_genes = gc[gc["class"] == "SecB-specific"]["gene"].tolist()
    # Also get NMF loadings for weighted scoring
    secA_weights = dict(zip(gc[gc["class"]=="SecA-specific"]["gene"],
                            gc[gc["class"]=="SecA-specific"]["secA_loading"]))
    secB_weights = dict(zip(gc[gc["class"]=="SecB-specific"]["gene"],
                            gc[gc["class"]=="SecB-specific"]["secB_loading"]))
    print(f"    SecA genes: {len(secA_genes)}")
    print(f"    SecB genes: {len(secB_genes)}")

    # ================================================================
    # 2. LOAD TCGA EXPRESSION
    # ================================================================
    print("\n[2] Loading TCGA expression...")
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    print(f"    Expression: {expr.shape[0]} genes x {expr.shape[1]} samples")

    # Check gene overlap
    secA_in = [g for g in secA_genes if g in expr.index]
    secB_in = [g for g in secB_genes if g in expr.index]
    print(f"    SecA genes in TCGA: {len(secA_in)}/{len(secA_genes)}")
    print(f"    SecB genes in TCGA: {len(secB_in)}/{len(secB_genes)}")

    # ================================================================
    # 3. SCORE SAMPLES
    # ================================================================
    print("\n[3] Computing signature scores...")

    # Z-score normalize each gene across samples
    expr_z = expr.apply(lambda row: (row - row.mean()) / (row.std() + 1e-10), axis=1)

    # Mean z-score of gene sets
    secA_score = expr_z.loc[secA_in].mean(axis=0)
    secB_score = expr_z.loc[secB_in].mean(axis=0)

    # Weighted score (NMF loading-weighted mean z-score)
    secA_w = np.array([secA_weights[g] for g in secA_in])
    secB_w = np.array([secB_weights[g] for g in secB_in])
    secA_score_w = (expr_z.loc[secA_in].T * secA_w).T.sum(axis=0) / secA_w.sum()
    secB_score_w = (expr_z.loc[secB_in].T * secB_w).T.sum(axis=0) / secB_w.sum()

    # Polarization: SecA - SecB (positive = more SecA-like)
    polarization = secA_score - secB_score
    polarization_w = secA_score_w - secB_score_w

    scores = pd.DataFrame({
        "secA_score": secA_score,
        "secB_score": secB_score,
        "secA_score_weighted": secA_score_w,
        "secB_score_weighted": secB_score_w,
        "polarization": polarization,
        "polarization_weighted": polarization_w,
    })

    print(f"    SecA score: mean={secA_score.mean():.4f}, sd={secA_score.std():.4f}")
    print(f"    SecB score: mean={secB_score.mean():.4f}, sd={secB_score.std():.4f}")
    print(f"    Polarization (SecA-SecB): mean={polarization.mean():.4f}, sd={polarization.std():.4f}")
    print(f"    Correlation(SecA, SecB): r={np.corrcoef(secA_score, secB_score)[0,1]:.4f}")

    # ================================================================
    # 4. LOAD FRACTIONS & CLINICAL, JOIN
    # ================================================================
    print("\n[4] Loading fractions & clinical data...")

    frac = pd.read_csv(FRACTIONS, sep="\t", index_col=0)
    qc_cols = ["P-value", "Correlation", "RMSE"]
    frac = frac[[c for c in frac.columns if c not in qc_cols]]

    # Sum epithelial types
    epi_cols = [c for c in frac.columns if "epithelium" in c.lower()]
    frac["epi_fraction"] = frac[epi_cols].sum(axis=1)
    print(f"    Epithelial fraction: mean={frac['epi_fraction'].mean():.3f}, "
          f"median={frac['epi_fraction'].median():.3f}")

    clinical = pd.read_csv(CLINICAL, index_col=0)

    # Match IDs via truncation
    frac_t = {truncate_id(s): s for s in frac.index}
    clin_t = {truncate_id(s): s for s in clinical["PATIENT.NUMBER"]}
    score_t = {truncate_id(s): s for s in scores.index}

    common = set(frac_t) & set(clin_t) & set(score_t)
    print(f"    Three-way match: {len(common)} samples")

    # Build joined dataframe
    rows = []
    for tid in common:
        r = {}
        r["sample_id"] = tid

        # Scores
        sid_s = score_t[tid]
        for col in scores.columns:
            r[col] = scores.loc[sid_s, col]

        # Fractions
        sid_f = frac_t[tid]
        r["epi_fraction"] = frac.loc[sid_f, "epi_fraction"]

        # Clinical
        sid_c = clin_t[tid]
        crow = clinical[clinical["PATIENT.NUMBER"] == sid_c].iloc[0]
        r["os_months"] = pd.to_numeric(crow.get("OS..Months."), errors="coerce")
        r["os_event"] = 1 if "DECEASED" in str(crow.get("OS", "")) else 0
        r["pfs_months"] = pd.to_numeric(crow.get("PFS..Months."), errors="coerce")
        r["pfs_event"] = 1 if "PROGRESSION" in str(crow.get("PFS", "")) else 0
        r["stage_coded"] = pd.to_numeric(crow.get("Stage..Coded."), errors="coerce")
        r["age"] = pd.to_numeric(crow.get("Age"), errors="coerce")
        r["platinum_coded"] = pd.to_numeric(crow.get("Platinum..Coded."), errors="coerce")
        rows.append(r)

    df = pd.DataFrame(rows).set_index("sample_id")
    df = df[~df.index.duplicated(keep="first")]

    n_os = df[["os_months","os_event"]].dropna().shape[0]
    e_os = int(df.dropna(subset=["os_months","os_event"])["os_event"].sum())
    n_pfs = df[["pfs_months","pfs_event"]].dropna().shape[0]
    e_pfs = int(df.dropna(subset=["pfs_months","pfs_event"])["pfs_event"].sum())
    print(f"    OS: {n_os} patients ({e_os} events)")
    print(f"    PFS: {n_pfs} patients ({e_pfs} events)")

    # Create grouping variables (median splits)
    for col, label_lo, label_hi in [
        ("secA_score", "SecA-low", "SecA-high"),
        ("secB_score", "SecB-low", "SecB-high"),
        ("polarization", "SecB-like", "SecA-like"),
        ("polarization_weighted", "SecB-like (w)", "SecA-like (w)"),
    ]:
        med = df[col].median()
        df[f"{col}_group"] = pd.cut(df[col], bins=[-np.inf, med, np.inf],
                                     labels=[label_lo, label_hi])

    # Tertile for polarization
    t1, t2 = df["polarization"].quantile([1/3, 2/3])
    df["polar_tertile"] = pd.cut(df["polarization"],
                                  bins=[-np.inf, t1, t2, np.inf],
                                  labels=["SecB-like", "Intermediate", "SecA-like"])

    # Save scores
    df.to_csv(os.path.join(OUT_DIR, "22d_signature_scores.csv"))
    print(f"    Saved: 22d_signature_scores.csv")


if __name__ == "__main__":
    main()
