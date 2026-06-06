#!/usr/bin/env python3
"""
CNV — statistical independence of subclone x SecA/SecB epitype
==============================================================
HGSC malignant-states atlas backend.

Tests whether SecA/SecB epitype is independent of CopyKAT subclone within each
sample. If SecB is transcriptional (not clonal), subclone identity should not
predict epitype. Per sample: chi-square / Fisher's exact on subclone x 3-bin
epitype, multinomial logistic (5-fold CV AUROC), and Cramer's V effect size.

INPUTS (output_root/05_cnv/):
  - per_sample/*/copykat_subclones.csv, per_sample/*/barcodes.csv
  - tables/sample_manifest.csv

OUTPUTS (output_root/05_cnv/):
  - tables/chisq_fisher_results.csv, logistic_regression_results.csv,
    summary_statistics.csv
  - figs/19_cnv_{pvalue_histogram,auroc_violin,cramers_v}.svg/pdf

MANUSCRIPT PANELS: supporting statistics for the CNV-independence claim
  (Fig 1J / SF7 narrative).

RUNTIME TIER: moderate.

SEEDING: LogisticRegression(random_state=SEED); AUROC strip jitter uses
  np.random.default_rng(SEED) for determinism.

Usage:
    python 04_cnv_independence.py
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import OneHotEncoder
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ── Cook Lab style v1.2 ─────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         8,
    "axes.titlesize":    9,
    "axes.labelsize":    8,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "legend.fontsize":   6,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# ── Paths ────────────────────────────────────────────────────
CNV_DIR    = path("output_root", "05_cnv")
PER_SAMPLE = os.path.join(CNV_DIR, "per_sample")
TBL_DIR    = os.path.join(CNV_DIR, "tables")
FIG_DIR    = os.path.join(CNV_DIR, "figs")
os.makedirs(TBL_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

MIN_CELLS_PER_GROUP = 20
MIN_SUBCLONES = 2
EPITYPES_3BIN = ["SecA", "Intermediate", "SecB"]


def cramers_v(contingency_table):
    chi2 = chi2_contingency(contingency_table)[0]
    n = contingency_table.values.sum()
    min_dim = min(contingency_table.shape) - 1
    if min_dim == 0 or n == 0:
        return 0.0
    return np.sqrt(chi2 / (n * min_dim))


def run_chisq_fisher(merged, sample_id):
    dat = merged[merged["epitype"].isin(EPITYPES_3BIN)].copy()
    if len(dat) < MIN_CELLS_PER_GROUP * 2:
        return None
    ct = pd.crosstab(dat["subclone"], dat["epitype"]).reindex(
        columns=EPITYPES_3BIN, fill_value=0)
    ct = ct[ct.sum(axis=1) >= MIN_CELLS_PER_GROUP]
    if ct.shape[0] < MIN_SUBCLONES:
        return None
    ct = ct.loc[:, ct.sum(axis=0) > 0]
    if ct.shape[1] < 2:
        return None
    n_total = ct.values.sum()
    if ct.shape == (2, 2):
        odds_ratio, p_value = fisher_exact(ct.values)
        test_used = "fisher"
    else:
        _chi2, p_value, _dof, _expected = chi2_contingency(ct.values)
        odds_ratio = np.nan
        test_used = "chi2"
    return {
        "sample_id": sample_id, "test": test_used, "p_value": p_value,
        "cramers_v": cramers_v(ct), "odds_ratio": odds_ratio,
        "n_subclones": ct.shape[0], "n_cells": n_total,
        "n_secA": int(ct["SecA"].sum()) if "SecA" in ct.columns else 0,
        "n_intermediate": int(ct["Intermediate"].sum()) if "Intermediate" in ct.columns else 0,
        "n_secB": int(ct["SecB"].sum()) if "SecB" in ct.columns else 0,
    }


def run_logistic(merged, sample_id):
    dat = merged[merged["epitype"].isin(EPITYPES_3BIN)].copy()
    if len(dat) < MIN_CELLS_PER_GROUP * 2:
        return None
    subclone_counts = dat["subclone"].value_counts()
    valid_clones = subclone_counts[subclone_counts >= MIN_CELLS_PER_GROUP].index
    if len(valid_clones) < MIN_SUBCLONES:
        return None
    dat = dat[dat["subclone"].isin(valid_clones)]
    y = dat["epitype"].values
    classes_present = [c for c in EPITYPES_3BIN if (y == c).sum() >= MIN_CELLS_PER_GROUP]
    if len(classes_present) < 2:
        return None
    dat = dat[dat["epitype"].isin(classes_present)]
    y = dat["epitype"].values

    enc = OneHotEncoder(sparse_output=False, drop="first")
    X = enc.fit_transform(dat[["subclone"]])
    if X.shape[1] == 0:
        return None
    n_folds = min(5, int(pd.Series(y).value_counts().min()))
    if n_folds < 2:
        return None
    clf = LogisticRegression(max_iter=1000, solver="lbfgs", random_state=SEED)
    scoring = "roc_auc" if len(classes_present) == 2 else "roc_auc_ovr"
    try:
        aurocs = cross_val_score(clf, X, y, cv=n_folds, scoring=scoring)
    except ValueError:
        return None
    return {
        "sample_id": sample_id, "auroc_mean": aurocs.mean(),
        "auroc_std": aurocs.std(), "scoring": scoring, "n_folds": n_folds,
        "n_cells": len(dat), "n_classes": len(classes_present),
        "n_subclones": len(valid_clones), "secB_fraction": float((y == "SecB").mean()),
    }


def main():
    print("=" * 65)
    print("  CNV - Epitype (3-bin) x CNV Subclone Independence Tests")
    print("=" * 65)

    completed = [s for s in sorted(os.listdir(PER_SAMPLE))
                 if os.path.isdir(os.path.join(PER_SAMPLE, s))
                 and os.path.exists(os.path.join(PER_SAMPLE, s, "DONE.txt"))]
    print(f"\n  Completed CopyKAT samples: {len(completed)}")

    chisq_results, logistic_results, skipped = [], [], 0
    for i, sample_id in enumerate(completed):
        sample_dir = os.path.join(PER_SAMPLE, sample_id)
        subclone_path = os.path.join(sample_dir, "copykat_subclones.csv")
        if not os.path.exists(subclone_path):
            skipped += 1
            continue
        subclones = pd.read_csv(subclone_path)
        barcodes = pd.read_csv(os.path.join(sample_dir, "barcodes.csv"))
        epi = barcodes[~barcodes["is_reference"]].copy()
        merged = epi.merge(subclones, on="barcode", how="inner")
        merged_aneuploid = merged[merged["subclone"] != "diploid"]
        if len(merged_aneuploid) < MIN_CELLS_PER_GROUP * 2:
            skipped += 1
            continue
        r1 = run_chisq_fisher(merged_aneuploid, sample_id)
        if r1:
            chisq_results.append(r1)
        r2 = run_logistic(merged_aneuploid, sample_id)
        if r2:
            logistic_results.append(r2)
        if (i + 1) % 10 == 0:
            print(f"  [{i+1}/{len(completed)}] processed")

    print(f"\n  Chi-square/Fisher tests: {len(chisq_results)} samples")
    print(f"  Logistic regression: {len(logistic_results)} samples")
    print(f"  Skipped (insufficient cells/subclones): {skipped}")

    chisq_df = pd.DataFrame(chisq_results) if chisq_results else None
    if chisq_results:
        from statsmodels.stats.multitest import multipletests
        _, chisq_df["p_adj"], _, _ = multipletests(chisq_df["p_value"], method="fdr_bh")
        chisq_df.to_csv(os.path.join(TBL_DIR, "chisq_fisher_results.csv"), index=False)
        print("\n  Saved chisq_fisher_results.csv")
        n_sig = (chisq_df["p_adj"] < 0.05).sum()
        print(f"  Significant (FDR < 0.05): {n_sig}/{len(chisq_df)} "
              f"({100*n_sig/len(chisq_df):.1f}%)")
        print(f"  Median Cramer's V: {chisq_df['cramers_v'].median():.3f}")

    logistic_df = pd.DataFrame(logistic_results) if logistic_results else None
    if logistic_results:
        logistic_df.to_csv(os.path.join(TBL_DIR, "logistic_regression_results.csv"),
                           index=False)
        print("  Saved logistic_regression_results.csv")
        print(f"  Median AUROC: {logistic_df['auroc_mean'].median():.3f}")

    summary = {
        "n_samples_completed": len(completed),
        "n_samples_tested_chisq": len(chisq_results),
        "n_samples_tested_logistic": len(logistic_results),
        "n_significant_fdr05": int((chisq_df["p_adj"] < 0.05).sum()) if chisq_results else 0,
        "median_cramers_v": chisq_df["cramers_v"].median() if chisq_results else np.nan,
        "median_auroc": logistic_df["auroc_mean"].median() if logistic_results else np.nan,
        "mean_auroc": logistic_df["auroc_mean"].mean() if logistic_results else np.nan,
    }
    pd.DataFrame([summary]).to_csv(os.path.join(TBL_DIR, "summary_statistics.csv"),
                                   index=False)

    print("\n  Generating figures...")
    if chisq_results:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.hist(chisq_df["p_value"], bins=20, color="#E6A141", alpha=0.8,
                edgecolor="#333333", linewidth=0.4)
        ax.axhline(len(chisq_df) / 20, color="red", ls="--", lw=1,
                   label="Uniform expectation")
        ax.set_xlabel("p-value (chi-square / Fisher)")
        ax.set_ylabel("Number of samples")
        ax.set_title("SecA/SecB x CNV Subclone Independence", fontweight="bold")
        ax.legend(fontsize=6)
        for fmt in ["svg", "pdf"]:
            fig.savefig(os.path.join(FIG_DIR, f"19_cnv_pvalue_histogram.{fmt}"),
                        bbox_inches="tight")
        plt.close()
        print("  Saved pvalue_histogram")

        fig, ax = plt.subplots(figsize=(4, 3))
        ax.hist(chisq_df["cramers_v"], bins=20, color="#B8741A", alpha=0.8,
                edgecolor="#333333", linewidth=0.4)
        ax.axvline(0.1, color="red", ls="--", lw=1, label="Small effect (0.1)")
        ax.set_xlabel("Cramer's V"); ax.set_ylabel("Number of samples")
        ax.set_title("Effect Size: Subclone x Epitype Association", fontweight="bold")
        ax.legend(fontsize=6)
        for fmt in ["svg", "pdf"]:
            fig.savefig(os.path.join(FIG_DIR, f"19_cnv_cramers_v.{fmt}"),
                        bbox_inches="tight")
        plt.close()
        print("  Saved cramers_v")

    if logistic_results:
        fig, ax = plt.subplots(figsize=(3, 4))
        parts = ax.violinplot([logistic_df["auroc_mean"].values],
                              showmedians=True, showextrema=False)
        parts["bodies"][0].set_facecolor("#E6A141")
        parts["bodies"][0].set_alpha(0.7)
        parts["cmedians"].set_color("black")
        ax.axhline(0.5, color="red", ls="--", lw=1, label="Chance (0.5)")
        ax.set_ylabel("AUROC (5-fold CV)")
        ax.set_xticks([1]); ax.set_xticklabels(["Subclone -> SecB"])
        ax.set_title("Logistic Regression:\nCNV Subclone Predicting Epitype",
                     fontweight="bold")
        ax.legend(fontsize=6)
        jitter = np.random.default_rng(SEED).uniform(-0.05, 0.05, len(logistic_df))
        ax.scatter(1 + jitter, logistic_df["auroc_mean"], s=8, c="#7D4E4E",
                   alpha=0.6, edgecolors="none", zorder=5)
        for fmt in ["svg", "pdf"]:
            fig.savefig(os.path.join(FIG_DIR, f"19_cnv_auroc_violin.{fmt}"),
                        bbox_inches="tight")
        plt.close()
        print("  Saved auroc_violin")

    print(f"\n{'=' * 65}\n  DONE\n{'=' * 65}")


if __name__ == "__main__":
    main()
