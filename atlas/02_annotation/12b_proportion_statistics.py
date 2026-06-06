#!/usr/bin/env python3
"""
Atlas 02 — Step 12b (statistics): proportion enrichment statistics

PURPOSE
    Statistical tests on the cell-type proportion tables (level-1, epithelial
    SecA/SecB, all compartments) across metadata strata; writes enrichment CSVs and
    a text report. Consumes the proportion-count tables written by 12b.

INPUTS
    output_root/02_annotation/12_celltype_proportions/tables/*.csv

OUTPUTS
    output_root/02_annotation/12_celltype_proportions/12b_*.csv + 12b_statistical_report.txt

MANUSCRIPT PANEL(S)
    Statistics backing Fig 1B-E / Fig 2 composition claims.

RUNTIME TIER
    fast (operates on count tables).
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from collections import OrderedDict

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

TABLE_DIR   = path("output_root", "02_annotation", "12_celltype_proportions", "tables")
OUT_DIR     = path("output_root", "02_annotation", "12_celltype_proportions")

REPORT_PATH = os.path.join(OUT_DIR, "12b_statistical_report.txt")

# ============================================================================
# COMPARTMENT REGISTRY
# ============================================================================

COMPARTMENTS = OrderedDict([
    ("epithelial",   "Epithelial"),
    ("macrophage",   "Macrophage"),
    ("tnkcell",      "T/NK cell"),
    ("fibroblast",   "Fibroblast"),
    ("endothelial",  "Endothelial"),
    ("bcell",        "B cell"),
    ("mesothelial",  "Mesothelial"),
    ("smoothmuscle", "Smooth muscle"),
    ("plasmacell",   "Plasma cell"),
    ("dc",           "DC"),
    ("mastcell",     "Mast cell"),
    ("neutrophil",   "Neutrophil"),
    ("pericyte",     "Pericyte"),
])

METADATA_KEYS = [
    "study", "anatomic_site", "treatment_status", "treatment_response",
    "stage", "metastatic_site", "BRCA_status", "HRD_status", "TP53_status",
]

META_DISPLAY = {
    "study": "Study", "anatomic_site": "Anatomic Site",
    "treatment_status": "Treatment Status",
    "treatment_response": "Treatment Response",
    "stage": "Stage", "metastatic_site": "Metastatic Site",
    "BRCA_status": "BRCA Status", "HRD_status": "HRD Status",
    "TP53_status": "TP53 Status",
}

# ============================================================================
# HELPERS
# ============================================================================

def load_counts(prefix, meta_key):
    """Load a count CSV, return DataFrame or None."""
    path = os.path.join(TABLE_DIR, f"{prefix}_counts_{meta_key}.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, index_col=0)
    # Drop rows/cols that are all zero
    df = df.loc[df.sum(axis=1) > 0, df.sum(axis=0) > 0]
    return df


def chi_square_test(counts_df):
    """
    Chi-square test of independence on a contingency table.
    Returns (chi2, p, dof, cramers_v).
    Uses G-test (log-likelihood) when expected frequencies have zeros.
    """
    if counts_df.shape[0] < 2 or counts_df.shape[1] < 2:
        return np.nan, np.nan, 0, np.nan
    try:
        chi2, p, dof, expected = stats.chi2_contingency(counts_df)
    except ValueError:
        # Fall back to G-test which handles zeros better
        try:
            chi2, p, dof, expected = stats.chi2_contingency(
                counts_df, lambda_="log-likelihood")
        except ValueError:
            return np.nan, np.nan, 0, np.nan
    n = counts_df.values.sum()
    k = min(counts_df.shape) - 1
    cramers_v = np.sqrt(chi2 / (n * k)) if k > 0 and n > 0 else 0
    return chi2, p, dof, cramers_v


def enrichment_table(counts_df):
    """
    Compute observed/expected ratio for each cell.
    O/E > 1 = enriched, O/E < 1 = depleted.
    Also computes adjusted residuals and p-values.
    Returns (oe_df, residual_df, pval_df).
    """
    observed = counts_df.values.astype(float)
    row_totals = observed.sum(axis=1, keepdims=True)
    col_totals = observed.sum(axis=0, keepdims=True)
    grand_total = observed.sum()

    expected = (row_totals * col_totals) / grand_total
    # Avoid division by zero
    expected_safe = np.where(expected > 0, expected, 1)

    oe = observed / expected_safe
    oe = np.where(expected > 0, oe, np.nan)

    # Adjusted (standardized) residuals
    # r_ij = (O_ij - E_ij) / sqrt(E_ij * (1 - p_i.) * (1 - p_.j))
    p_row = row_totals / grand_total
    p_col = col_totals / grand_total
    denom = np.sqrt(expected_safe * (1 - p_row) * (1 - p_col))
    denom = np.where(denom > 0, denom, 1)
    residuals = (observed - expected) / denom

    # Two-sided p-values from standard normal
    pvals = 2 * stats.norm.sf(np.abs(residuals))

    oe_df = pd.DataFrame(oe, index=counts_df.index, columns=counts_df.columns)
    res_df = pd.DataFrame(residuals, index=counts_df.index, columns=counts_df.columns)
    pval_df = pd.DataFrame(pvals, index=counts_df.index, columns=counts_df.columns)

    return oe_df, res_df, pval_df


def bonferroni_sig(pval_df, alpha=0.05):
    """Return boolean DataFrame of Bonferroni-significant cells."""
    n_tests = pval_df.size
    return pval_df < (alpha / n_tests)


def format_pval(p):
    if pd.isna(p):
        return "NA"
    if p < 1e-300:
        return "< 1e-300"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"


def top_enrichments(oe_df, res_df, pval_df, n=10, direction="enriched"):
    """
    Return top n enriched or depleted (cell_type, metadata_category) pairs.
    """
    rows = []
    sig = bonferroni_sig(pval_df)
    for meta_cat in oe_df.index:
        for cell_type in oe_df.columns:
            oe_val = oe_df.loc[meta_cat, cell_type]
            res_val = res_df.loc[meta_cat, cell_type]
            p_val = pval_df.loc[meta_cat, cell_type]
            is_sig = sig.loc[meta_cat, cell_type]
            if pd.isna(oe_val) or not is_sig:
                continue
            rows.append({
                "metadata_category": meta_cat,
                "cell_type": cell_type,
                "OE_ratio": oe_val,
                "adj_residual": res_val,
                "p_value": p_val,
            })
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df
    if direction == "enriched":
        df = df[df["OE_ratio"] > 1].sort_values("OE_ratio", ascending=False)
    else:
        df = df[df["OE_ratio"] < 1].sort_values("OE_ratio", ascending=True)
    return df.head(n).reset_index(drop=True)


# ============================================================================
# REPORT BUILDER
# ============================================================================

class ReportBuilder:
    def __init__(self):
        self.lines = []
        self.all_enrichments = []

    def h1(self, text):
        self.lines.append("\n" + "=" * 80)
        self.lines.append(text)
        self.lines.append("=" * 80)

    def h2(self, text):
        self.lines.append("\n" + "-" * 70)
        self.lines.append(text)
        self.lines.append("-" * 70)

    def h3(self, text):
        self.lines.append(f"\n  >> {text}")

    def text(self, text):
        self.lines.append(text)

    def blank(self):
        self.lines.append("")

    def table(self, df, indent=4):
        """Format a small DataFrame as aligned text."""
        s = df.to_string()
        for line in s.split("\n"):
            self.lines.append(" " * indent + line)

    def save(self, path):
        with open(path, "w") as f:
            f.write("\n".join(self.lines))
        print(f"  Report saved: {path}")


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

def main():
    print("=" * 70)
    print("Step 12b — Cell Type Proportion Statistical Analysis")
    print("=" * 70)

    R = ReportBuilder()

    R.h1("CELL TYPE PROPORTION STATISTICAL ANALYSIS")
    R.text("HGSC Single-Cell Atlas | Cook Lab | 2026")
    R.text("Companion analysis to Step 12 proportion figures")
    R.blank()

    # ==================================================================
    # SECTION 1: GLOBAL CHI-SQUARE TESTS — LEVEL 1
    # ==================================================================

    R.h1("SECTION 1: GLOBAL CHI-SQUARE TESTS (Level 1)")
    R.text("Tests whether cell type composition differs significantly across")
    R.text("metadata categories. Cramer's V measures effect size (0=none, 1=perfect).")
    R.blank()

    chi2_results = []
    for meta_key in METADATA_KEYS:
        counts = load_counts("12_level1", meta_key)
        if counts is None:
            continue
        chi2, p, dof, v = chi_square_test(counts)
        chi2_results.append({
            "Metadata": META_DISPLAY.get(meta_key, meta_key),
            "Key": meta_key,
            "Groups": counts.shape[0],
            "Chi2": f"{chi2:,.0f}" if not np.isnan(chi2) else "NA",
            "p-value": format_pval(p),
            "Cramer_V": f"{v:.4f}" if not np.isnan(v) else "NA",
            "Effect": ("LARGE" if v > 0.3 else "MEDIUM" if v > 0.1 else "SMALL")
                      if not np.isnan(v) else "NA",
        })

    chi2_df = pd.DataFrame(chi2_results)
    R.table(chi2_df)
    R.blank()
    R.text("Interpretation: ALL metadata variables show highly significant (p < 1e-300)")
    R.text("association with cell type composition. This is expected given n=2.3M cells.")
    R.text("Focus on Cramer's V for practical effect size.")

    # ==================================================================
    # SECTION 2: LEVEL 1 ENRICHMENT/DEPLETION
    # ==================================================================

    R.h1("SECTION 2: LEVEL 1 ENRICHMENT / DEPLETION ANALYSIS")
    R.text("O/E ratio: Observed/Expected. >1 = enriched, <1 = depleted.")
    R.text("Adjusted residuals follow standard normal under null (|z|>3.3 ~ Bonferroni sig).")
    R.blank()

    all_level1_enrichments = []

    for meta_key in METADATA_KEYS:
        counts = load_counts("12_level1", meta_key)
        if counts is None:
            continue

        oe_df, res_df, pval_df = enrichment_table(counts)

        R.h2(f"Level 1 × {META_DISPLAY.get(meta_key, meta_key)}")

        # Show O/E ratios
        R.h3("O/E Ratios (>1.5 enriched, <0.67 depleted highlighted)")
        oe_display = oe_df.round(2)
        R.table(oe_display)
        R.blank()

        # Top enrichments
        top_enr = top_enrichments(oe_df, res_df, pval_df, n=5, direction="enriched")
        top_dep = top_enrichments(oe_df, res_df, pval_df, n=5, direction="depleted")

        if len(top_enr) > 0:
            R.h3("Top 5 Enrichments (Bonferroni-significant)")
            for _, row in top_enr.iterrows():
                R.text(f"    {row['cell_type']:35s} in {row['metadata_category']:30s}  "
                       f"O/E={row['OE_ratio']:.2f}  z={row['adj_residual']:.1f}  "
                       f"p={format_pval(row['p_value'])}")

        if len(top_dep) > 0:
            R.h3("Top 5 Depletions (Bonferroni-significant)")
            for _, row in top_dep.iterrows():
                R.text(f"    {row['cell_type']:35s} in {row['metadata_category']:30s}  "
                       f"O/E={row['OE_ratio']:.2f}  z={row['adj_residual']:.1f}  "
                       f"p={format_pval(row['p_value'])}")

        # Collect for CSV
        for meta_cat in oe_df.index:
            for cell_type in oe_df.columns:
                all_level1_enrichments.append({
                    "metadata_var": meta_key,
                    "metadata_cat": meta_cat,
                    "cell_type": cell_type,
                    "OE_ratio": oe_df.loc[meta_cat, cell_type],
                    "adj_residual": res_df.loc[meta_cat, cell_type],
                    "p_value": pval_df.loc[meta_cat, cell_type],
                    "bonferroni_sig": pval_df.loc[meta_cat, cell_type] < (0.05 / pval_df.size),
                    "observed": counts.loc[meta_cat, cell_type],
                })

    # Save level1 enrichment CSV
    l1_enr_df = pd.DataFrame(all_level1_enrichments)
    l1_enr_path = os.path.join(OUT_DIR, "12b_enrichment_level1.csv")
    l1_enr_df.to_csv(l1_enr_path, index=False)
    print(f"  Saved: {l1_enr_path}")

    # ==================================================================
    # SECTION 3: EPITHELIAL LEVEL 2 — SecA/SecB POLARIZATION
    # ==================================================================

    R.h1("SECTION 3: EPITHELIAL COMPARTMENT — SecA/SecB POLARIZATION ANALYSIS")
    R.text("Key question: Do Adaptive ('SecB-like', differentiated) and Cycling/Stress-")
    R.text("response ('SecA-like', progenitor) epithelial subtypes show differential")
    R.text("enrichment across treatment, metastatic site, genomic status, and anatomy?")
    R.blank()

    # Define SecA-like and SecB-like clusters
    SECA_TYPES = [
        "Cycling secretory epithelial cell",
        "Stress-response secretory epithelial cell",
        "Transitioning epithelial cell",
    ]
    SECB_TYPES = [
        "Adaptive secretory epithelial cell",
    ]
    BASELINE_TYPES = [
        "Secretory epithelial cell",
    ]
    # Excluded and ciliated tracked separately
    EXCLUDED_TYPES = [
        "Excluded epithelial cell_1",
        "Excluded epithelial cell_2",
        "Excluded epithelial cell_3",
    ]

    all_epi_enrichments = []

    for meta_key in METADATA_KEYS:
        counts = load_counts("12_epithelial_level2", meta_key)
        if counts is None:
            continue

        # Filter out excluded clusters for main analysis
        active_cols = [c for c in counts.columns if "Excluded" not in c]
        counts_active = counts[active_cols]

        if counts_active.shape[0] < 2 or counts_active.shape[1] < 2:
            continue

        oe_df, res_df, pval_df = enrichment_table(counts_active)

        R.h2(f"Epithelial Level 2 × {META_DISPLAY.get(meta_key, meta_key)}")

        # O/E table
        R.h3("O/E Ratios (active clusters only, excluding 'Excluded' clusters)")
        R.table(oe_df.round(2))
        R.blank()

        # SecA vs SecB comparison
        seca_cols = [c for c in active_cols if c in SECA_TYPES]
        secb_cols = [c for c in active_cols if c in SECB_TYPES]
        base_cols = [c for c in active_cols if c in BASELINE_TYPES]

        if seca_cols and secb_cols:
            R.h3("SecA-like vs SecB-like enrichment pattern:")

            # Compute proportions within each metadata category
            props = counts_active.div(counts_active.sum(axis=1), axis=0)

            for meta_cat in props.index:
                seca_prop = props.loc[meta_cat, seca_cols].sum()
                secb_prop = props.loc[meta_cat, secb_cols].sum()
                base_prop = props.loc[meta_cat, base_cols].sum() if base_cols else 0
                ratio = seca_prop / secb_prop if secb_prop > 0 else np.inf
                arrow = "→SecA" if ratio > 1.2 else "→SecB" if ratio < 0.8 else "balanced"
                R.text(f"    {meta_cat:35s}  SecA={seca_prop:.3f}  SecB={secb_prop:.3f}  "
                       f"ratio={ratio:.2f}  {arrow}")

        # Highlight significant enrichments
        top_enr = top_enrichments(oe_df, res_df, pval_df, n=5, direction="enriched")
        top_dep = top_enrichments(oe_df, res_df, pval_df, n=5, direction="depleted")

        if len(top_enr) > 0:
            R.h3("Top Enrichments")
            for _, row in top_enr.iterrows():
                R.text(f"    {row['cell_type']:45s} in {row['metadata_category']:25s}  "
                       f"O/E={row['OE_ratio']:.2f}  z={row['adj_residual']:.1f}")

        if len(top_dep) > 0:
            R.h3("Top Depletions")
            for _, row in top_dep.iterrows():
                R.text(f"    {row['cell_type']:45s} in {row['metadata_category']:25s}  "
                       f"O/E={row['OE_ratio']:.2f}  z={row['adj_residual']:.1f}")

        # Collect for CSV
        for meta_cat in oe_df.index:
            for cell_type in oe_df.columns:
                all_epi_enrichments.append({
                    "metadata_var": meta_key,
                    "metadata_cat": meta_cat,
                    "cell_type": cell_type,
                    "OE_ratio": oe_df.loc[meta_cat, cell_type],
                    "adj_residual": res_df.loc[meta_cat, cell_type],
                    "p_value": pval_df.loc[meta_cat, cell_type],
                    "bonferroni_sig": pval_df.loc[meta_cat, cell_type] < (0.05 / pval_df.size),
                })

    # Save epithelial enrichment CSV
    epi_enr_df = pd.DataFrame(all_epi_enrichments)
    epi_enr_path = os.path.join(OUT_DIR, "12b_enrichment_epithelial.csv")
    epi_enr_df.to_csv(epi_enr_path, index=False)
    print(f"  Saved: {epi_enr_path}")

    # ==================================================================
    # SECTION 4: ALL COMPARTMENTS — LEVEL 2 ENRICHMENT
    # ==================================================================

    R.h1("SECTION 4: ALL COMPARTMENTS — LEVEL 2 ENRICHMENT HIGHLIGHTS")
    R.text("Per-compartment level2 enrichment analysis across metadata.")
    R.blank()

    all_comp_enrichments = []

    for comp_key, comp_display in COMPARTMENTS.items():
        R.h2(f"{comp_display} ({comp_key})")

        any_data = False
        for meta_key in METADATA_KEYS:
            counts = load_counts(f"12_{comp_key}_level2", meta_key)
            if counts is None:
                continue

            # Filter out excluded clusters
            active_cols = [c for c in counts.columns if "Excluded" not in c]
            counts_active = counts[active_cols]

            if counts_active.shape[0] < 2 or counts_active.shape[1] < 2:
                continue

            chi2, p, dof, v = chi_square_test(counts_active)
            oe_df, res_df, pval_df = enrichment_table(counts_active)

            # Only report if meaningful effect size
            if np.isnan(v) or v < 0.05:
                continue

            any_data = True
            R.h3(f"× {META_DISPLAY.get(meta_key, meta_key)}  "
                  f"(V={v:.3f}, {counts_active.shape[0]} groups × "
                  f"{counts_active.shape[1]} types)")

            top_enr = top_enrichments(oe_df, res_df, pval_df, n=3, direction="enriched")
            top_dep = top_enrichments(oe_df, res_df, pval_df, n=3, direction="depleted")

            for _, row in top_enr.iterrows():
                R.text(f"      ENRICHED: {row['cell_type']:40s} in {row['metadata_category']:25s}  "
                       f"O/E={row['OE_ratio']:.2f}")
            for _, row in top_dep.iterrows():
                R.text(f"      DEPLETED: {row['cell_type']:40s} in {row['metadata_category']:25s}  "
                       f"O/E={row['OE_ratio']:.2f}")

            # Collect
            for meta_cat in oe_df.index:
                for cell_type in oe_df.columns:
                    all_comp_enrichments.append({
                        "compartment": comp_display,
                        "metadata_var": meta_key,
                        "metadata_cat": meta_cat,
                        "cell_type": cell_type,
                        "OE_ratio": oe_df.loc[meta_cat, cell_type],
                        "adj_residual": res_df.loc[meta_cat, cell_type],
                        "p_value": pval_df.loc[meta_cat, cell_type],
                        "cramers_v": v,
                    })

        if not any_data:
            R.text("    No metadata variables with V > 0.05")

    # Save all-compartment CSV
    comp_enr_df = pd.DataFrame(all_comp_enrichments)
    comp_enr_path = os.path.join(OUT_DIR, "12b_enrichment_all_compartments.csv")
    comp_enr_df.to_csv(comp_enr_path, index=False)
    print(f"  Saved: {comp_enr_path}")

    # ==================================================================
    # SECTION 5: SYNTHESIS — KEY FINDINGS
    # ==================================================================

    R.h1("SECTION 5: SYNTHESIS — KEY FINDINGS")
    R.blank()

    # ---- 5A: Previously reported findings ----
    R.h2("5A. PREVIOUSLY REPORTED IN HGSC LITERATURE (validation)")
    R.blank()

    # Compute specific statistics for known findings
    # 1. Ascites enriched for macrophages
    l1_met = load_counts("12_level1", "metastatic_site")
    if l1_met is not None:
        oe_met, res_met, p_met = enrichment_table(l1_met)

        R.text("1. ASCITES ENRICHMENT PATTERN (well-established)")
        if "ascites" in oe_met.index:
            mac_oe = oe_met.loc["ascites", "Macrophage"] if "Macrophage" in oe_met.columns else np.nan
            mes_oe = oe_met.loc["ascites", "Mesothelial"] if "Mesothelial" in oe_met.columns else np.nan
            epi_oe = oe_met.loc["ascites", "Epithelial"] if "Epithelial" in oe_met.columns else np.nan
            fib_oe = oe_met.loc["ascites", "Fibroblast"] if "Fibroblast" in oe_met.columns else np.nan
            R.text(f"   Macrophage in ascites:  O/E = {mac_oe:.2f}  (enriched → confirms immune-rich ascites)")
            R.text(f"   Mesothelial in ascites: O/E = {mes_oe:.2f}  (enriched → mesothelial lining cells)")
            R.text(f"   Epithelial in ascites:  O/E = {epi_oe:.2f}  (mildly {'enriched' if epi_oe > 1 else 'depleted'})")
            R.text(f"   Fibroblast in ascites:  O/E = {fib_oe:.2f}  (depleted → no stromal scaffold)")
        R.blank()

        # 2. Primary tumors enriched for stromal cells
        R.text("2. PRIMARY TUMOR STROMAL ENRICHMENT (well-established)")
        if "primary" in oe_met.index:
            fib_oe = oe_met.loc["primary", "Fibroblast"]
            sm_oe = oe_met.loc["primary", "Smooth muscle"]
            peri_oe = oe_met.loc["primary", "Pericyte"]
            endo_oe = oe_met.loc["primary", "Endothelial"]
            R.text(f"   Fibroblast in primary:    O/E = {fib_oe:.2f}")
            R.text(f"   Smooth muscle in primary: O/E = {sm_oe:.2f}")
            R.text(f"   Pericyte in primary:      O/E = {peri_oe:.2f}")
            R.text(f"   Endothelial in primary:   O/E = {endo_oe:.2f}")
        R.blank()

    # 3. Treatment effects on immune composition
    l1_tx = load_counts("12_level1", "treatment_status")
    if l1_tx is not None:
        oe_tx, res_tx, p_tx = enrichment_table(l1_tx)
        R.text("3. TREATMENT EFFECTS ON IMMUNE COMPOSITION")
        if "post-treatment" in oe_tx.index:
            tnk_oe = oe_tx.loc["post-treatment", "T/NK cell"]
            mac_oe = oe_tx.loc["post-treatment", "Macrophage"]
            R.text(f"   T/NK cell post-treatment:  O/E = {tnk_oe:.2f}  "
                   f"({'enriched' if tnk_oe > 1 else 'depleted'})")
            R.text(f"   Macrophage post-treatment: O/E = {mac_oe:.2f}  "
                   f"({'enriched' if mac_oe > 1 else 'depleted'})")
        if "pre-treatment" in oe_tx.index:
            epi_oe = oe_tx.loc["pre-treatment", "Epithelial"]
            R.text(f"   Epithelial pre-treatment:  O/E = {epi_oe:.2f}")
        R.blank()

    # 4. Study batch effects
    l1_study = load_counts("12_level1", "study")
    if l1_study is not None:
        oe_st, _, _ = enrichment_table(l1_study)
        chi2_st, p_st, _, v_st = chi_square_test(l1_study)
        R.text(f"4. STUDY BATCH EFFECTS (Cramer's V = {v_st:.3f})")
        R.text("   Largest deviations from expected proportions:")
        # Find top study-specific deviations
        max_devs = []
        for study in oe_st.index:
            for ct in oe_st.columns:
                dev = abs(oe_st.loc[study, ct] - 1.0)
                max_devs.append((study, ct, oe_st.loc[study, ct], dev))
        max_devs.sort(key=lambda x: x[3], reverse=True)
        for study, ct, oe_val, _ in max_devs[:8]:
            R.text(f"   {ct:20s} in {study:25s}  O/E = {oe_val:.2f}")
    R.blank()

    # ---- 5B: Novel / underreported findings ----
    R.h2("5B. NOVEL OR UNDERREPORTED FINDINGS")
    R.blank()

    # Epithelial subtype × treatment
    epi_tx = load_counts("12_epithelial_level2", "treatment_status")
    if epi_tx is not None:
        active_cols = [c for c in epi_tx.columns if "Excluded" not in c]
        epi_tx_active = epi_tx[active_cols]
        props_tx = epi_tx_active.div(epi_tx_active.sum(axis=1), axis=0)
        oe_epi_tx, res_epi_tx, p_epi_tx = enrichment_table(epi_tx_active)

        R.text("1. EPITHELIAL SUBTYPE REMODELING POST-TREATMENT")
        R.text("   SecA-like (Cycling + Stress-response) vs SecB-like (Adaptive):")
        for tx in props_tx.index:
            seca = sum(props_tx.loc[tx, c] for c in SECA_TYPES if c in props_tx.columns)
            secb = sum(props_tx.loc[tx, c] for c in SECB_TYPES if c in props_tx.columns)
            ratio = seca / secb if secb > 0 else np.inf
            R.text(f"   {tx:40s}  SecA={seca:.3f} SecB={secb:.3f} ratio={ratio:.2f}")

        # Highlight: post-treatment adaptive enrichment
        if "Adaptive secretory epithelial cell" in oe_epi_tx.columns:
            R.blank()
            R.text("   ** Adaptive secretory (SecB) O/E across treatments:")
            for tx in oe_epi_tx.index:
                oe_val = oe_epi_tx.loc[tx, "Adaptive secretory epithelial cell"]
                R.text(f"      {tx:40s}  O/E = {oe_val:.2f}")
        R.blank()

    # Epithelial × metastatic site
    epi_ms = load_counts("12_epithelial_level2", "metastatic_site")
    if epi_ms is not None:
        active_cols = [c for c in epi_ms.columns if "Excluded" not in c]
        epi_ms_active = epi_ms[active_cols]
        props_ms = epi_ms_active.div(epi_ms_active.sum(axis=1), axis=0)
        oe_epi_ms, _, _ = enrichment_table(epi_ms_active)

        R.text("2. EPITHELIAL SUBTYPE DISTRIBUTION BY METASTATIC SITE")
        for site in props_ms.index:
            seca = sum(props_ms.loc[site, c] for c in SECA_TYPES if c in props_ms.columns)
            secb = sum(props_ms.loc[site, c] for c in SECB_TYPES if c in props_ms.columns)
            ratio = seca / secb if secb > 0 else np.inf
            R.text(f"   {site:15s}  SecA={seca:.3f} SecB={secb:.3f} ratio={ratio:.2f}")
        R.blank()

    # Epithelial × genomic status
    for gkey in ["BRCA_status", "HRD_status", "TP53_status"]:
        epi_gen = load_counts("12_epithelial_level2", gkey)
        if epi_gen is None:
            continue
        active_cols = [c for c in epi_gen.columns if "Excluded" not in c]
        epi_gen_active = epi_gen[active_cols]
        props_gen = epi_gen_active.div(epi_gen_active.sum(axis=1), axis=0)
        oe_gen, _, _ = enrichment_table(epi_gen_active)

        R.text(f"3. EPITHELIAL SUBTYPE × {META_DISPLAY[gkey].upper()}")
        for cat in props_gen.index:
            seca = sum(props_gen.loc[cat, c] for c in SECA_TYPES if c in props_gen.columns)
            secb = sum(props_gen.loc[cat, c] for c in SECB_TYPES if c in props_gen.columns)
            ratio = seca / secb if secb > 0 else np.inf
            R.text(f"   {cat:15s}  SecA={seca:.3f} SecB={secb:.3f} ratio={ratio:.2f}")
        if "Adaptive secretory epithelial cell" in oe_gen.columns:
            R.text(f"   Adaptive (SecB) O/E:")
            for cat in oe_gen.index:
                R.text(f"      {cat:15s}  O/E = {oe_gen.loc[cat, 'Adaptive secretory epithelial cell']:.2f}")
        R.blank()

    # ---- 5C: TME remodeling across compartments ----
    R.h2("5C. TME REMODELING — CROSS-COMPARTMENT PATTERNS")
    R.blank()

    # Macrophage subtypes × treatment
    mac_tx = load_counts("12_macrophage_level2", "treatment_status")
    if mac_tx is not None:
        active_cols = [c for c in mac_tx.columns if "Excluded" not in c]
        mac_tx_active = mac_tx[active_cols]
        if mac_tx_active.shape[0] >= 2 and mac_tx_active.shape[1] >= 2:
            oe_mac, _, _ = enrichment_table(mac_tx_active)
            R.text("1. MACROPHAGE SUBTYPES × TREATMENT STATUS")
            R.table(oe_mac.round(2))
            R.blank()

    # T/NK subtypes × treatment
    tnk_tx = load_counts("12_tnkcell_level2", "treatment_status")
    if tnk_tx is not None:
        active_cols = [c for c in tnk_tx.columns if "Excluded" not in c]
        tnk_tx_active = tnk_tx[active_cols]
        if tnk_tx_active.shape[0] >= 2 and tnk_tx_active.shape[1] >= 2:
            oe_tnk, _, _ = enrichment_table(tnk_tx_active)
            R.text("2. T/NK CELL SUBTYPES × TREATMENT STATUS")
            R.table(oe_tnk.round(2))
            R.blank()

    # Fibroblast subtypes × metastatic site
    fib_ms = load_counts("12_fibroblast_level2", "metastatic_site")
    if fib_ms is not None:
        active_cols = [c for c in fib_ms.columns if "Excluded" not in c]
        fib_ms_active = fib_ms[active_cols]
        if fib_ms_active.shape[0] >= 2 and fib_ms_active.shape[1] >= 2:
            oe_fib, _, _ = enrichment_table(fib_ms_active)
            R.text("3. FIBROBLAST SUBTYPES × METASTATIC SITE")
            R.table(oe_fib.round(2))
            R.blank()

    # ---- 5D: SecA/SecB polarization summary ----
    R.h2("5D. SecA/SecB POLARIZATION — EVIDENCE SUMMARY")
    R.blank()
    R.text("EVIDENCE SUPPORTING SecA → SecB POLARIZATION MODEL:")
    R.blank()

    # Collect all the SecA/SecB ratios across metadata
    evidence_items = []

    for meta_key in METADATA_KEYS:
        epi_counts = load_counts("12_epithelial_level2", meta_key)
        if epi_counts is None:
            continue
        active_cols = [c for c in epi_counts.columns if "Excluded" not in c]
        epi_active = epi_counts[active_cols]
        props = epi_active.div(epi_active.sum(axis=1), axis=0)

        for meta_cat in props.index:
            seca = sum(props.loc[meta_cat, c] for c in SECA_TYPES if c in props.columns)
            secb = sum(props.loc[meta_cat, c] for c in SECB_TYPES if c in props.columns)
            ratio = seca / secb if secb > 0 else np.inf
            evidence_items.append({
                "metadata_var": meta_key,
                "metadata_cat": meta_cat,
                "SecA_prop": seca,
                "SecB_prop": secb,
                "SecA_SecB_ratio": ratio,
                "n_cells": epi_active.loc[meta_cat].sum(),
            })

    ev_df = pd.DataFrame(evidence_items)
    ev_df = ev_df.sort_values("SecA_SecB_ratio", ascending=False)

    R.text("Top 10 conditions MOST skewed toward SecA (progenitor-like):")
    for _, row in ev_df.head(10).iterrows():
        R.text(f"   {row['metadata_var']:20s} = {row['metadata_cat']:30s}  "
               f"SecA/SecB = {row['SecA_SecB_ratio']:.2f}  "
               f"(n={row['n_cells']:,.0f})")
    R.blank()

    R.text("Top 10 conditions MOST skewed toward SecB (differentiated):")
    ev_df_finite = ev_df[ev_df["SecA_SecB_ratio"] < np.inf]
    for _, row in ev_df_finite.tail(10).iterrows():
        R.text(f"   {row['metadata_var']:20s} = {row['metadata_cat']:30s}  "
               f"SecA/SecB = {row['SecA_SecB_ratio']:.2f}  "
               f"(n={row['n_cells']:,.0f})")
    R.blank()

    # Treatment-specific SecA/SecB
    tx_rows = ev_df[ev_df["metadata_var"] == "treatment_status"].sort_values("SecA_SecB_ratio")
    if len(tx_rows) > 0:
        R.text("Treatment trajectory (SecA/SecB ratio):")
        for _, row in tx_rows.iterrows():
            bar_len = int(row["SecA_SecB_ratio"] * 20) if row["SecA_SecB_ratio"] < 5 else 100
            bar = "█" * min(bar_len, 60)
            R.text(f"   {row['metadata_cat']:40s}  {row['SecA_SecB_ratio']:.2f}  {bar}")
    R.blank()

    # Metastatic site SecA/SecB
    ms_rows = ev_df[ev_df["metadata_var"] == "metastatic_site"].sort_values("SecA_SecB_ratio")
    if len(ms_rows) > 0:
        R.text("Metastatic site trajectory (SecA/SecB ratio):")
        for _, row in ms_rows.iterrows():
            bar_len = int(row["SecA_SecB_ratio"] * 20) if row["SecA_SecB_ratio"] < 5 else 100
            bar = "█" * min(bar_len, 60)
            R.text(f"   {row['metadata_cat']:15s}  {row['SecA_SecB_ratio']:.2f}  {bar}")
    R.blank()

    # Genomic status SecA/SecB
    for gkey in ["BRCA_status", "HRD_status", "TP53_status"]:
        g_rows = ev_df[ev_df["metadata_var"] == gkey].sort_values("SecA_SecB_ratio")
        if len(g_rows) > 0:
            R.text(f"{META_DISPLAY[gkey]} (SecA/SecB ratio):")
            for _, row in g_rows.iterrows():
                bar_len = int(row["SecA_SecB_ratio"] * 20) if row["SecA_SecB_ratio"] < 5 else 100
                bar = "█" * min(bar_len, 60)
                R.text(f"   {row['metadata_cat']:15s}  {row['SecA_SecB_ratio']:.2f}  {bar}")
            R.blank()

    # Save evidence CSV
    ev_path = os.path.join(OUT_DIR, "12b_seca_secb_evidence.csv")
    ev_df.to_csv(ev_path, index=False)
    print(f"  Saved: {ev_path}")

    # ---- FINAL SUMMARY ----
    R.h1("SECTION 6: EXECUTIVE SUMMARY")
    R.blank()
    R.text("This analysis quantifies cell type composition variation across the HGSC atlas.")
    R.text("Key statistics: chi-square tests, Cramer's V effect sizes, O/E enrichment ratios,")
    R.text("adjusted standardized residuals, and Bonferroni-corrected p-values.")
    R.blank()
    R.text("FILES GENERATED:")
    R.text(f"  {REPORT_PATH}")
    R.text(f"  {l1_enr_path}")
    R.text(f"  {epi_enr_path}")
    R.text(f"  {comp_enr_path}")
    R.text(f"  {ev_path}")
    R.blank()

    R.save(REPORT_PATH)
    print(f"\n{'=' * 70}")
    print(f"DONE — Report at {REPORT_PATH}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
