#!/usr/bin/env python3
"""
Atlas 02 — Step 12d (epithelial violins): epithelial composition violins by metadata

PURPOSE
    Per-sample epithelial subtype proportion violins across metadata variables,
    using the SecA -> Intermediate -> SecB -> other display ordering.

INPUTS
    obj("atlas_celltype_l2")  = hgsc_atlas_celltype_level2.h5ad

OUTPUTS
    output_root/02_annotation/12d_epithelial_proportion_violins/*

MANUSCRIPT PANEL(S)
    Epithelial composition substrate for Fig 2D / SF9.

RUNTIME TIER
    moderate (backed read of obs; violins).
"""

import argparse
import gc
import os
import time
import warnings
from collections import OrderedDict
from datetime import datetime

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================================
# PATHS
# ============================================================================

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

ATLAS_H5AD  = obj("atlas_celltype_l2")
OUT_DIR     = path("output_root", "02_annotation", "12d_epithelial_proportion_violins")
FIG_DIR     = os.path.join(OUT_DIR, "figs")
TABLE_DIR   = os.path.join(OUT_DIR, "tables")
PDF_PATH    = os.path.join(OUT_DIR, "12d_epithelial_proportion_violins_report.pdf")

for d in [OUT_DIR, FIG_DIR, TABLE_DIR]:
    os.makedirs(d, exist_ok=True)

# ============================================================================
# COOK LAB v1.2 STYLE
# ============================================================================

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.dpi":      450,
    "savefig.dpi":     450,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
    "savefig.bbox":    "tight",
})

# ============================================================================
# COLOUR PALETTES
# ============================================================================

# Epithelial-specific palette — SecA/SecB polarization axis (from 9e)
EPITHELIAL_PALETTE = {
    "Adaptive secretory epithelial cell":        "#B8741A",
    "Ciliated epithelial cell":                  "#E05A2C",
    "Cycling secretory epithelial cell":         "#F6D28B",
    "Secretory epithelial cell":                 "#E6A141",
    "Stress-response secretory epithelial cell": "#D9C5A2",
    "Transitioning epithelial cell":             "#7D4E4E",
}

# Display order (by biological axis: SecA → intermediate → SecB → other)
EPITHELIAL_ORDER = [
    "Cycling secretory epithelial cell",
    "Secretory epithelial cell",
    "Adaptive secretory epithelial cell",
    "Stress-response secretory epithelial cell",
    "Ciliated epithelial cell",
    "Transitioning epithelial cell",
]

# Short labels for subplot titles
EPITHELIAL_SHORT = {
    "Cycling secretory epithelial cell":         "Cycling secretory",
    "Secretory epithelial cell":                 "Secretory",
    "Adaptive secretory epithelial cell":        "Adaptive secretory",
    "Stress-response secretory epithelial cell": "Stress-response",
    "Ciliated epithelial cell":                  "Ciliated",
    "Transitioning epithelial cell":             "Intermediate",
}

NA_COLOR = "#999999"

# ============================================================================
# METADATA VARIABLE REGISTRY
# ============================================================================

METADATA_VARS = OrderedDict([
    ("study",              {"display": "Study"}),
    ("anatomic_site",      {"display": "Anatomic Site"}),
    ("treatment_status",   {"display": "Treatment Status"}),
    ("treatment_response", {"display": "Treatment Response"}),
    ("stage",              {"display": "Stage"}),
    ("metastatic_site",    {"display": "Metastatic Site"}),
    ("BRCA_status",        {"display": "BRCA Status"}),
    ("HRD_status",         {"display": "HRD Status"}),
    ("TP53_status",        {"display": "TP53 Status"}),
])

BOOLEAN_LIKE_COLS = {"BRCA_status", "HRD_status", "TP53_status"}


# ============================================================================
# DATA LOADING
# ============================================================================

def load_obs_metadata():
    """Load obs metadata, filter to epithelial only, exclude 'Excluded' cells."""
    print(f"\nLoading (backed): {ATLAS_H5AD}", flush=True)
    t0 = time.time()

    adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
    n_cells, n_genes = adata.shape
    print(f"  Shape: {n_cells:,} cells × {n_genes:,} genes")

    cols = ["celltype_level1", "celltype_level2", "patient_id"] + list(METADATA_VARS.keys())
    obs = adata.obs[cols].copy()

    adata.file.close()
    del adata
    gc.collect()

    # Filter out "Excluded" cells
    n_before = len(obs)
    mask_excl = (
        obs["celltype_level1"].astype(str).str.startswith("Excluded") |
        obs["celltype_level2"].astype(str).str.startswith("Excluded")
    )
    n_excluded = mask_excl.sum()
    obs = obs.loc[~mask_excl].copy()
    print(f"  Removed {n_excluded:,} 'Excluded' cells "
          f"({n_before:,} → {len(obs):,})")

    # Filter to epithelial only
    mask_epi = obs["celltype_level1"] == "Epithelial"
    obs_epi = obs.loc[mask_epi].copy()
    n_epi = len(obs_epi)
    n_patients = obs_epi["patient_id"].astype(str).nunique()
    n_subtypes = obs_epi["celltype_level2"].nunique()
    print(f"  Epithelial subset: {n_epi:,} cells, {n_subtypes} subtypes, "
          f"{n_patients} patients")

    print(f"  Loaded in {time.time()-t0:.0f}s\n")
    return obs_epi


# ============================================================================
# NA HANDLING
# ============================================================================

def clean_metadata(obs, col):
    """Unify NA values → 'NA/Unknown'."""
    s = obs[col].astype(str).copy()
    na_vals = {"nan", "NA", "None", ""}
    if col in BOOLEAN_LIKE_COLS:
        na_vals.add("False")
    s = s.where(~s.isin(na_vals), "NA/Unknown")
    return s


# ============================================================================
# PER-PATIENT PROPORTIONS
# ============================================================================

def compute_patient_proportions(obs_epi):
    """
    Compute per-patient epithelial level2 proportions.

    Returns a DataFrame: rows=patient_id, cols=epithelial subtypes,
    values in [0,1] summing to 1 per row.
    """
    counts = pd.crosstab(obs_epi["patient_id"].astype(str),
                         obs_epi["celltype_level2"])
    props = counts.div(counts.sum(axis=1), axis=0)

    # Ensure all epithelial subtypes are present as columns
    for st in EPITHELIAL_ORDER:
        if st not in props.columns:
            props[st] = 0.0

    # Order columns
    present = [c for c in EPITHELIAL_ORDER if c in props.columns]
    extra = sorted([c for c in props.columns if c not in EPITHELIAL_ORDER])
    props = props[present + extra]

    return props, counts


# ============================================================================
# PLOTTING — VIOLIN GRID
# ============================================================================

def plot_violin_grid(patient_props, patient_meta, meta_key, display,
                     subtypes, n_epi, n_patients):
    """
    Create a 2×3 grid of violin plots (one per epithelial subtype).

    patient_props: DataFrame (patient_id × subtype proportions)
    patient_meta:  Series (patient_id → metadata category)
    subtypes:      list of subtype names to plot
    """
    n_subtypes = len(subtypes)
    nrows = 2
    ncols = 3

    # Determine categories and their order
    categories = sorted(patient_meta.unique())
    # Move NA/Unknown to end
    if "NA/Unknown" in categories:
        categories = [c for c in categories if c != "NA/Unknown"] + ["NA/Unknown"]
    n_cats = len(categories)

    # Compute per-category n= and s= for x-axis labels
    cat_labels = []
    for cat in categories:
        pids_in_cat = patient_meta[patient_meta == cat].index
        n_cells_cat = len(pids_in_cat)  # number of patients
        cat_labels.append(f"{cat}\n(s={n_cells_cat})")

    fig_w = max(10, n_cats * 0.8 + 4)
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, 8),
                             squeeze=False)

    fig.suptitle(f"Epithelial Level 2 Proportions × {display}"
                 f"  (n={n_epi:,}, s={n_patients})",
                 fontsize=11, fontweight="bold", y=1.02)

    for idx, st in enumerate(subtypes):
        row = idx // ncols
        col = idx % ncols
        ax = axes[row, col]

        color = EPITHELIAL_PALETTE.get(st, "#CCCCCC")
        short = EPITHELIAL_SHORT.get(st, st)

        # Collect data per category
        violin_data = []
        positions = []
        for i, cat in enumerate(categories):
            pids = patient_meta[patient_meta == cat].index
            pids_in_props = [p for p in pids if p in patient_props.index]
            if len(pids_in_props) > 0 and st in patient_props.columns:
                vals = patient_props.loc[pids_in_props, st].values
            else:
                vals = np.array([0.0])
            violin_data.append(vals)
            positions.append(i)

        # Draw violins
        if any(len(v) > 1 for v in violin_data):
            parts = ax.violinplot(violin_data, positions=positions,
                                  showmeans=False, showmedians=False,
                                  showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(color)
                pc.set_edgecolor("none")
                pc.set_alpha(0.7)

        # Overlay box plots (thin)
        bp = ax.boxplot(violin_data, positions=positions, widths=0.15,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", linewidth=1),
                        boxprops=dict(facecolor="white", edgecolor="black",
                                      linewidth=0.5),
                        whiskerprops=dict(linewidth=0.5),
                        capprops=dict(linewidth=0.5))

        # Strip plot (individual patients)
        rng = np.random.RandomState(SEED)
        for i, vals in enumerate(violin_data):
            if len(vals) > 0:
                jitter = rng.uniform(-0.12, 0.12, size=len(vals))
                ax.scatter(positions[i] + jitter, vals,
                           c=color, s=3, alpha=0.4, edgecolors="none",
                           zorder=3, rasterized=True)

        ax.set_xticks(positions)
        ax.set_xticklabels(cat_labels, fontsize=5.5, rotation=45,
                           ha="right")
        ax.set_ylabel("Proportion")
        ax.set_title(short, fontsize=8, fontweight="bold", color=color)
        ax.set_ylim(-0.02, 1.02)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Hide unused subplots
    for idx in range(n_subtypes, nrows * ncols):
        row = idx // ncols
        col = idx % ncols
        axes[row, col].set_visible(False)

    fig.tight_layout()
    return fig


def plot_summary_violins(patient_props, subtypes, n_epi, n_patients):
    """
    Summary figure: one violin per epithelial subtype (all patients pooled).
    """
    n_subtypes = len(subtypes)
    fig, ax = plt.subplots(figsize=(max(7, n_subtypes * 1.2), 5))

    violin_data = []
    colors = []
    labels = []
    for st in subtypes:
        if st in patient_props.columns:
            vals = patient_props[st].dropna().values
        else:
            vals = np.array([0.0])
        violin_data.append(vals)
        colors.append(EPITHELIAL_PALETTE.get(st, "#CCCCCC"))
        labels.append(EPITHELIAL_SHORT.get(st, st))

    positions = np.arange(n_subtypes)

    # Violins
    parts = ax.violinplot(violin_data, positions=positions,
                          showmeans=False, showmedians=False,
                          showextrema=False)
    for pc, color in zip(parts["bodies"], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor("none")
        pc.set_alpha(0.7)

    # Box overlay
    ax.boxplot(violin_data, positions=positions, widths=0.18,
               patch_artist=True, showfliers=False,
               medianprops=dict(color="black", linewidth=1),
               boxprops=dict(facecolor="white", edgecolor="black",
                             linewidth=0.5),
               whiskerprops=dict(linewidth=0.5),
               capprops=dict(linewidth=0.5))

    # Strip plot
    rng = np.random.RandomState(SEED)
    for i, (vals, color) in enumerate(zip(violin_data, colors)):
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(positions[i] + jitter, vals,
                   c=color, s=5, alpha=0.4, edgecolors="none",
                   zorder=3, rasterized=True)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
    ax.set_ylabel("Per-patient proportion\n(within epithelial)")
    ax.set_title(f"Epithelial Level 2 Proportions — All Patients"
                 f"  (n={n_epi:,}, s={n_patients})",
                 fontsize=10, fontweight="bold", pad=10)
    ax.set_ylim(-0.02, 1.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    return fig


# ============================================================================
# PDF HELPER
# ============================================================================

def save_figure(fig, pdf, fig_name):
    """Save figure to master PDF and as individual PDF in figs/."""
    pdf.savefig(fig, bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, fig_name), bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Step 12d — Epithelial Level 2 Per-Patient Proportion Violins")
    print("=" * 70)
    t_start = time.time()

    # Load data
    obs_epi = load_obs_metadata()
    n_epi = len(obs_epi)
    n_patients = obs_epi["patient_id"].astype(str).nunique()

    # Compute per-patient proportions
    print("Computing per-patient epithelial proportions...")
    patient_props, patient_counts = compute_patient_proportions(obs_epi)
    print(f"  {len(patient_props)} patients × "
          f"{len(patient_props.columns)} subtypes")

    # Save proportions table
    patient_props.to_csv(os.path.join(TABLE_DIR,
                                       "12d_patient_epithelial_proportions.csv"))
    patient_counts.to_csv(os.path.join(TABLE_DIR,
                                        "12d_patient_epithelial_counts.csv"))
    print("  Saved proportion and count tables")

    # Determine subtypes to plot (present in data, in preferred order)
    subtypes = [s for s in EPITHELIAL_ORDER
                if s in patient_props.columns and patient_props[s].sum() > 0]
    print(f"  Subtypes to plot: {len(subtypes)}")

    # Build per-patient metadata lookup (one row per patient)
    # Use the first occurrence for each patient (metadata is patient-level)
    patient_meta_df = (obs_epi.groupby(obs_epi["patient_id"].astype(str))
                       .first()
                       [list(METADATA_VARS.keys())])

    print(f"\nOutput PDF: {PDF_PATH}")

    with PdfPages(PDF_PATH) as pdf:

        # --- Title page ---
        fig = plt.figure(figsize=(11, 8.5))
        fig.patch.set_facecolor("white")
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.5, 0.90, "Epithelial Level 2 — Per-Patient Proportion Violins",
                transform=ax.transAxes, fontsize=18, fontweight="bold",
                ha="center", va="top")
        ax.text(0.5, 0.85, "HGSC Atlas | Cook Lab | Step 12d",
                transform=ax.transAxes, fontsize=12, ha="center", va="top",
                color="#555555")
        summary = (
            f"Epithelial cells: {n_epi:,}\n"
            f"Patients: {n_patients}\n"
            f"Subtypes: {len(subtypes)}\n"
            f"Excluded cells removed: yes\n"
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}"
        )
        ax.text(0.5, 0.70, summary, transform=ax.transAxes, fontsize=10,
                ha="center", va="top", family="monospace",
                bbox=dict(boxstyle="round,pad=0.5", facecolor="#F5F5F5",
                          edgecolor="#CCCCCC"))
        # Subtype colour legend
        for i, st in enumerate(subtypes):
            color = EPITHELIAL_PALETTE.get(st, "#CCCCCC")
            short = EPITHELIAL_SHORT.get(st, st)
            ax.text(0.30, 0.48 - i * 0.045, f"■  {short}  ({st})",
                    transform=ax.transAxes, fontsize=8, color=color,
                    fontweight="bold", va="top")

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig_count = 0

        # --- Summary violin (all patients pooled) ---
        print("\n  Summary violin plot (all patients)...")
        fig_sum = plot_summary_violins(patient_props, subtypes,
                                        n_epi, n_patients)
        save_figure(fig_sum, pdf, "12d_summary_violins.pdf")
        fig_count += 1

        # --- Per metadata variable ---
        for meta_key, meta_info in METADATA_VARS.items():
            display = meta_info["display"]
            print(f"  Processing: {display} ({meta_key})")
            t0 = time.time()

            # Get per-patient metadata
            meta_col = clean_metadata(patient_meta_df, meta_key)
            # meta_col is indexed by patient_id

            n_cats = meta_col.nunique()
            if n_cats == 0:
                print(f"    No categories — skipping")
                continue

            print(f"    {n_cats} categories")

            fig_v = plot_violin_grid(patient_props, meta_col, meta_key,
                                     display, subtypes, n_epi, n_patients)
            fig_name = f"12d_violins_vs_{meta_key}.pdf"
            save_figure(fig_v, pdf, fig_name)
            fig_count += 1

            print(f"    Done ({time.time()-t0:.1f}s)")

    elapsed = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"DONE — {fig_count} figures in {elapsed:.0f}s")
    print(f"  PDF:    {PDF_PATH}")
    print(f"  Figs:   {FIG_DIR}/")
    print(f"  Tables: {TABLE_DIR}/")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
