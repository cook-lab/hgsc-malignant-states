#!/usr/bin/env python3
"""
Atlas 02 — Step 12c (counts): absolute cell-type counts + figures

PURPOSE
    Absolute (not proportional) cell-type counts across metadata strata, with the
    same epithelial breakdown as 12b; renders count figures + tables.

INPUTS
    obj("atlas_celltype_l2")  = hgsc_atlas_celltype_level2.h5ad

OUTPUTS
    output_root/02_annotation/12c_celltype_counts/{figs,tables}/*, report PDF

MANUSCRIPT PANEL(S)
    Supporting counts for Fig 1B-E composition.

RUNTIME TIER
    moderate (backed read of obs; tallies + plots).
"""

import argparse
import gc
import os
import sys
import time
import warnings
from collections import OrderedDict
from datetime import datetime

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

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
OUT_DIR     = path("output_root", "02_annotation", "12c_celltype_counts")
FIG_DIR     = os.path.join(OUT_DIR, "figs")
TABLE_DIR   = os.path.join(OUT_DIR, "tables")
PDF_PATH    = os.path.join(OUT_DIR, "12c_celltype_counts_report.pdf")

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

CELLTYPE_PALETTE = {
    "Epithelial":    "#E6A141",
    "Mesothelial":   "#A8A298",
    "Fibroblast":    "#DDD5CA",
    "Smooth muscle": "#D14E6C",
    "Pericyte":      "#B87A7A",
    "Endothelial":   "#7D4E4E",
    "T/NK cell":     "#87CEFA",
    "B cell":        "#5665B6",
    "Plasma cell":   "#8A5DAF",
    "Macrophage":    "#8FBC8F",
    "DC":            "#2E8B57",
    "Neutrophil":    "#6B8E23",
    "Mast cell":     "#8B9B6B",
}

CELLTYPE_ORDER = [
    "Epithelial", "Mesothelial",
    "Fibroblast", "Smooth muscle",
    "Pericyte", "Endothelial",
    "T/NK cell", "B cell", "Plasma cell",
    "Macrophage", "DC", "Neutrophil", "Mast cell",
]

KELLY_22 = [
    "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032",
    "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
    "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
    "#F2F3F4", "#222222",
]

STUDY_ORDER = [
    "denisenko_2024", "geistlinger_2020", "hornburg_2021",
    "loret_2022",     "luo_2024",         "nath_2021",
    "olalekan_2021",  "olbrecht_2021",    "regner_2021",
    "vazquez_garcia_2022", "xu_2022",     "zhang_2022",
    "zheng_2023",
]
STUDY_PALETTE = {s: KELLY_22[i] for i, s in enumerate(STUDY_ORDER)}

TREATMENT_PALETTE = {
    "pre-treatment":               "#DDD5CA",
    "post-treatment":              "#56B4E9",
    "post-chemotherapy":           "#009E73",
    "post-chemotherapy_niraparib": "#E69F00",
    "post-chemotherapy_olaparib":  "#0072B2",
    "post-chemotherapy_pembro":    "#D55E00",
    "post-niraparib":              "#CC79A7",
    "NA/Unknown":                  "#999999",
}

METASTATIC_PALETTE = {
    "primary":    "#7A9EBF",
    "metastasis": "#B07AA1",
    "ascites":    "#8FAC8C",
    "healthy":    "#C2956B",
}

LEVEL2_COLOR_SEQUENCE = [
    "#E6A141",   # golden orange   (Epithelial)
    "#5665B6",   # purple-blue     (B cell)
    "#D14E6C",   # rose pink       (Smooth muscle)
    "#2E8B57",   # forest green    (DC)
    "#87CEFA",   # sky blue        (T/NK cell)
    "#7D4E4E",   # dark brown      (Endothelial)
    "#8A5DAF",   # purple          (Plasma cell)
    "#8FBC8F",   # medium green    (Macrophage)
    "#D4A574",   # warm sand       (Mesothelial)
    "#6B8E23",   # olive drab      (Neutrophil)
    "#B87A7A",   # mauve           (Pericyte)
    "#C4B9A8",   # warm taupe      (Fibroblast)
    "#8B9B6B",   # sage green      (Mast cell)
]

EPITHELIAL_PALETTE = {
    "Adaptive secretory epithelial cell":        "#B8741A",
    "Ciliated epithelial cell":                  "#E05A2C",
    "Cycling secretory epithelial cell":         "#F6D28B",
    "Secretory epithelial cell":                 "#E6A141",
    "Stress-response secretory epithelial cell": "#D9C5A2",
    "Transitioning epithelial cell":             "#7D4E4E",
}

NA_COLOR = "#999999"

# ============================================================================
# COMPARTMENT REGISTRY
# ============================================================================

COMPARTMENTS = OrderedDict()
COMPARTMENTS["mastcell"]     = {"level1": "Mast cell"}
COMPARTMENTS["neutrophil"]   = {"level1": "Neutrophil"}
COMPARTMENTS["pericyte"]     = {"level1": "Pericyte"}
COMPARTMENTS["dc"]           = {"level1": "DC"}
COMPARTMENTS["plasmacell"]   = {"level1": "Plasma cell"}
COMPARTMENTS["smoothmuscle"] = {"level1": "Smooth muscle"}
COMPARTMENTS["bcell"]        = {"level1": "B cell"}
COMPARTMENTS["mesothelial"]  = {"level1": "Mesothelial"}
COMPARTMENTS["endothelial"]  = {"level1": "Endothelial"}
COMPARTMENTS["fibroblast"]   = {"level1": "Fibroblast"}
COMPARTMENTS["macrophage"]   = {"level1": "Macrophage"}
COMPARTMENTS["tnkcell"]      = {"level1": "T/NK cell"}
COMPARTMENTS["epithelial"]   = {"level1": "Epithelial"}

# ============================================================================
# METADATA VARIABLE REGISTRY
# ============================================================================

METADATA_VARS = OrderedDict([
    ("study",              {"display": "Study",              "palette": "STUDY_PALETTE"}),
    ("dataset",            {"display": "Dataset",            "palette": "auto"}),
    ("anatomic_site",      {"display": "Anatomic Site",      "palette": "auto"}),
    ("treatment_status",   {"display": "Treatment Status",   "palette": "TREATMENT_PALETTE"}),
    ("treatment_response", {"display": "Treatment Response", "palette": "auto"}),
    ("stage",              {"display": "Stage",              "palette": "auto"}),
    ("metastatic_site",    {"display": "Metastatic Site",    "palette": "METASTATIC_PALETTE"}),
    ("BRCA_status",        {"display": "BRCA Status",        "palette": "auto"}),
    ("HRD_status",         {"display": "HRD Status",         "palette": "auto"}),
    ("TP53_status",        {"display": "TP53 Status",        "palette": "auto"}),
    ("patient_id",         {"display": "Patient ID",         "palette": None}),
])

BOOLEAN_LIKE_COLS = {"BRCA_status", "HRD_status", "TP53_status"}


# ============================================================================
# DATA LOADING
# ============================================================================

def load_obs_metadata():
    """Load obs metadata from h5ad in backed mode. No gene expression needed.
    Excludes any cells labelled 'Excluded' at level1 or level2."""
    print(f"\nLoading (backed): {ATLAS_H5AD}", flush=True)
    t0 = time.time()

    adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
    n_cells, n_genes = adata.shape
    print(f"  Shape: {n_cells:,} cells × {n_genes:,} genes")

    cols = ["celltype_level1", "celltype_level2"] + list(METADATA_VARS.keys())
    obs = adata.obs[cols].copy()

    adata.file.close()
    del adata
    gc.collect()

    # Filter out "Excluded" cells (labels like "Excluded epithelial cell_1", etc.)
    n_before = len(obs)
    mask_excl = (
        obs["celltype_level1"].astype(str).str.startswith("Excluded") |
        obs["celltype_level2"].astype(str).str.startswith("Excluded")
    )
    n_excluded = mask_excl.sum()
    obs = obs.loc[~mask_excl].copy()
    print(f"  Removed {n_excluded:,} 'Excluded' cells "
          f"({n_before:,} → {len(obs):,})")

    print(f"  Extracted {len(obs):,} cells in {time.time()-t0:.0f}s\n")
    return obs


# ============================================================================
# NA HANDLING
# ============================================================================

def clean_metadata(obs, col):
    """
    Unify NA values: 'nan', 'NA', 'None', '' → 'NA/Unknown'.
    For boolean-like columns, also treat 'False' as 'NA/Unknown'.
    """
    s = obs[col].astype(str).copy()
    na_vals = {"nan", "NA", "None", ""}
    if col in BOOLEAN_LIKE_COLS:
        na_vals.add("False")
    s = s.where(~s.isin(na_vals), "NA/Unknown")
    return s


# ============================================================================
# COUNT COMPUTATION
# ============================================================================

def compute_counts(obs, group_col, category_col, category_order=None):
    """
    Compute cell type counts per group.

    Returns counts_df — rows=groups, cols=categories, values are integers.
    """
    counts = pd.crosstab(obs[group_col], obs[category_col])

    # Remove rows with zero total cells
    nonzero = counts.sum(axis=1) > 0
    counts = counts.loc[nonzero]

    # Order columns
    if category_order is not None:
        present = [c for c in category_order if c in counts.columns]
        extra = sorted([c for c in counts.columns if c not in category_order])
        col_order = present + extra
    else:
        col_order = counts.sum().sort_values(ascending=False).index.tolist()

    counts = counts[col_order]
    return counts


# ============================================================================
# PALETTE HELPERS
# ============================================================================

def get_metadata_palette(meta_key, categories):
    """Return a dict mapping category → color for a metadata variable."""
    info = METADATA_VARS[meta_key]
    pal_name = info["palette"]

    if pal_name == "STUDY_PALETTE":
        pal = dict(STUDY_PALETTE)
    elif pal_name == "TREATMENT_PALETTE":
        pal = dict(TREATMENT_PALETTE)
    elif pal_name == "METASTATIC_PALETTE":
        pal = dict(METASTATIC_PALETTE)
    elif pal_name == "auto":
        cats_sorted = sorted([c for c in categories if c != "NA/Unknown"])
        pal = {c: KELLY_22[i % len(KELLY_22)] for i, c in enumerate(cats_sorted)}
    else:
        pal = {}

    if "NA/Unknown" in categories:
        pal["NA/Unknown"] = NA_COLOR
    return pal


def get_level2_palette(compartment_level1, categories):
    """Return a dict mapping level2 categories → 9e-consistent colors."""
    if compartment_level1 == "Epithelial":
        pal = {c: EPITHELIAL_PALETTE.get(c, "#CCCCCC") for c in categories
               if c != "NA/Unknown"}
    else:
        cats_sorted = sorted([c for c in categories if c != "NA/Unknown"])
        pal = {c: LEVEL2_COLOR_SEQUENCE[i % len(LEVEL2_COLOR_SEQUENCE)]
               for i, c in enumerate(cats_sorted)}
    if "NA/Unknown" in categories:
        pal["NA/Unknown"] = NA_COLOR
    return pal


def abbreviate_label(label, max_len=30):
    """Truncate a label to max_len chars with ellipsis."""
    if len(str(label)) > max_len:
        return str(label)[:max_len - 1] + "…"
    return str(label)


# ============================================================================
# PLOTTING — GROUPED HORIZONTAL BAR
# ============================================================================

def plot_grouped_bar(counts, palette, title, xlabel="Cell count",
                     ylabel=None, annotations=None):
    """
    Grouped horizontal bar chart (side-by-side bars per group).

    counts: DataFrame, rows=groups, cols=categories (integer counts)
    palette: dict mapping category → color
    annotations: optional dict mapping group → subtitle string
    Returns matplotlib Figure.
    """
    n_groups = len(counts)
    n_cats = len(counts.columns)

    if n_cats == 0 or n_groups == 0:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    has_annot = annotations is not None and len(annotations) > 0
    bar_height = 0.8 / n_cats
    group_spacing = max(1.0, 0.8 + n_cats * 0.05)
    fig_h = max(4, n_groups * group_spacing * 0.45)
    fig_w = max(8, n_groups * 0.15 + 6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    group_positions = np.arange(n_groups) * group_spacing
    offsets = np.linspace(-0.4 + bar_height / 2,
                          0.4 - bar_height / 2,
                          n_cats)

    for j, cat in enumerate(counts.columns):
        color = palette.get(cat, "#CCCCCC")
        vals = counts[cat].values
        y_pos = group_positions + offsets[j]
        ax.barh(y_pos, vals, height=bar_height, color=color,
                edgecolor="white", linewidth=0.2, label=cat)

    ax.set_yticks(group_positions)
    tick_fs = max(5, min(7, 200 // max(n_groups, 1)))
    if has_annot:
        labels = []
        for g in counts.index:
            lbl = abbreviate_label(g)
            if g in annotations:
                lbl = f"{lbl}\n{annotations[g]}"
            labels.append(lbl)
        ax.set_yticklabels(labels, fontsize=tick_fs, linespacing=1.4)
    else:
        labels = [abbreviate_label(g) for g in counts.index]
        ax.set_yticklabels(labels, fontsize=tick_fs)

    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10, fontweight="bold", pad=10)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Comma-separated x-axis tick labels
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda x, p: f"{int(x):,}"))

    # Legend to the right
    handles, labels_leg = ax.get_legend_handles_labels()
    ax.legend(handles, labels_leg, bbox_to_anchor=(1.02, 1), loc="upper left",
              fontsize=6, frameon=False,
              ncol=1 if n_cats <= 15 else 2)

    fig.tight_layout()
    return fig


# ============================================================================
# PLOTTING — PATIENT HEATMAP (log counts)
# ============================================================================

def plot_patient_heatmap(counts, title, celltype_colors=None):
    """
    Clustered heatmap for patient_id × cell type absolute counts (log10+1).
    Returns a seaborn ClusterGrid figure.
    """
    log_counts = np.log10(counts + 1)

    col_colors = None
    if celltype_colors:
        col_colors = pd.Series(
            {c: celltype_colors.get(c, "#CCCCCC") for c in counts.columns},
            name="Cell type"
        )

    n_patients = len(counts)
    n_cols = len(counts.columns)
    fig_h = max(8, n_patients * 0.08 + 2)

    do_row_cluster = n_patients > 1
    do_col_cluster = n_cols > 1

    g = sns.clustermap(
        log_counts,
        cmap="YlOrRd",
        figsize=(10, min(fig_h, 30)),
        col_colors=col_colors,
        row_cluster=do_row_cluster,
        col_cluster=do_col_cluster,
        method="ward",
        metric="euclidean",
        linewidths=0,
        xticklabels=True,
        yticklabels=True if n_patients <= 60 else False,
        cbar_kws={"label": "log₁₀(count + 1)", "shrink": 0.5},
        dendrogram_ratio=(0.1, 0.08),
    )
    g.fig.suptitle(title, fontsize=11, fontweight="bold", y=1.02)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("Patient ID" if n_patients <= 60 else
                            f"Patient ID (n={n_patients})")

    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right", fontsize=6)
    if n_patients <= 60:
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=5)

    return g.fig


def plot_patient_top30_bar(counts, palette, title):
    """Grouped bar for top-30 patients by cell count (supplemental)."""
    if len(counts) <= 30:
        return plot_grouped_bar(counts, palette, title)
    top30 = counts.iloc[:30]
    return plot_grouped_bar(top30, palette,
                            title + " (Top 30 by cell count)")


# ============================================================================
# PDF HELPER
# ============================================================================

def save_figure(fig, pdf, fig_name):
    """Save figure to master PDF and as individual PDF in figs/."""
    pdf.savefig(fig, bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, fig_name), bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# TITLE PAGE
# ============================================================================

def render_title_page(pdf, obs):
    """Render title page with atlas summary stats."""
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    ax.axis("off")

    ax.text(0.5, 0.90, "Cell Type Absolute Count Analysis",
            transform=ax.transAxes, fontsize=20, fontweight="bold",
            ha="center", va="top")
    ax.text(0.5, 0.85, "HGSC Atlas | Cook Lab | Step 12c",
            transform=ax.transAxes, fontsize=12, ha="center", va="top",
            color="#555555")

    n_cells = len(obs)
    n_level1 = obs["celltype_level1"].nunique()
    n_level2 = obs["celltype_level2"].nunique()
    n_patients = obs["patient_id"].astype(str).nunique()
    n_studies = obs["study"].astype(str).nunique()

    summary_text = (
        f"Total cells: {n_cells:,}\n"
        f"Level 1 cell types: {n_level1}\n"
        f"Level 2 cell types: {n_level2}\n"
        f"Studies: {n_studies}\n"
        f"Patients: {n_patients}\n"
        f"Excluded cells removed: yes\n"
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}"
    )
    ax.text(0.5, 0.72, summary_text, transform=ax.transAxes, fontsize=10,
            ha="center", va="top", family="monospace",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#F5F5F5",
                      edgecolor="#CCCCCC"))

    # Metadata variable summary table
    table_data = []
    for key, info in METADATA_VARS.items():
        cleaned = clean_metadata(obs, key)
        cats = sorted(cleaned.unique())
        n_cat = len(cats)
        n_na = (cleaned == "NA/Unknown").sum()
        pct_na = n_na / n_cells * 100
        table_data.append([
            info["display"], key, str(n_cat),
            f"{n_na:,} ({pct_na:.1f}%)"
        ])

    col_labels = ["Variable", "Key", "Categories", "NA/Unknown"]
    table = ax.table(
        cellText=table_data, colLabels=col_labels,
        cellLoc="left", colLoc="left", loc="center",
        bbox=[0.1, 0.05, 0.8, 0.45],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    for j in range(len(col_labels)):
        cell = table[0, j]
        cell.set_facecolor("#2C3E50")
        cell.set_text_props(color="white", fontweight="bold")
    for i in range(1, len(table_data) + 1):
        for j in range(len(col_labels)):
            table[i, j].set_facecolor("#F9F9F9" if i % 2 == 0 else "white")

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    summary_df = pd.DataFrame(table_data, columns=col_labels)
    summary_df.to_csv(os.path.join(TABLE_DIR, "12c_metadata_summary.csv"),
                      index=False)
    print("  Saved metadata summary table")


# ============================================================================
# PART 1 — LEVEL 1 COUNTS
# ============================================================================

def render_part1(pdf, obs, skip_dataset=False):
    """Part 1: celltype_level1 absolute counts across each metadata variable."""
    print("\n" + "=" * 70)
    print("PART 1: Level 1 cell type absolute counts")
    print("=" * 70)

    # Section header page
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.55, "Part 1: celltype_level1 Absolute Counts",
            transform=ax.transAxes, fontsize=18, fontweight="bold",
            ha="center", va="center")
    ax.text(0.5, 0.45, f"{len(CELLTYPE_ORDER)} cell types × "
            f"{len(METADATA_VARS)} metadata variables",
            transform=ax.transAxes, fontsize=12, ha="center", va="center",
            color="#777777")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    fig_count = 0

    for meta_key, meta_info in METADATA_VARS.items():
        display = meta_info["display"]

        if meta_key == "dataset" and skip_dataset:
            print(f"  Skipping {display} (identical to study)")
            continue

        print(f"  Processing: {display} ({meta_key})")
        t0 = time.time()

        work = obs[["celltype_level1"]].copy()
        work["_group"] = clean_metadata(obs, meta_key)

        counts = compute_counts(work, "_group", "celltype_level1",
                                category_order=CELLTYPE_ORDER)

        if len(counts) == 0:
            print(f"    No data — skipping")
            continue

        csv_name = f"12c_level1_counts_{meta_key}.csv"
        counts.to_csv(os.path.join(TABLE_DIR, csv_name))

        if meta_key == "patient_id":
            cell_totals = counts.sum(axis=1).sort_values(ascending=False)
            counts_sorted = counts.loc[cell_totals.index]

            # Heatmap (log counts)
            fig_hm = plot_patient_heatmap(
                counts, f"Level 1 Counts × {display}",
                celltype_colors=CELLTYPE_PALETTE
            )
            hm_name = f"12c_level1_vs_{meta_key}_heatmap.pdf"
            save_figure(fig_hm, pdf, hm_name)
            fig_count += 1

            # Top-30 grouped bar
            fig_bar = plot_patient_top30_bar(
                counts_sorted, CELLTYPE_PALETTE,
                f"Level 1 Counts × {display}"
            )
            bar_name = f"12c_level1_vs_{meta_key}_top30.pdf"
            save_figure(fig_bar, pdf, bar_name)
            fig_count += 1

        else:
            fig_bar = plot_grouped_bar(
                counts, CELLTYPE_PALETTE,
                f"Level 1 Counts × {display}",
                ylabel=display
            )
            fig_name = f"12c_level1_vs_{meta_key}.pdf"
            save_figure(fig_bar, pdf, fig_name)
            fig_count += 1

        elapsed = time.time() - t0
        print(f"    {len(counts)} groups, {len(counts.columns)} types "
              f"({elapsed:.1f}s)")

    print(f"\n  Part 1 complete: {fig_count} figures")
    return fig_count


# ============================================================================
# PART 2 — LEVEL 2 COUNTS PER COMPARTMENT
# ============================================================================

def render_part2(pdf, obs, compartment_key=None):
    """
    Part 2: celltype_level2 absolute counts across metadata, split by compartment.
    Includes n= (cell count) and s= (sample count) annotations on y-axis labels.
    """
    print("\n" + "=" * 70)
    print("PART 2: Level 2 cell type absolute counts per compartment")
    print("=" * 70)

    # Section header page
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)
    ax.axis("off")
    ax.text(0.5, 0.55, "Part 2: celltype_level2 Absolute Counts by Compartment",
            transform=ax.transAxes, fontsize=18, fontweight="bold",
            ha="center", va="center")
    comp_label = compartment_key if compartment_key else "all compartments"
    ax.text(0.5, 0.45, f"Scope: {comp_label}",
            transform=ax.transAxes, fontsize=12, ha="center", va="center",
            color="#777777")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    fig_count = 0

    if compartment_key:
        if compartment_key not in COMPARTMENTS:
            print(f"  ERROR: Unknown compartment '{compartment_key}'")
            print(f"  Available: {', '.join(COMPARTMENTS.keys())}")
            return 0
        comp_iter = [(compartment_key, COMPARTMENTS[compartment_key])]
    else:
        comp_iter = list(COMPARTMENTS.items())

    for comp_key, comp_info in comp_iter:
        level1_name = comp_info["level1"]
        print(f"\n  Compartment: {level1_name} ({comp_key})")

        mask = obs["celltype_level1"] == level1_name
        obs_comp = obs.loc[mask].copy()
        n_comp = len(obs_comp)

        if n_comp == 0:
            print(f"    Empty — skipping")
            continue

        n_types = obs_comp["celltype_level2"].nunique()
        n_samples = obs_comp["patient_id"].astype(str).nunique()
        print(f"    {n_comp:,} cells, {n_types} level2 types, "
              f"{n_samples} patients")

        for meta_key, meta_info in METADATA_VARS.items():
            display = meta_info["display"]

            work = obs_comp[["celltype_level2"]].copy()
            work["_group"] = clean_metadata(obs_comp, meta_key)

            counts = compute_counts(work, "_group", "celltype_level2")

            if len(counts) == 0 or len(counts.columns) == 0:
                continue

            csv_name = f"12c_{comp_key}_level2_counts_{meta_key}.csv"
            counts.to_csv(os.path.join(TABLE_DIR, csv_name))

            pal = get_level2_palette(level1_name, counts.columns.tolist())

            title = (f"{level1_name} — Level 2 Counts × {display}"
                     f"  (n={n_comp:,}, s={n_samples})")

            if meta_key == "patient_id":
                if len(counts.columns) < 2:
                    continue
                fig_hm = plot_patient_heatmap(counts, title)
                hm_name = f"12c_{comp_key}_level2_vs_{meta_key}_heatmap.pdf"
                save_figure(fig_hm, pdf, hm_name)
                fig_count += 1
            else:
                # Per-group n= and s= annotations
                group_cells = counts.sum(axis=1)
                group_samples = (
                    obs_comp.assign(_group=clean_metadata(obs_comp, meta_key))
                    .groupby("_group")["patient_id"]
                    .nunique()
                )
                annot = {
                    g: f"(n={group_cells.get(g, 0):,}, s={group_samples.get(g, 0)})"
                    for g in counts.index
                }

                fig_bar = plot_grouped_bar(counts, pal, title, ylabel=display,
                                           annotations=annot)
                fig_name = f"12c_{comp_key}_level2_vs_{meta_key}.pdf"
                save_figure(fig_bar, pdf, fig_name)
                fig_count += 1

        print(f"    → {n_types} types plotted across metadata vars")

    print(f"\n  Part 2 complete: {fig_count} figures")
    return fig_count


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Step 12c — Cell Type Absolute Count Analysis"
    )
    parser.add_argument("--level1-only", action="store_true",
                        help="Only run Part 1 (skip per-compartment level2)")
    parser.add_argument("--compartment", type=str, default=None,
                        help="Run Part 2 for a single compartment "
                             f"({', '.join(COMPARTMENTS.keys())})")
    args = parser.parse_args()

    print("=" * 70)
    print("Step 12c — Cell Type Absolute Count Analysis")
    print("=" * 70)
    t_start = time.time()

    obs = load_obs_metadata()

    skip_dataset = False
    if "dataset" in obs.columns and "study" in obs.columns:
        ds = clean_metadata(obs, "dataset")
        st = clean_metadata(obs, "study")
        if ds.equals(st):
            print("  NOTE: 'dataset' is identical to 'study' — will skip dataset")
            skip_dataset = True
        else:
            n_diff = (ds != st).sum()
            print(f"  NOTE: 'dataset' differs from 'study' in {n_diff:,} cells")

    print(f"\nOutput PDF: {PDF_PATH}")

    with PdfPages(PDF_PATH) as pdf:
        render_title_page(pdf, obs)

        n1 = render_part1(pdf, obs, skip_dataset=skip_dataset)

        n2 = 0
        if not args.level1_only:
            n2 = render_part2(pdf, obs, compartment_key=args.compartment)
        else:
            print("\n  Skipping Part 2 (--level1-only)")

    total_figs = n1 + n2
    elapsed = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"DONE — {total_figs} figures in {elapsed:.0f}s")
    print(f"  PDF:    {PDF_PATH}")
    print(f"  Figs:   {FIG_DIR}/")
    print(f"  Tables: {TABLE_DIR}/")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
