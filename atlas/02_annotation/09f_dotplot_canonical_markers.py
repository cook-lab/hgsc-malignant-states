#!/usr/bin/env python3
"""
Atlas 02 — Step 09f: canonical lineage marker dot plot (Level-1 cell types)

PURPOSE
    Custom dot plot of canonical lineage markers across the 13 level-1 cell types
    (dot size = fraction expressing; colour = per-gene normalised mean expression).
    Reads the final atlas in backed mode and extracts only the marker genes.

INPUTS
    obj("atlas_final")  = hgsc_atlas_final.h5ad

OUTPUTS
    output_root/02_annotation/09f_dotplot_canonical_markers/dotplot_stats.csv
    output_root/02_annotation/09f_dotplot_canonical_markers/lilra4_stats.csv   (LILRA4/pDC row for SF4B)
    output_root/02_annotation/09f_dotplot_canonical_markers/dotplot_canonical_markers.{svg,pdf}

MANUSCRIPT PANEL(S)
    SF4B (canonical marker dotplot; the figure script reads dotplot_stats.csv +
    lilra4_stats.csv).

RUNTIME TIER
    moderate (backed read; marker-gene subset into memory).
"""

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

# ── Cook Lab style v1.2 ─────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         8,
    "axes.titlesize":    9,
    "axes.labelsize":    8,
    "xtick.labelsize":   7,
    "ytick.labelsize":   8,
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
H5AD    = obj("atlas_final")
OUT_DIR = path("output_root", "02_annotation", "09f_dotplot_canonical_markers")
os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================================
# COLOUR PALETTES  (Cook Lab v1.2)
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

EXPR_CMAP = LinearSegmentedColormap.from_list(
    "atlas_expr", ["#5665B6", "#F7F7F7", "#D14E6C"], N=256
)


# ============================================================================
# CANONICAL MARKERS  (2 per compartment, ordered to match Y-axis)
# ============================================================================

MARKERS = {
    "Epithelial":    ["EPCAM", "PAX8"],
    "Mesothelial":   ["MSLN", "UPK3B"],
    "Fibroblast":    ["COL1A1", "DCN"],
    "Smooth muscle": ["MYH11", "ACTA2"],
    "Pericyte":      ["PDGFRB", "RGS5"],
    "Endothelial":   ["CDH5", "VWF"],
    "T/NK cell":     ["CD3E", "NKG7"],
    "B cell":        ["MS4A1", "CD79A"],
    "Plasma cell":   ["MZB1", "XBP1"],
    "Macrophage":    ["CD14", "C1QA"],
    "DC":            ["FCER1A", "CD1C"],
    "Neutrophil":    ["CSF3R", "S100A8"],
    "Mast cell":     ["KIT", "TPSAB1"],
}

GENE_ORDER = [g for ct in CELLTYPE_ORDER for g in MARKERS[ct]]

# LILRA4 (pDC marker) is not part of the celltype_level1 MARKERS panel, but SF4B
# appends it as a separate row; emitted to its own lilra4_stats.csv (same schema).
LILRA4_GENES = ["LILRA4"]


def load_marker_data():
    print("Opening atlas in backed mode ...")
    t0 = time.time()
    adata = sc.read_h5ad(H5AD, backed="r")
    print(f"  Loaded metadata in {time.time() - t0:.0f}s  —  "
          f"{adata.n_obs:,} cells × {adata.n_vars:,} genes")
    present = [g for g in GENE_ORDER if g in adata.var_names]
    missing = set(GENE_ORDER) - set(present)
    if missing:
        print(f"  Missing genes (will skip): {missing}")
    extra = [g for g in LILRA4_GENES if g in adata.var_names]
    load_genes = present + extra
    print(f"  Reading {len(load_genes)} marker genes into memory ...")
    small = adata[:, load_genes].to_memory()
    adata.file.close()
    return small, present


def compute_dotplot_stats(adata, genes, groupby="celltype_level1"):
    groups = adata.obs[groupby]
    X = adata[:, genes].X
    if sparse.issparse(X):
        X = X.toarray()
    records = []
    for ct in CELLTYPE_ORDER:
        mask = (groups == ct).values
        n_cells = mask.sum()
        if n_cells == 0:
            continue
        sub = X[mask, :]
        for j, gene in enumerate(genes):
            vals = sub[:, j]
            records.append({
                "celltype": ct, "gene": gene,
                "frac_expressing": np.mean(vals > 0),
                "mean_expression": np.mean(vals),
                "n_cells": n_cells,
            })
    df = pd.DataFrame(records)
    for gene in genes:
        mask = df["gene"] == gene
        vals = df.loc[mask, "mean_expression"]
        vmin, vmax = vals.min(), vals.max()
        df.loc[mask, "norm_expression"] = (vals - vmin) / (vmax - vmin) if vmax > vmin else 0.0
    return df


def plot_dotplot(df, genes, out_dir):
    n_genes = len(genes)
    n_types = len(CELLTYPE_ORDER)
    fig = plt.figure(figsize=(1.4 + n_genes * 0.34, 0.6 + n_types * 0.19))
    gs = gridspec.GridSpec(1, 2, width_ratios=[n_genes, 4], wspace=0.05, figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    ax_legend = fig.add_subplot(gs[0, 1])
    ax_legend.axis("off")

    SIZE_MIN, SIZE_MAX = 2, 90
    for _, row in df.iterrows():
        xi = genes.index(row["gene"])
        yi = CELLTYPE_ORDER.index(row["celltype"])
        size = SIZE_MIN + row["frac_expressing"] * (SIZE_MAX - SIZE_MIN)
        ax.scatter(xi, yi, s=size, c=[EXPR_CMAP(row["norm_expression"])],
                   edgecolors="#333333", linewidths=0.3, zorder=3)

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(genes, rotation=90, ha="center", fontsize=7, fontstyle="italic")
    ax.set_yticks(range(n_types))
    ylabels = ax.set_yticklabels(CELLTYPE_ORDER, fontsize=8)
    for lbl, ct in zip(ylabels, CELLTYPE_ORDER):
        lbl.set_color(CELLTYPE_PALETTE[ct])
        lbl.set_fontweight("bold")
    ax.set_xlim(-0.6, n_genes - 0.4)
    ax.set_ylim(-0.6, n_types - 0.4)
    ax.invert_yaxis()
    for y in range(n_types):
        ax.axhline(y, color="#E8E8E8", linewidth=0.4, zorder=0)
    for x in range(n_genes):
        ax.axvline(x, color="#E8E8E8", linewidth=0.4, zorder=0)
    x_pos = 0
    for ct in CELLTYPE_ORDER:
        x_pos += len(MARKERS[ct])
        if x_pos < n_genes:
            ax.axvline(x_pos - 0.5, color="#AAAAAA", linewidth=0.7, linestyle="--", zorder=1)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(left=False, bottom=False)

    frac_ticks = [0.0, 0.25, 0.50, 0.75, 1.0]
    size_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#CCCCCC",
               markeredgecolor="#333333", markeredgewidth=0.3,
               markersize=np.sqrt(SIZE_MIN + f * (SIZE_MAX - SIZE_MIN)) * 0.72,
               linestyle="none")
        for f in frac_ticks
    ]
    leg1 = ax_legend.legend(size_handles, [f"{f:.0%}" for f in frac_ticks],
                            title="Fraction\nexpressing", loc="upper left",
                            bbox_to_anchor=(0.05, 0.95), frameon=True, framealpha=0.9,
                            edgecolor="#CCCCCC", labelspacing=1.2, handletextpad=0.6,
                            title_fontsize=7, fontsize=6)
    ax_legend.add_artist(leg1)

    bbox = ax.get_position()
    cbar_ax = fig.add_axes([bbox.x1 + 0.12, bbox.y0 + 0.02, 0.015, 0.25])
    sm = plt.cm.ScalarMappable(cmap=EXPR_CMAP, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_ticks([0, 0.5, 1])
    cbar.set_ticklabels(["Low", "", "High"])
    cbar.ax.tick_params(labelsize=6, length=2)
    cbar.set_label("Normalised\nmean expression", fontsize=7, labelpad=4)
    cbar.outline.set_linewidth(0.4)

    fig.savefig(os.path.join(out_dir, "dotplot_canonical_markers.svg"), dpi=450, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, "dotplot_canonical_markers.pdf"), dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    t_start = time.time()
    print("=" * 60)
    print("Step 09f — Canonical Lineage Marker Dot Plot")
    print("=" * 60)

    stats_path = os.path.join(OUT_DIR, "dotplot_stats.csv")
    lilra4_path = os.path.join(OUT_DIR, "lilra4_stats.csv")
    if os.path.exists(stats_path) and os.path.exists(lilra4_path):
        print(f"Loading cached stats from {stats_path} ...")
        df = pd.read_csv(stats_path)
        present_genes = [g for g in GENE_ORDER if g in df["gene"].unique()]
    else:
        adata, present_genes = load_marker_data()
        print("\nComputing dot-plot statistics ...")
        df = compute_dotplot_stats(adata, present_genes)
        df.to_csv(stats_path, index=False)
        print(f"  Saved  {stats_path}")
        # LILRA4 / pDC marker row — consumed (concatenated) by SF4B. Same schema,
        # computed at celltype_level1 like the main panel.
        if "LILRA4" in adata.var_names:
            compute_dotplot_stats(adata, LILRA4_GENES).to_csv(lilra4_path, index=False)
            print(f"  Saved  {lilra4_path}")
        else:
            print("  [warn] LILRA4 absent from atlas var_names — lilra4_stats.csv not written")

    print("\nRendering dot plot ...")
    plot_dotplot(df, [g for g in GENE_ORDER if g in present_genes], OUT_DIR)
    print(f"\nDone in {time.time() - t_start:.0f}s")
