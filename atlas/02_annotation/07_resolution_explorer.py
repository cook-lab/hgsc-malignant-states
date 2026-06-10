#!/usr/bin/env python3
"""
Atlas 02 — Step 07 (resolution explorer): Leiden resolution sweep on level-1 atlas

PURPOSE
    Subset compartments from the level-1 atlas, recompute neighbours/UMAP on the
    scANVI latent space, sweep Leiden resolutions, and score clusters (silhouette,
    SecA/SecB polarization markers, top markers) to choose the annotation
    resolution. Exploratory annotation aid; emits stats CSVs + an HTML report.

INPUTS
    DATA_ROOT/2026_final_atlas/hgsc_atlas_celltype_level1.h5ad

OUTPUTS
    output_root/02_annotation/07_resolution_explorer/{csvs,figs}/*, *_report.html

MANUSCRIPT PANEL(S)
    Annotation backend (resolution selection); no panel rendered directly.

RUNTIME TIER
    heavy (neighbours + UMAP + multi-resolution Leiden on large subsets).
"""

import base64
import gc
import io
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
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy as scipy_entropy
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    silhouette_score,
)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================================
# PATHS
# ============================================================================

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

# Canonical SecA/SecB polarization signatures (noBCAM 7-gene sets) loaded from the
# shared single source of truth — do NOT inline divergent gene lists.
import yaml as _yaml
with open(Path(__file__).resolve().parents[2] / "shared" / "signatures.yml") as _fh:
    _SIGS = _yaml.safe_load(_fh)
SECA_SIGNATURE = list(_SIGS["SecA"])
SECB_SIGNATURE = list(_SIGS["SecB"])

np.random.seed(SEED)

ATLAS_H5AD  = path("data_root", "2026_final_atlas", "hgsc_atlas_celltype_level1.h5ad")
OUT_DIR     = path("output_root", "02_annotation", "07_resolution_explorer")
CSV_DIR     = os.path.join(OUT_DIR, "csvs")
FIG_DIR     = os.path.join(OUT_DIR, "figs")

for d in [OUT_DIR, CSV_DIR, FIG_DIR]:
    os.makedirs(d, exist_ok=True)

# ============================================================================
# COOK LAB v1.2 STYLE
# ============================================================================

DPI   = 450
ALPHA = 0.6

plt.rcParams.update({
    "font.family":     "Arial",
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "figure.dpi":      DPI,
    "savefig.dpi":     DPI,
    "savefig.bbox":    "tight",
    "svg.fonttype":    "none",
})

KELLY_22 = [
    "#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032",
    "#C2B280", "#848482", "#008856", "#E68FAC", "#0067A5",
    "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
    "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
    "#F2F3F4", "#222222",
]


def _adaptive_pt(n_cells):
    """Adaptive dot size for UMAP scatter plots."""
    return max(0.8, min(12.0, 30000 / max(n_cells, 1)))


# ============================================================================
# LITERATURE-DERIVED GENE SETS (canonical markers, NOT DE-derived)
# ============================================================================

GENE_SETS = {
    "Epithelial": OrderedDict([
        ("SecA (Progenitor)",     SECA_SIGNATURE),
        ("SecB (Differentiated)", SECB_SIGNATURE),
        ("Ciliated",              ["FOXJ1", "CAPS", "TPPP3", "RSPH1", "ZMYND10", "DNAH5"]),
        ("Proliferative",         ["MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CENPF"]),
        ("EMT",                   ["VIM", "CDH2", "SNAI1", "ZEB1", "FN1", "TWIST1"]),
        ("Hypoxia",               ["VEGFA", "BNIP3", "SLC2A1", "LDHA", "CA9", "HIF1A"]),
    ]),
    "Macrophage": OrderedDict([
        ("M1 Pro-inflammatory",    ["TNF", "IL1B", "IL6", "CXCL10", "NOS2", "CD80", "IDO1"]),
        ("M2 Anti-inflammatory",   ["CD163", "MRC1", "MSR1", "TGFB1", "ARG1", "CCL18"]),
        ("TAM Lipid-associated",   ["TREM2", "SPP1", "APOE", "FABP5", "PLIN2", "LIPA"]),
        ("C1Q+ Complement",        ["C1QA", "C1QB", "C1QC", "C3", "VSIG4"]),
        ("Classical Monocyte",     ["CD14", "FCN1", "S100A8", "S100A9", "VCAN"]),
        ("Non-classical Monocyte", ["FCGR3A", "CX3CR1", "CDKN1C", "MS4A7"]),
        ("Proliferative",          ["MKI67", "TOP2A", "STMN1", "TYMS", "CENPF"]),
    ]),
    "T/NK cell": OrderedDict([
        ("CD4 Naive/CM",              ["CCR7", "LEF1", "TCF7", "IL7R", "SELL"]),
        ("CD4 Treg",                  ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF4", "BATF"]),
        ("CD8 Effector",              ["GZMA", "GZMB", "PRF1", "IFNG", "NKG7", "GNLY"]),
        ("CD8 Exhausted/Tumor-react", ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "CXCL13"]),
        ("NK CD56bright",             ["NCAM1", "KLRC1", "XCL1", "XCL2", "AREG"]),
        ("NK CD56dim Cytotoxic",      ["FCGR3A", "FGFBP2", "GZMB", "SPON2", "PRF1"]),
        ("Proliferative",             ["MKI67", "STMN1", "TOP2A", "TYMS", "BIRC5"]),
    ]),
    "Fibroblast": OrderedDict([
        ("iCAF (Inflammatory)",       ["IL6", "CXCL12", "CFD", "PDPN", "C3", "HAS1"]),
        ("myCAF (Myofibroblastic)",   ["ACTA2", "TAGLN", "MMP11", "COL1A1", "THY1", "POSTN"]),
        ("apCAF (Antigen-presenting)", ["CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1"]),
        ("Quiescent/Normal",          ["DCN", "LUM", "GSN", "CFD", "FBLN1"]),
        ("Proliferative",             ["MKI67", "TOP2A", "STMN1", "TYMS"]),
    ]),
    "B cell": OrderedDict([
        ("Naive",           ["IGHD", "TCL1A", "FCER2", "IL4R"]),
        ("Memory",          ["CD27", "IGHG1", "IGHG2", "IGHA1", "AIM2"]),
        ("Germinal Center", ["AICDA", "BCL6", "RGS13", "NEIL1", "MEF2B"]),
        ("Activated",       ["CD83", "NFKBIA", "EGR1", "MYC", "CD69"]),
        ("IFN-stimulated",  ["ISG15", "MX1", "IFI44L", "STAT1", "IFIT3"]),
        ("Proliferative",   ["MKI67", "TOP2A", "STMN1"]),
    ]),
    "Plasma cell": OrderedDict([
        ("Long-lived/Mature", ["XBP1", "IRF4", "PRDM1", "SDC1"]),
        ("Cycling",           ["MKI67", "STMN1", "PCNA", "TYMS"]),
        ("IFN-stimulated",    ["ISG15", "MX1", "IFIT1", "STAT1", "IFI6"]),
    ]),
    "Mesothelial": OrderedDict([
        ("Mesothelial identity", ["MSLN", "WT1", "CALB2", "UPK3B", "KRT19"]),
        ("EMT/Activated",        ["VIM", "SNAI1", "CDH2", "FN1", "MMP2", "LOX"]),
        ("Complement",           ["C3", "CFB", "CFD", "SERPING1", "C1S"]),
    ]),
    "Smooth muscle": OrderedDict([
        ("Contractile", ["MYH11", "ACTA2", "TAGLN", "CNN1", "LMOD1"]),
        ("Synthetic",   ["COL1A2", "COL3A1", "FN1", "SPARC", "IGFBP7"]),
    ]),
    "Pericyte": OrderedDict([
        ("Pericyte identity",  ["RGS5", "PDGFRB", "NOTCH3", "ACTA2", "MCAM"]),
        ("Basement membrane",  ["COL4A1", "COL4A2", "LAMB1", "FN1", "SPARC"]),
    ]),
    "Endothelial": OrderedDict([
        ("Arterial",        ["EFNB2", "GJA5", "GJA4", "SEMA3G", "HEY1"]),
        ("Venous",          ["ACKR1", "VWF", "NR2F2", "PLVAP"]),
        ("Lymphatic",       ["PROX1", "LYVE1", "CCL21", "FLT4", "PDPN"]),
        ("Tip/Angiogenic",  ["ESM1", "APLN", "PGF", "KDR", "CXCR4"]),
    ]),
    "DC": OrderedDict([
        ("cDC1",              ["CLEC9A", "XCR1", "BATF3", "IRF8", "IDO1"]),
        ("cDC2",              ["CD1C", "CLEC10A", "FCER1A", "IRF4"]),
        ("pDC",               ["LILRA4", "IRF7", "TCF4", "GZMB", "BCL11A"]),
        ("Mature/Migratory",  ["CCR7", "LAMP3", "FSCN1", "CD83", "BIRC3"]),
    ]),
    "Neutrophil": OrderedDict([
        ("Classical/Mature", ["S100A8", "S100A9", "CSF3R", "FCGR3B", "CXCR2"]),
        ("IFN-stimulated",   ["ISG15", "IFIT3", "MX1", "RSAD2", "STAT1"]),
        ("TAN",              ["VEGFA", "MMP9", "ARG1", "CCL2", "CXCL8"]),
    ]),
    "Mast cell": OrderedDict([
        ("Mast identity",         ["TPSAB1", "TPSB2", "CPA3", "KIT", "FCER1A", "HPGDS"]),
        ("Activated/Degranulating", ["IL1B", "TNF", "CCL3", "CCL4", "NFKBIA"]),
        ("Proliferative",          ["MKI67", "TOP2A", "STMN1"]),
    ]),
}

# ============================================================================
# COMPARTMENT CONFIGURATIONS
# ============================================================================

# Adaptive resolution ranges by size tier
RES_SMALL  = [0.05, 0.1, 0.15, 0.2, 0.3]    # <15K cells
RES_MEDIUM = [0.1, 0.2, 0.3, 0.4, 0.5]       # 15K-100K
RES_LARGE  = [0.1, 0.2, 0.3, 0.5, 0.8]       # >100K

# Maps celltype_level1 value → config
COMPARTMENTS = OrderedDict([
    ("mastcell",     {"level1": "Mast cell",     "resolutions": RES_SMALL,  "purity": ["Mast"]}),
    ("neutrophil",   {"level1": "Neutrophil",    "resolutions": RES_SMALL,  "purity": ["Neutrophil"]}),
    ("dc",           {"level1": "DC",            "resolutions": RES_SMALL,  "purity": ["DC"]}),
    ("pericyte",     {"level1": "Pericyte",      "resolutions": RES_MEDIUM, "purity": ["Pericyte"]}),
    ("smoothmuscle", {"level1": "Smooth muscle", "resolutions": RES_MEDIUM, "purity": ["Smooth_Muscle"]}),
    ("bcell",        {"level1": "B cell",        "resolutions": RES_MEDIUM, "purity": ["B_cell"]}),
    ("plasmacell",   {"level1": "Plasma cell",   "resolutions": RES_MEDIUM, "purity": ["Plasma_cell"]}),
    ("mesothelial",  {"level1": "Mesothelial",   "resolutions": RES_MEDIUM, "purity": ["Mesothelial"]}),
    ("endothelial",  {"level1": "Endothelial",   "resolutions": RES_MEDIUM, "purity": ["Endothelial"]}),
    ("fibroblast",   {"level1": "Fibroblast",    "resolutions": RES_LARGE,  "purity": ["Fibroblast"]}),
    ("macrophage",   {"level1": "Macrophage",    "resolutions": RES_LARGE,  "purity": ["Macrophage"]}),
    ("tnkcell",      {"level1": "T/NK cell",     "resolutions": [0.1, 0.3, 0.5, 0.7],  "purity": ["T_cell", "NK_cell"]}),
    ("epithelial",   {"level1": "Epithelial",    "resolutions": [0.1, 0.3, 0.5, 0.7],  "purity": ["Epithelial"],
                      "split_ciliated": True}),
])

# Run order: smallest → largest
RUN_ORDER = list(COMPARTMENTS.keys())


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _svg_to_b64(fig):
    """Convert matplotlib figure to base64-encoded SVG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", dpi=DPI, bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return encoded


def _fig_tag(b64, width="100%"):
    """HTML img tag from base64 SVG."""
    return f'<img src="data:image/svg+xml;base64,{b64}" style="width:{width}; max-width:1200px;" />'


def _save_svg(fig, name):
    """Save figure as SVG to FIG_DIR and return base64."""
    path = os.path.join(FIG_DIR, f"{name}.svg")
    fig.savefig(path, format="svg", dpi=DPI, bbox_inches="tight")
    b64 = _svg_to_b64(fig)
    return b64


# ---------- UMAP ----------

def plot_umap_cat(adata, color_key, title, pt_size=None, palette=None):
    """Categorical UMAP. Returns base64 SVG."""
    if pt_size is None:
        pt_size = _adaptive_pt(len(adata))
    fig, ax = plt.subplots(figsize=(6, 5))
    cats = adata.obs[color_key].cat.categories if hasattr(adata.obs[color_key], "cat") else sorted(adata.obs[color_key].unique())
    n_cats = len(cats)
    if palette is None:
        palette_list = (KELLY_22 * ((n_cats // len(KELLY_22)) + 1))[:n_cats]
        palette = {str(c): palette_list[i] for i, c in enumerate(cats)}
    for i, cat in enumerate(cats):
        mask = adata.obs[color_key] == cat
        n = mask.sum()
        coords = adata.obsm["X_umap"][mask.values]
        color = palette.get(str(cat), "#999999")
        ax.scatter(coords[:, 0], coords[:, 1], s=pt_size, alpha=ALPHA,
                   c=color, label=f"{cat} ({n:,})", rasterized=True, edgecolors="none")
    ax.set_title(title, fontsize=9, pad=8)
    ax.set_xlabel("UMAP1", fontsize=7)
    ax.set_ylabel("UMAP2", fontsize=7)
    ax.tick_params(labelsize=6)
    # Legend outside
    if n_cats <= 20:
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=5,
                  markerscale=3, frameon=False, handletextpad=0.3)
    else:
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=4,
                  markerscale=2, frameon=False, ncol=2, handletextpad=0.2)
    fig.tight_layout()
    return fig


def plot_umap_cont(adata, gene, title=None, pt_size=None):
    """Continuous gene expression UMAP. Returns base64 SVG."""
    if pt_size is None:
        pt_size = _adaptive_pt(len(adata))
    if title is None:
        title = gene
    fig, ax = plt.subplots(figsize=(5, 4.5))
    if gene in adata.var_names:
        gi = list(adata.var_names).index(gene)
        expr = adata.X[:, gi]
        if hasattr(expr, "toarray"):
            expr = expr.toarray().ravel()
        else:
            expr = np.asarray(expr).ravel()
    else:
        expr = np.zeros(len(adata))
    order = np.argsort(expr)
    sc_plot = ax.scatter(
        adata.obsm["X_umap"][order, 0], adata.obsm["X_umap"][order, 1],
        c=expr[order], cmap="magma", s=pt_size, alpha=ALPHA,
        rasterized=True, edgecolors="none", vmin=0,
    )
    ax.set_title(title, fontsize=9, pad=6)
    ax.set_xlabel("UMAP1", fontsize=7)
    ax.set_ylabel("UMAP2", fontsize=7)
    ax.tick_params(labelsize=6)
    cbar = fig.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
    cbar.ax.tick_params(labelsize=5)
    fig.tight_layout()
    return fig


# ---------- QC VIOLINS ----------

def plot_qc_violins(adata, cluster_key, metrics, title_prefix=""):
    """Violin plots of QC metrics per cluster. Returns base64 SVG."""
    n_metrics = len(metrics)
    fig, axes = plt.subplots(1, n_metrics, figsize=(3.5 * n_metrics, 4))
    if n_metrics == 1:
        axes = [axes]
    for ax, (col, label) in zip(axes, metrics):
        cats = sorted(adata.obs[cluster_key].unique(), key=lambda x: int(x) if x.isdigit() else x)
        data_list = []
        positions = []
        for i, cat in enumerate(cats):
            mask = adata.obs[cluster_key] == cat
            vals = adata.obs.loc[mask, col].dropna().values.astype(float)
            data_list.append(vals)
            positions.append(i)
        if all(len(d) > 0 for d in data_list):
            vp = ax.violinplot(data_list, positions=positions, showmedians=True, widths=0.8)
            for body in vp.get("bodies", []):
                body.set_alpha(0.6)
        ax.set_xticks(positions)
        ax.set_xticklabels([str(c) for c in cats], fontsize=6, rotation=45, ha="right")
        ax.set_title(f"{title_prefix}{label}", fontsize=8, pad=6)
        ax.set_xlabel("Cluster", fontsize=7)
        ax.tick_params(labelsize=6)
    fig.tight_layout()
    return fig


# ---------- STACKED BARS ----------

def plot_stacked_bar(adata, cluster_key, meta_key, title):
    """Stacked bar chart of metadata composition per cluster. Returns base64 SVG."""
    ct = pd.crosstab(adata.obs[cluster_key], adata.obs[meta_key], normalize="index")
    cats = sorted(ct.index, key=lambda x: int(x) if str(x).isdigit() else x)
    ct = ct.loc[cats]
    n_cols = len(ct.columns)
    colors = (KELLY_22 * ((n_cols // len(KELLY_22)) + 1))[:n_cols]
    fig, ax = plt.subplots(figsize=(max(6, len(cats) * 0.6), 4))
    ct.plot.bar(stacked=True, ax=ax, color=colors, width=0.85, edgecolor="none")
    ax.set_title(title, fontsize=9, pad=8)
    ax.set_xlabel("Cluster", fontsize=7)
    ax.set_ylabel("Proportion", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=5,
              frameon=False, handletextpad=0.3)
    ax.set_ylim(0, 1)
    fig.tight_layout()
    return fig


# ---------- DOTPLOT ----------

def plot_marker_dotplot(adata, gene_sets, cluster_key, title=""):
    """Dotplot of literature marker genes vs clusters. Returns base64 SVG."""
    # Build ordered gene list + group annotations
    dot_genes = []
    dot_positions = []
    dot_labels = []
    start = 0
    var_set = set(adata.var_names)
    for group_name, genes in gene_sets.items():
        valid = [g for g in genes if g in var_set]
        if not valid:
            continue
        dot_positions.append((start, start + len(valid) - 1))
        dot_labels.append(group_name)
        dot_genes.extend(valid)
        start += len(valid)

    if not dot_genes:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No marker genes found in dataset", ha="center", va="center")
        ax.set_axis_off()
        return fig

    # Order clusters numerically
    cats = sorted(adata.obs[cluster_key].unique(), key=lambda x: int(x) if str(x).isdigit() else x)
    cluster_labels = [str(c) for c in cats]

    try:
        dp = sc.pl.dotplot(
            adata,
            var_names=dot_genes,
            groupby=cluster_key,
            var_group_positions=dot_positions,
            var_group_labels=dot_labels,
            categories_order=cluster_labels,
            vmax=5,
            figsize=(max(12, len(dot_genes) * 0.35), max(3, len(cats) * 0.35)),
            show=False,
            return_fig=True,
        )
        dp.style(cmap="Reds", dot_edge_color="none", dot_edge_lw=0, size_exponent=1.5)
        dp.make_figure()
        fig = dp.fig
        if title:
            fig.suptitle(title, fontsize=10, y=1.02)
    except Exception as e:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, f"Dotplot error: {e}", ha="center", va="center", fontsize=8)
        ax.set_axis_off()
    return fig


# ---------- STATISTICS ----------

def compute_resolution_stats(adata, resolutions, max_cells_silhouette=50000):
    """Compute clustering quality metrics for each resolution.
    Uses X_scanvi (scANVI latent space) for all distance-based metrics."""
    latent = adata.obsm["X_scanvi"]
    results = []

    for res in resolutions:
        key = f"leiden_r{res}"
        labels = adata.obs[key].values.astype(str)
        unique_labels = np.unique(labels)
        n_clusters = len(unique_labels)

        # Cluster sizes
        sizes = pd.Series(labels).value_counts().values
        min_size = int(sizes.min())
        max_size = int(sizes.max())
        size_ent = float(scipy_entropy(sizes / sizes.sum()))

        row = {
            "resolution": res,
            "n_clusters": n_clusters,
            "min_cluster_size": min_size,
            "max_cluster_size": max_size,
            "size_entropy": round(size_ent, 4),
        }

        # Metrics that need >1 cluster
        if n_clusters > 1:
            # Subsample for silhouette if needed
            n = len(adata)
            if n > max_cells_silhouette:
                idx = np.random.choice(n, max_cells_silhouette, replace=False)
                latent_sub = latent[idx]
                labels_sub = labels[idx]
            else:
                latent_sub = latent
                labels_sub = labels

            row["silhouette"] = round(float(silhouette_score(latent_sub, labels_sub)), 4)
            row["calinski_harabasz"] = round(float(calinski_harabasz_score(latent, labels)), 2)
            row["davies_bouldin"] = round(float(davies_bouldin_score(latent, labels)), 4)

            # Batch mixing: ARI between leiden and study
            if "study" in adata.obs.columns:
                row["study_ari"] = round(float(
                    adjusted_rand_score(labels, adata.obs["study"].values.astype(str))
                ), 4)
        else:
            row["silhouette"] = None
            row["calinski_harabasz"] = None
            row["davies_bouldin"] = None
            row["study_ari"] = None

        results.append(row)

    return pd.DataFrame(results)


def plot_resolution_stats(stats_df):
    """Line plots of resolution comparison metrics. Returns base64 SVG."""
    metrics = [
        ("silhouette", "Silhouette Score (↑)", True),
        ("calinski_harabasz", "Calinski-Harabasz (↑)", True),
        ("davies_bouldin", "Davies-Bouldin (↓)", False),
        ("n_clusters", "N Clusters", None),
        ("size_entropy", "Size Entropy (↑)", True),
    ]
    # Only plot metrics that exist
    valid = [(m, l, d) for m, l, d in metrics if m in stats_df.columns and stats_df[m].notna().any()]
    n = len(valid)
    if n == 0:
        fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, "No metrics to plot", ha="center", va="center")
        return fig

    fig, axes = plt.subplots(1, n, figsize=(3.5 * n, 3))
    if n == 1:
        axes = [axes]

    for ax, (metric, label, higher_better) in zip(axes, valid):
        vals = stats_df[metric].values
        res_vals = stats_df["resolution"].values
        valid_mask = pd.notna(vals)
        ax.plot(res_vals[valid_mask], vals[valid_mask], "o-", color="#0067A5", markersize=5)
        ax.set_xlabel("Resolution", fontsize=7)
        ax.set_title(label, fontsize=8, pad=6)
        ax.tick_params(labelsize=6)
        if higher_better is not None:
            best_idx = np.nanargmax(vals) if higher_better else np.nanargmin(vals)
            if valid_mask[best_idx]:
                ax.axvline(res_vals[best_idx], color="#BE0032", alpha=0.3, ls="--", lw=1)

    fig.tight_layout()
    return fig


# ---------- RECOMMENDATION ----------

def generate_recommendation(stats_df):
    """Auto-generate a resolution recommendation."""
    lines = []
    if stats_df is None or len(stats_df) == 0:
        return "No statistics available."

    # Find best by silhouette
    if "silhouette" in stats_df.columns and stats_df["silhouette"].notna().any():
        best_sil = stats_df.loc[stats_df["silhouette"].idxmax()]
        lines.append(f"**Best silhouette score**: r={best_sil['resolution']} "
                     f"(score={best_sil['silhouette']:.4f}, {int(best_sil['n_clusters'])} clusters)")

    # Find best by CH
    if "calinski_harabasz" in stats_df.columns and stats_df["calinski_harabasz"].notna().any():
        best_ch = stats_df.loc[stats_df["calinski_harabasz"].idxmax()]
        lines.append(f"**Best Calinski-Harabasz**: r={best_ch['resolution']} "
                     f"(score={best_ch['calinski_harabasz']:.1f})")

    # Find best by DB (lower is better)
    if "davies_bouldin" in stats_df.columns and stats_df["davies_bouldin"].notna().any():
        best_db = stats_df.loc[stats_df["davies_bouldin"].idxmin()]
        lines.append(f"**Best Davies-Bouldin**: r={best_db['resolution']} "
                     f"(score={best_db['davies_bouldin']:.4f})")

    # Warn about small clusters
    for _, row in stats_df.iterrows():
        if row["min_cluster_size"] < 200:
            lines.append(f"⚠️ r={row['resolution']}: smallest cluster has only "
                         f"{int(row['min_cluster_size'])} cells")

    # Consensus
    if "silhouette" in stats_df.columns and stats_df["silhouette"].notna().any():
        # Rank each metric, average ranks
        rank_df = pd.DataFrame({"resolution": stats_df["resolution"]})
        if stats_df["silhouette"].notna().any():
            rank_df["sil_rank"] = stats_df["silhouette"].rank(ascending=False)
        if "calinski_harabasz" in stats_df.columns and stats_df["calinski_harabasz"].notna().any():
            rank_df["ch_rank"] = stats_df["calinski_harabasz"].rank(ascending=False)
        if "davies_bouldin" in stats_df.columns and stats_df["davies_bouldin"].notna().any():
            rank_df["db_rank"] = stats_df["davies_bouldin"].rank(ascending=True)
        rank_cols = [c for c in rank_df.columns if c.endswith("_rank")]
        if rank_cols:
            rank_df["avg_rank"] = rank_df[rank_cols].mean(axis=1)
            best = rank_df.loc[rank_df["avg_rank"].idxmin()]
            lines.append(f"\n**Consensus recommendation**: r={best['resolution']} "
                         f"(avg rank {best['avg_rank']:.1f} across metrics)")

    if not lines:
        return "Insufficient data for recommendation."

    lines.append("\n*This is a statistical suggestion. Expert biological review should "
                 "determine the final resolution based on marker separation and biological interpretability.*")
    return "\n".join(lines)


# ============================================================================
# MAIN PIPELINE: PROCESS ONE COMPARTMENT
# ============================================================================

def process_compartment(atlas, comp_key, config):
    """Run the full resolution exploration pipeline for one compartment."""
    t0 = time.time()
    level1_name = config["level1"]
    resolutions = config["resolutions"]
    purity_types = config["purity"]
    split_ciliated = config.get("split_ciliated", False)

    prefix = comp_key
    print(f"\n{'='*70}")
    print(f"  Processing: {level1_name} ({comp_key})")
    print(f"{'='*70}\n")

    # ------------------------------------------------------------------
    # 1. Subset
    # ------------------------------------------------------------------
    print("1. Subsetting from atlas...", flush=True)
    mask = atlas.obs["celltype_level1"] == level1_name
    adata = atlas[mask].copy()
    print(f"   {len(adata):,} cells x {adata.n_vars:,} genes")

    # Handle epithelial ciliated split
    ciliated_adata = None
    if split_ciliated:
        print("   Splitting ciliated cells...", flush=True)
        cil_markers = ["FOXJ1", "CAPS", "TPPP3", "RSPH1"]
        valid_markers = [g for g in cil_markers if g in adata.var_names]
        if valid_markers and "celltype_pred" in adata.obs.columns:
            cil_mask = adata.obs["celltype_pred"] == "Ciliated"
            n_cil = cil_mask.sum()
            if n_cil > 0:
                ciliated_adata = adata[cil_mask].copy()
                adata = adata[~cil_mask].copy()
                print(f"   Ciliated: {n_cil:,} cells (separated)")
                print(f"   Secretory: {len(adata):,} cells (for clustering)")
            else:
                print("   No ciliated cells found by celltype_pred, using all cells")
        else:
            print("   Cannot split ciliated (missing markers/pred), using all cells")

    n_cells = len(adata)

    # ------------------------------------------------------------------
    # 2. Neighbors + UMAP from scANVI embedding
    # ------------------------------------------------------------------
    print("\n2. Computing neighbors → UMAP from X_scanvi...", flush=True)
    if "X_scanvi" not in adata.obsm:
        print("   ERROR: X_scanvi not found in obsm. Available:", list(adata.obsm.keys()))
        raise ValueError("X_scanvi embedding not found — cannot proceed without batch-corrected latent space")
    n_scanvi = adata.obsm["X_scanvi"].shape[1]
    print(f"   X_scanvi: {n_scanvi}-dimensional latent space", flush=True)
    sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=15)
    print("   Neighbors done (from X_scanvi)", flush=True)
    sc.tl.umap(adata, random_state=SEED)
    print("   UMAP done", flush=True)

    # ------------------------------------------------------------------
    # 3. Multi-resolution Leiden
    # ------------------------------------------------------------------
    print(f"\n3. Running Leiden at {len(resolutions)} resolutions...", flush=True)
    for res in resolutions:
        key = f"leiden_r{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=SEED)
        n_cl = adata.obs[key].nunique()
        print(f"   r={res}: {n_cl} clusters", flush=True)

    # ------------------------------------------------------------------
    # 4. MALAT1/NEAT1 low-quality score
    # ------------------------------------------------------------------
    print("\n4. Computing lq_score (MALAT1 + NEAT1)...", flush=True)
    malat1 = np.zeros(n_cells)
    neat1  = np.zeros(n_cells)
    if "MALAT1" in adata.var_names:
        gi = list(adata.var_names).index("MALAT1")
        v = adata.X[:, gi]
        malat1 = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()
    if "NEAT1" in adata.var_names:
        gi = list(adata.var_names).index("NEAT1")
        v = adata.X[:, gi]
        neat1 = v.toarray().ravel() if hasattr(v, "toarray") else np.asarray(v).ravel()

    total = adata.obs["total_counts"].values.astype(float)
    total = np.where(total > 0, total, 1)
    adata.obs["lq_score"] = (malat1 + neat1) / total
    lq_95 = float(np.percentile(adata.obs["lq_score"], 95))
    lq_median = float(np.median(adata.obs["lq_score"]))
    print(f"   Median lq_score: {lq_median:.4f}")
    print(f"   95th percentile: {lq_95:.4f}")

    # ------------------------------------------------------------------
    # 5. Resolution statistics
    # ------------------------------------------------------------------
    print("\n5. Computing resolution statistics...", flush=True)
    stats_df = compute_resolution_stats(adata, resolutions)
    stats_csv = os.path.join(CSV_DIR, f"{prefix}_resolution_stats.csv")
    stats_df.to_csv(stats_csv, index=False)
    print(f"   Saved: {stats_csv}")
    print(stats_df.to_string(index=False))

    recommendation = generate_recommendation(stats_df)

    # ------------------------------------------------------------------
    # 6. Generate all figures
    # ------------------------------------------------------------------
    print("\n6. Generating figures...", flush=True)
    html_sections = []
    pt = _adaptive_pt(n_cells)

    # Get gene sets for this compartment
    gene_set_key = level1_name
    if gene_set_key not in GENE_SETS:
        # Try variations
        for k in GENE_SETS:
            if k.lower().replace(" ", "").replace("/", "") in gene_set_key.lower().replace(" ", "").replace("/", ""):
                gene_set_key = k
                break
    gene_sets = GENE_SETS.get(gene_set_key, {})

    # QC metrics to violin-plot
    qc_metrics = []
    if "n_genes_by_counts" in adata.obs.columns:
        qc_metrics.append(("n_genes_by_counts", "Genes Detected"))
    if "total_counts" in adata.obs.columns:
        qc_metrics.append(("total_counts", "UMI Counts"))
    if "doublet_score_scrublet" in adata.obs.columns:
        qc_metrics.append(("doublet_score_scrublet", "Doublet Score"))
    qc_metrics.append(("lq_score", "LQ Score (MALAT1+NEAT1)"))

    # --- Overview UMAP (before clustering) ---
    print("   Overview UMAP...", flush=True)
    overview_res = resolutions[min(2, len(resolutions) - 1)]
    fig_overview = plot_umap_cat(adata, "leiden_r" + str(overview_res),
                                  f"{level1_name} — Fresh UMAP", pt_size=pt)
    b64_overview = _save_svg(fig_overview, f"{prefix}_overview_umap")

    # Compute CellAssign purity per cluster for each resolution
    has_celltype_pred = "celltype_pred" in adata.obs.columns

    # --- Per-resolution analysis ---
    resolution_htmls = []
    for res_idx, res in enumerate(resolutions):
        key = f"leiden_r{res}"
        n_cl = adata.obs[key].nunique()
        print(f"\n   --- Resolution {res} ({n_cl} clusters) ---", flush=True)
        res_figs = {}

        # 6a. UMAP colored by leiden
        print(f"   [{res}] UMAP leiden...", flush=True)
        fig = plot_umap_cat(adata, key, f"{level1_name} — Leiden r={res}", pt_size=pt)
        res_figs["umap_leiden"] = _save_svg(fig, f"{prefix}_r{res}_umap_leiden")

        # 6b. Marker dotplot
        if gene_sets:
            print(f"   [{res}] Dotplot...", flush=True)
            fig = plot_marker_dotplot(adata, gene_sets, key,
                                      title=f"{level1_name} — Literature Markers (r={res})")
            res_figs["dotplot"] = _save_svg(fig, f"{prefix}_r{res}_dotplot")

        # 6b2. Gene UMAPs (key markers — first gene from each group)
        print(f"   [{res}] Gene UMAPs...", flush=True)
        gene_umap_b64s = []
        key_genes = []
        for group_name, genes in gene_sets.items():
            for g in genes:
                if g in adata.var_names:
                    key_genes.append((group_name, g))
                    break
        # Plot up to 8 representative genes in a 2x4 grid
        genes_to_plot = key_genes[:8]
        if genes_to_plot:
            n_genes = len(genes_to_plot)
            ncols = min(4, n_genes)
            nrows = (n_genes + ncols - 1) // ncols
            fig_genes, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4 * nrows))
            if nrows == 1 and ncols == 1:
                axes = np.array([[axes]])
            elif nrows == 1:
                axes = axes[np.newaxis, :]
            elif ncols == 1:
                axes = axes[:, np.newaxis]
            for idx, (group, gene) in enumerate(genes_to_plot):
                r, c = divmod(idx, ncols)
                ax = axes[r, c]
                gi = list(adata.var_names).index(gene)
                expr = adata.X[:, gi]
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().ravel()
                else:
                    expr = np.asarray(expr).ravel()
                order = np.argsort(expr)
                ax.scatter(adata.obsm["X_umap"][order, 0], adata.obsm["X_umap"][order, 1],
                          c=expr[order], cmap="magma", s=pt * 0.5, alpha=ALPHA,
                          rasterized=True, edgecolors="none", vmin=0)
                ax.set_title(f"{gene} ({group})", fontsize=7, pad=4)
                ax.tick_params(labelsize=5)
            # Hide empty axes
            for idx in range(len(genes_to_plot), nrows * ncols):
                r, c = divmod(idx, ncols)
                axes[r, c].set_visible(False)
            fig_genes.suptitle(f"{level1_name} — Key Markers (r={res})", fontsize=9, y=1.01)
            fig_genes.tight_layout()
            res_figs["gene_umaps"] = _save_svg(fig_genes, f"{prefix}_r{res}_gene_umaps")

        # 6c. Top 50 DEGs
        print(f"   [{res}] Wilcoxon DEGs...", flush=True)
        try:
            sc.tl.rank_genes_groups(adata, groupby=key, method="wilcoxon",
                                     n_genes=50, use_raw=False)
            # Extract markers
            marker_rows = []
            for cl in sorted(adata.obs[key].unique(), key=lambda x: int(x) if x.isdigit() else x):
                try:
                    names = adata.uns["rank_genes_groups"]["names"][cl][:50]
                    scores = adata.uns["rank_genes_groups"]["scores"][cl][:50]
                    logfc = adata.uns["rank_genes_groups"]["logfoldchanges"][cl][:50]
                    pvals = adata.uns["rank_genes_groups"]["pvals_adj"][cl][:50]
                    for rank, (n, s, l, p) in enumerate(zip(names, scores, logfc, pvals), 1):
                        marker_rows.append({
                            "cluster": cl, "rank": rank, "gene": n,
                            "score": round(float(s), 3),
                            "logFC": round(float(l), 3),
                            "pval_adj": float(p),
                        })
                except (KeyError, IndexError):
                    pass

            markers_df = pd.DataFrame(marker_rows)
            csv_path = os.path.join(CSV_DIR, f"{prefix}_r{res}_markers_top50.csv")
            markers_df.to_csv(csv_path, index=False)
            res_figs["degs_csv"] = csv_path

            # Build top 10 table for HTML
            deg_table_html = _deg_table_html(markers_df, key, adata)

        except Exception as e:
            print(f"   WARNING: DEG computation failed: {e}")
            deg_table_html = f"<p>DEG computation failed: {e}</p>"
            markers_df = pd.DataFrame()

        res_figs["deg_table"] = deg_table_html

        # 6d. QC violins
        print(f"   [{res}] QC violins...", flush=True)
        fig = plot_qc_violins(adata, key, qc_metrics, title_prefix=f"r={res} — ")
        res_figs["qc_violins"] = _save_svg(fig, f"{prefix}_r{res}_qc_violins")

        # CellAssign purity table
        purity_html = ""
        if has_celltype_pred:
            purity_rows = []
            for cl in sorted(adata.obs[key].unique(), key=lambda x: int(x) if x.isdigit() else x):
                cl_mask = adata.obs[key] == cl
                n_cl = cl_mask.sum()
                pred = adata.obs.loc[cl_mask, "celltype_pred"]
                purity_pct = 0
                for pt_name in purity_types:
                    purity_pct += (pred == pt_name).sum()
                purity_pct = 100.0 * purity_pct / max(n_cl, 1)
                # Top predicted types
                top_pred = pred.value_counts().head(3)
                top_str = ", ".join([f"{v}:{c}" for v, c in top_pred.items()])
                purity_rows.append({
                    "Cluster": cl, "N": f"{n_cl:,}",
                    "Purity %": f"{purity_pct:.1f}%",
                    "Top predictions": top_str,
                    "Med genes": f"{adata.obs.loc[cl_mask, 'n_genes_by_counts'].median():.0f}",
                    "Med lq": f"{adata.obs.loc[cl_mask, 'lq_score'].median():.4f}",
                })
            purity_df = pd.DataFrame(purity_rows)
            purity_html = purity_df.to_html(index=False, classes="styled-table", border=0)

        res_figs["purity_table"] = purity_html

        # 6e. Stacked bars
        print(f"   [{res}] Stacked bars...", flush=True)
        bar_figs = {}
        for meta_key, meta_label in [("treatment_status", "Treatment Status"),
                                       ("anatomic_site", "Anatomic Site"),
                                       ("study", "Study")]:
            if meta_key in adata.obs.columns:
                fig = plot_stacked_bar(adata, key, meta_key,
                                        f"{level1_name} r={res} — {meta_label}")
                bar_figs[meta_key] = _save_svg(fig, f"{prefix}_r{res}_bar_{meta_key}")

        res_figs["stacked_bars"] = bar_figs

        resolution_htmls.append((res, n_cl, res_figs))

    # ------------------------------------------------------------------
    # 7. Stats plots
    # ------------------------------------------------------------------
    print("\n   Resolution stats plot...", flush=True)
    fig_stats = plot_resolution_stats(stats_df)
    b64_stats = _save_svg(fig_stats, f"{prefix}_resolution_stats")

    # ------------------------------------------------------------------
    # 8. Build HTML
    # ------------------------------------------------------------------
    print("\n7. Building HTML report...", flush=True)
    html = _build_html(
        comp_key=comp_key,
        level1_name=level1_name,
        n_cells=n_cells,
        resolutions=resolutions,
        stats_df=stats_df,
        recommendation=recommendation,
        b64_overview=b64_overview,
        b64_stats=b64_stats,
        resolution_htmls=resolution_htmls,
        lq_95=lq_95,
        lq_median=lq_median,
        ciliated_adata=ciliated_adata,
    )

    html_path = os.path.join(OUT_DIR, f"{prefix}_resolution_report.html")
    with open(html_path, "w") as f:
        f.write(html)
    print(f"   Saved: {html_path}")

    elapsed = time.time() - t0
    print(f"\n   ✓ {level1_name} complete in {elapsed/60:.1f} min")

    # Cleanup
    del adata
    if ciliated_adata is not None:
        del ciliated_adata
    gc.collect()

    return html_path


def _deg_table_html(markers_df, cluster_key, adata):
    """Build HTML table of top 10 DEGs per cluster."""
    if markers_df.empty:
        return "<p>No DEGs computed.</p>"

    rows = []
    clusters = sorted(markers_df["cluster"].unique(),
                       key=lambda x: int(x) if str(x).isdigit() else x)
    for cl in clusters:
        cl_df = markers_df[markers_df["cluster"] == cl].head(10)
        genes_str = ", ".join(cl_df["gene"].tolist())
        n_cl = (adata.obs[cluster_key] == cl).sum()
        rows.append(f"<tr><td><b>{cl}</b></td><td>{n_cl:,}</td><td>{genes_str}</td></tr>")

    return (
        '<table class="styled-table" border="0">'
        "<tr><th>Cluster</th><th>N cells</th><th>Top 10 DEGs</th></tr>"
        + "".join(rows)
        + "</table>"
    )


# ============================================================================
# HTML BUILDER
# ============================================================================

def _build_html(comp_key, level1_name, n_cells, resolutions, stats_df,
                recommendation, b64_overview, b64_stats, resolution_htmls,
                lq_95, lq_median, ciliated_adata=None):
    """Build self-contained HTML report."""

    css = """
    <style>
        body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; background: #fafafa; }
        h1 { color: #333; border-bottom: 3px solid #0067A5; padding-bottom: 10px; }
        h2 { color: #0067A5; margin-top: 40px; border-bottom: 1px solid #ddd; padding-bottom: 5px; }
        h3 { color: #604E97; margin-top: 25px; }
        .overview-box { background: #f0f4f8; padding: 15px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #0067A5; }
        .recommendation-box { background: #fff8e1; padding: 15px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #F6A600; }
        .warning { color: #BE0032; font-weight: bold; }
        .styled-table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }
        .styled-table th { background: #0067A5; color: white; padding: 8px 12px; text-align: left; }
        .styled-table td { padding: 6px 12px; border-bottom: 1px solid #eee; }
        .styled-table tr:nth-child(even) { background: #f8f9fa; }
        .fig-container { margin: 15px 0; text-align: center; }
        .fig-container img { border: 1px solid #eee; border-radius: 4px; }
        .section-divider { border-top: 2px solid #0067A5; margin: 40px 0; }
        .resolution-header { background: #e8f0fe; padding: 12px; border-radius: 6px; margin: 10px 0; }
        .stats-highlight { font-size: 14px; font-weight: bold; color: #0067A5; }
        .lq-box { background: #fce4ec; padding: 15px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #BE0032; }
        .grid-2col { display: grid; grid-template-columns: 1fr 1fr; gap: 15px; }
        @media (max-width: 900px) { .grid-2col { grid-template-columns: 1fr; } }
    </style>
    """

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    ciliated_note = ""
    if ciliated_adata is not None:
        ciliated_note = f"<p>Ciliated cells separated: <b>{len(ciliated_adata):,}</b> (not included in clustering)</p>"

    parts = [f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>{level1_name} — Resolution Explorer</title>
{css}
</head><body>
<h1>{level1_name} — Multi-Resolution Clustering Explorer</h1>

<div class="overview-box">
<p><b>Compartment:</b> {level1_name} | <b>Cells for clustering:</b> {n_cells:,} |
<b>Resolutions tested:</b> {', '.join(str(r) for r in resolutions)}</p>
{ciliated_note}
<p><b>Source:</b> hgsc_atlas_celltype_level1.h5ad (X_scanvi → neighbors → UMAP)</p>
<p><b>Generated:</b> {timestamp}</p>
</div>

<h2>1. Overview UMAP</h2>
<div class="fig-container">{_fig_tag(b64_overview)}</div>
"""]

    # --- Resolution comparison summary ---
    parts.append('<div class="section-divider"></div>')
    parts.append("<h2>2. Resolution Comparison Summary</h2>")

    # Stats table
    parts.append(stats_df.to_html(index=False, classes="styled-table", border=0,
                                    float_format=lambda x: f"{x:.4f}" if isinstance(x, float) else str(x)))

    parts.append(f'<div class="fig-container">{_fig_tag(b64_stats)}</div>')

    parts.append('<div class="recommendation-box">')
    parts.append("<h3>Statistical Recommendation</h3>")
    # Convert markdown-like bold to HTML
    rec_html = recommendation.replace("**", "<b>").replace("**", "</b>")
    # Fix: proper bold tags
    import re
    rec_html = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', recommendation)
    rec_html = rec_html.replace("\n", "<br>")
    parts.append(f"<p>{rec_html}</p></div>")

    # --- Per-resolution sections ---
    for sec_idx, (res, n_cl, figs) in enumerate(resolution_htmls, start=3):
        parts.append('<div class="section-divider"></div>')
        parts.append(f'<h2>{sec_idx}. Resolution {res} ({n_cl} clusters)</h2>')
        parts.append(f'<div class="resolution-header">'
                     f'Leiden r={res} → <span class="stats-highlight">{n_cl} clusters</span></div>')

        # UMAP
        parts.append("<h3>UMAP — Leiden Clusters</h3>")
        parts.append(f'<div class="fig-container">{_fig_tag(figs["umap_leiden"])}</div>')

        # Dotplot
        if "dotplot" in figs:
            parts.append("<h3>Dotplot — Literature Markers</h3>")
            parts.append(f'<div class="fig-container">{_fig_tag(figs["dotplot"])}</div>')

        # Gene UMAPs
        if "gene_umaps" in figs:
            parts.append("<h3>Gene UMAPs — Key Markers</h3>")
            parts.append(f'<div class="fig-container">{_fig_tag(figs["gene_umaps"])}</div>')

        # DEG table
        parts.append("<h3>Top 10 DEGs per Cluster (Wilcoxon)</h3>")
        parts.append(figs["deg_table"])
        if "degs_csv" in figs:
            parts.append(f'<p style="font-size:11px;color:#666;">Full top 50: {figs["degs_csv"]}</p>')

        # Purity table
        if figs.get("purity_table"):
            parts.append("<h3>CellAssign Purity & QC Summary</h3>")
            parts.append(figs["purity_table"])

        # QC violins
        parts.append("<h3>QC Violins</h3>")
        parts.append(f'<div class="fig-container">{_fig_tag(figs["qc_violins"])}</div>')

        # Stacked bars
        if figs.get("stacked_bars"):
            parts.append("<h3>Metadata Composition</h3>")
            parts.append('<div class="grid-2col">')
            for meta_key, b64 in figs["stacked_bars"].items():
                parts.append(f'<div class="fig-container">{_fig_tag(b64, width="100%")}</div>')
            parts.append("</div>")

    # --- Low quality assessment ---
    parts.append('<div class="section-divider"></div>')
    n_res = len(resolutions)
    parts.append(f"<h2>{n_res + 3}. Low Quality Assessment</h2>")
    parts.append(f"""<div class="lq-box">
<p><b>MALAT1 + NEAT1 Low-Quality Score</b></p>
<p>Median lq_score: <b>{lq_median:.4f}</b></p>
<p>95th percentile: <b>{lq_95:.4f}</b></p>
<p>Cells above 95th pct: flagged as potential low-quality / dissociation artifacts</p>
<p><i>This cutoff should be applied consistently across all 13 compartments.</i></p>
</div>""")

    parts.append("</body></html>")
    return "\n".join(parts)


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print(f"{'='*70}")
    print(f"  07_resolution_explorer.py — HGSC Atlas")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}")

    # Parse command-line args
    if len(sys.argv) > 1:
        requested = [a.lower().strip() for a in sys.argv[1:]]
        to_run = []
        for name in requested:
            if name in COMPARTMENTS:
                to_run.append(name)
            else:
                print(f"ERROR: Unknown compartment '{name}'")
                print(f"Available: {', '.join(COMPARTMENTS.keys())}")
                sys.exit(1)
    else:
        to_run = RUN_ORDER

    print(f"\nCompartments to run ({len(to_run)}): {', '.join(to_run)}")
    print(f"Loading atlas: {ATLAS_H5AD}")

    atlas = ad.read_h5ad(ATLAS_H5AD)
    print(f"Atlas loaded: {atlas.shape[0]:,} cells x {atlas.shape[1]:,} genes\n")

    results = {}
    t_total = time.time()

    for comp_key in to_run:
        config = COMPARTMENTS[comp_key]
        try:
            html_path = process_compartment(atlas, comp_key, config)
            results[comp_key] = ("SUCCESS", html_path)
        except Exception as e:
            import traceback
            traceback.print_exc()
            results[comp_key] = ("FAILED", str(e))
        gc.collect()

    # Summary
    elapsed_total = time.time() - t_total
    print(f"\n{'='*70}")
    print(f"  SUMMARY — {elapsed_total/60:.1f} min total")
    print(f"{'='*70}")
    for comp, (status, detail) in results.items():
        print(f"  {comp:15s}: {status}  {detail}")
    print()

    n_ok = sum(1 for s, _ in results.values() if s == "SUCCESS")
    n_fail = sum(1 for s, _ in results.values() if s == "FAILED")
    print(f"  {n_ok}/{len(results)} succeeded, {n_fail} failed")
