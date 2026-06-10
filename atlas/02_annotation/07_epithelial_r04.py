#!/usr/bin/env python3
"""
Atlas 02 — Step 07 (epithelial r0.4): finalize epithelial Leiden at resolution 0.4

PURPOSE
    Subset epithelial cells from the level-1 atlas, recompute neighbours/UMAP on the
    scANVI latent space, run Leiden at the chosen resolution (0.4), and score
    clusters (silhouette, SecA/SecB polarization markers, top markers). Emits the
    epithelial cluster stats CSVs + an HTML report used to assign level-2 epithelial
    labels.

INPUTS
    DATA_ROOT/2026_final_atlas/hgsc_atlas_celltype_level1.h5ad
    shared/signatures.yml  (canonical 7-gene noBCAM SecA/SecB sets)

OUTPUTS
    output_root/02_annotation/07_resolution_explorer/{csvs,figs}/*, *_report.html

MANUSCRIPT PANEL(S)
    Annotation backend (epithelial resolution); feeds level-2 epithelial labels.

RUNTIME TIER
    heavy (neighbours + UMAP + Leiden on epithelial subset).
"""

import base64
import gc
import io
import os
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
# RESOLUTION TO TEST
# ============================================================================

RESOLUTION = 0.4
LEIDEN_KEY = f"leiden_r{RESOLUTION}"

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

# ============================================================================
# LITERATURE-DERIVED GENE SETS (canonical epithelial markers)
# ============================================================================

GENE_SETS = OrderedDict([
    ("SecA (Progenitor)",     SECA_SIGNATURE),
    ("SecB (Differentiated)", SECB_SIGNATURE),
    ("Ciliated",              ["FOXJ1", "CAPS", "TPPP3", "RSPH1", "ZMYND10", "DNAH5"]),
    ("Proliferative",         ["MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1", "CENPF"]),
    ("EMT",                   ["VIM", "CDH2", "SNAI1", "ZEB1", "FN1", "TWIST1"]),
    ("Hypoxia",               ["VEGFA", "BNIP3", "SLC2A1", "LDHA", "CA9", "HIF1A"]),
])

PURITY_TYPES = ["Epithelial"]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _adaptive_pt(n_cells):
    return max(0.8, min(12.0, 30000 / max(n_cells, 1)))


def _svg_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", dpi=DPI, bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return encoded


def _fig_tag(b64, width="100%"):
    return f'<img src="data:image/svg+xml;base64,{b64}" style="width:{width}; max-width:1200px;" />'


def _save_svg(fig, name):
    path = os.path.join(FIG_DIR, f"{name}.svg")
    fig.savefig(path, format="svg", dpi=DPI, bbox_inches="tight")
    b64 = _svg_to_b64(fig)
    return b64


# ---------- UMAP ----------

def plot_umap_cat(adata, color_key, title, pt_size=None, palette=None):
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
    if n_cats <= 20:
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=5,
                  markerscale=3, frameon=False, handletextpad=0.3)
    else:
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=4,
                  markerscale=2, frameon=False, ncol=2, handletextpad=0.2)
    fig.tight_layout()
    return fig


def plot_umap_cont(adata, gene, title=None, pt_size=None):
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

def compute_stats(adata, max_cells_silhouette=50000):
    """Compute clustering quality metrics for r=0.4."""
    latent = adata.obsm["X_scanvi"]
    labels = adata.obs[LEIDEN_KEY].values.astype(str)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    sizes = pd.Series(labels).value_counts().values
    row = {
        "resolution": RESOLUTION,
        "n_clusters": n_clusters,
        "min_cluster_size": int(sizes.min()),
        "max_cluster_size": int(sizes.max()),
        "size_entropy": round(float(scipy_entropy(sizes / sizes.sum())), 4),
    }

    if n_clusters > 1:
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

        if "study" in adata.obs.columns:
            row["study_ari"] = round(float(
                adjusted_rand_score(labels, adata.obs["study"].values.astype(str))
            ), 4)
    else:
        row["silhouette"] = None
        row["calinski_harabasz"] = None
        row["davies_bouldin"] = None
        row["study_ari"] = None

    return pd.DataFrame([row])


# ---------- DEG TABLE HTML ----------

def _deg_table_html(markers_df, cluster_key, adata):
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
# MAIN PIPELINE
# ============================================================================

def main():
    t0 = time.time()
    prefix = "epithelial_r0.4"

    print(f"{'='*70}")
    print(f"  07_epithelial_r04.py — Epithelial r=0.4 Resolution Explorer")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}\n")

    # ------------------------------------------------------------------
    # 1. Load and subset
    # ------------------------------------------------------------------
    print("1. Loading atlas...", flush=True)
    atlas = ad.read_h5ad(ATLAS_H5AD)
    print(f"   Atlas: {atlas.shape[0]:,} x {atlas.shape[1]:,}")

    print("   Subsetting epithelial...", flush=True)
    mask = atlas.obs["celltype_level1"] == "Epithelial"
    adata = atlas[mask].copy()
    del atlas
    gc.collect()
    print(f"   Epithelial: {len(adata):,} cells")

    # ------------------------------------------------------------------
    # 2. Separate ciliated cells
    # ------------------------------------------------------------------
    print("\n2. Splitting ciliated cells...", flush=True)
    ciliated_adata = None
    n_cil = 0
    if "celltype_pred" in adata.obs.columns:
        cil_mask = adata.obs["celltype_pred"] == "Ciliated"
        n_cil = cil_mask.sum()
        if n_cil > 0:
            ciliated_adata = adata[cil_mask].copy()
            adata = adata[~cil_mask].copy()
            print(f"   Ciliated: {n_cil:,} cells (separated)")
            print(f"   Secretory/other: {len(adata):,} cells (for clustering)")
        else:
            print("   No ciliated cells found by celltype_pred")
    else:
        print("   No celltype_pred column — using all cells")

    n_cells = len(adata)

    # ------------------------------------------------------------------
    # 3. Neighbors + UMAP from scANVI
    # ------------------------------------------------------------------
    print("\n3. Computing neighbors → UMAP from X_scanvi...", flush=True)
    sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=15)
    print("   Neighbors done", flush=True)
    sc.tl.umap(adata, random_state=SEED)
    print("   UMAP done", flush=True)

    # ------------------------------------------------------------------
    # 4. Leiden r=0.4
    # ------------------------------------------------------------------
    print(f"\n4. Leiden clustering at r={RESOLUTION}...", flush=True)
    sc.tl.leiden(adata, resolution=RESOLUTION, key_added=LEIDEN_KEY, random_state=SEED)
    n_cl = adata.obs[LEIDEN_KEY].nunique()
    print(f"   r={RESOLUTION}: {n_cl} clusters")

    # ------------------------------------------------------------------
    # 5. MALAT1/NEAT1 low-quality score
    # ------------------------------------------------------------------
    print("\n5. Computing lq_score (MALAT1 + NEAT1)...", flush=True)
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
    # 6. Resolution statistics
    # ------------------------------------------------------------------
    print("\n6. Computing resolution statistics...", flush=True)
    stats_df = compute_stats(adata)
    stats_csv = os.path.join(CSV_DIR, f"{prefix}_stats.csv")
    stats_df.to_csv(stats_csv, index=False)
    print(f"   Saved: {stats_csv}")
    print(stats_df.to_string(index=False))

    # ------------------------------------------------------------------
    # 7. Generate figures
    # ------------------------------------------------------------------
    print(f"\n7. Generating figures for r={RESOLUTION} ({n_cl} clusters)...", flush=True)
    pt = _adaptive_pt(n_cells)
    figs = {}

    # 7a. UMAP colored by leiden
    print("   UMAP leiden...", flush=True)
    fig = plot_umap_cat(adata, LEIDEN_KEY, f"Epithelial — Leiden r={RESOLUTION}", pt_size=pt)
    figs["umap_leiden"] = _save_svg(fig, f"{prefix}_umap_leiden")

    # 7b. UMAP colored by study
    if "study" in adata.obs.columns:
        print("   UMAP study...", flush=True)
        fig = plot_umap_cat(adata, "study", f"Epithelial — Study (r={RESOLUTION})", pt_size=pt)
        figs["umap_study"] = _save_svg(fig, f"{prefix}_umap_study")

    # 7c. Marker dotplot
    print("   Dotplot...", flush=True)
    fig = plot_marker_dotplot(adata, GENE_SETS, LEIDEN_KEY,
                              title=f"Epithelial — Literature Markers (r={RESOLUTION})")
    figs["dotplot"] = _save_svg(fig, f"{prefix}_dotplot")

    # 7d. Gene UMAPs (first gene from each group, 2x4 grid)
    print("   Gene UMAPs...", flush=True)
    key_genes = []
    for group_name, genes in GENE_SETS.items():
        for g in genes:
            if g in adata.var_names:
                key_genes.append((group_name, g))
                break
    genes_to_plot = key_genes[:8]
    gene_umaps_b64 = None
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
        for idx in range(len(genes_to_plot), nrows * ncols):
            r, c = divmod(idx, ncols)
            axes[r, c].set_visible(False)
        fig_genes.suptitle(f"Epithelial — Key Markers (r={RESOLUTION})", fontsize=9, y=1.01)
        fig_genes.tight_layout()
        gene_umaps_b64 = _save_svg(fig_genes, f"{prefix}_gene_umaps")

    # 7e. Wilcoxon DEGs (top 50)
    print("   Wilcoxon DEGs...", flush=True)
    markers_df = pd.DataFrame()
    deg_table_html = "<p>DEG computation failed.</p>"
    try:
        sc.tl.rank_genes_groups(adata, groupby=LEIDEN_KEY, method="wilcoxon",
                                 n_genes=50, use_raw=False)
        marker_rows = []
        for cl in sorted(adata.obs[LEIDEN_KEY].unique(), key=lambda x: int(x) if x.isdigit() else x):
            try:
                names  = adata.uns["rank_genes_groups"]["names"][cl][:50]
                scores = adata.uns["rank_genes_groups"]["scores"][cl][:50]
                logfc  = adata.uns["rank_genes_groups"]["logfoldchanges"][cl][:50]
                pvals  = adata.uns["rank_genes_groups"]["pvals_adj"][cl][:50]
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
        csv_path = os.path.join(CSV_DIR, f"{prefix}_markers_top50.csv")
        markers_df.to_csv(csv_path, index=False)
        print(f"   Saved: {csv_path}")
        deg_table_html = _deg_table_html(markers_df, LEIDEN_KEY, adata)
    except Exception as e:
        print(f"   WARNING: DEG computation failed: {e}")

    # 7f. QC violins
    print("   QC violins...", flush=True)
    qc_metrics = []
    if "n_genes_by_counts" in adata.obs.columns:
        qc_metrics.append(("n_genes_by_counts", "Genes Detected"))
    if "total_counts" in adata.obs.columns:
        qc_metrics.append(("total_counts", "UMI Counts"))
    if "doublet_score_scrublet" in adata.obs.columns:
        qc_metrics.append(("doublet_score_scrublet", "Doublet Score"))
    qc_metrics.append(("lq_score", "LQ Score (MALAT1+NEAT1)"))

    fig = plot_qc_violins(adata, LEIDEN_KEY, qc_metrics, title_prefix=f"r={RESOLUTION} — ")
    figs["qc_violins"] = _save_svg(fig, f"{prefix}_qc_violins")

    # 7g. CellAssign purity table
    purity_html = ""
    if "celltype_pred" in adata.obs.columns:
        print("   Purity table...", flush=True)
        purity_rows = []
        for cl in sorted(adata.obs[LEIDEN_KEY].unique(), key=lambda x: int(x) if x.isdigit() else x):
            cl_mask = adata.obs[LEIDEN_KEY] == cl
            n_in_cl = cl_mask.sum()
            pred = adata.obs.loc[cl_mask, "celltype_pred"]
            purity_pct = sum((pred == pt_name).sum() for pt_name in PURITY_TYPES)
            purity_pct = 100.0 * purity_pct / max(n_in_cl, 1)
            top_pred = pred.value_counts().head(3)
            top_str = ", ".join([f"{v}:{c}" for v, c in top_pred.items()])
            purity_rows.append({
                "Cluster": cl, "N": f"{n_in_cl:,}",
                "Purity %": f"{purity_pct:.1f}%",
                "Top predictions": top_str,
                "Med genes": f"{adata.obs.loc[cl_mask, 'n_genes_by_counts'].median():.0f}",
                "Med lq": f"{adata.obs.loc[cl_mask, 'lq_score'].median():.4f}",
            })
        purity_df = pd.DataFrame(purity_rows)
        purity_html = purity_df.to_html(index=False, classes="styled-table", border=0)

    # 7h. Stacked bars
    print("   Stacked bars...", flush=True)
    bar_figs = {}
    for meta_key, meta_label in [("treatment_status", "Treatment Status"),
                                   ("anatomic_site", "Anatomic Site"),
                                   ("study", "Study")]:
        if meta_key in adata.obs.columns:
            fig = plot_stacked_bar(adata, LEIDEN_KEY, meta_key,
                                    f"Epithelial r={RESOLUTION} — {meta_label}")
            bar_figs[meta_key] = _save_svg(fig, f"{prefix}_bar_{meta_key}")

    # ------------------------------------------------------------------
    # 8. Build HTML report
    # ------------------------------------------------------------------
    print("\n8. Building HTML report...", flush=True)

    import re
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    ciliated_note = ""
    if ciliated_adata is not None:
        ciliated_note = f"<p>Ciliated cells separated: <b>{n_cil:,}</b> (not included in clustering)</p>"

    css = """
    <style>
        body { font-family: Arial, sans-serif; max-width: 1400px; margin: 0 auto; padding: 20px; background: #fafafa; }
        h1 { color: #333; border-bottom: 3px solid #0067A5; padding-bottom: 10px; }
        h2 { color: #0067A5; margin-top: 40px; border-bottom: 1px solid #ddd; padding-bottom: 5px; }
        h3 { color: #604E97; margin-top: 25px; }
        .overview-box { background: #f0f4f8; padding: 15px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #0067A5; }
        .stats-box { background: #e8f5e9; padding: 15px; border-radius: 8px; margin: 15px 0; border-left: 4px solid #008856; }
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
        .context-table { margin: 15px 0; }
        @media (max-width: 900px) { .grid-2col { grid-template-columns: 1fr; } }
    </style>
    """

    # Context: existing resolutions for comparison
    context_html = """
    <h3>Existing Resolution Results (for context)</h3>
    <table class="styled-table" border="0">
    <tr><th>Resolution</th><th>N Clusters</th><th>Min Size</th><th>Max Size</th>
        <th>Silhouette</th><th>Calinski-Harabasz</th><th>Davies-Bouldin</th><th>Study ARI</th></tr>
    <tr><td>0.1</td><td>4</td><td>410</td><td>610,836</td><td>-0.1017</td><td>37,244</td><td>1.5581</td><td>0.0261</td></tr>
    <tr><td>0.2</td><td>6</td><td>715</td><td>276,566</td><td>0.0025</td><td>54,515</td><td>1.9686</td><td>0.0183</td></tr>
    <tr><td>0.3</td><td>8</td><td>28</td><td>270,599</td><td>-0.0528</td><td>39,501</td><td>1.7999</td><td>0.0189</td></tr>
    <tr style="background:#fff8e1;font-weight:bold"><td>0.4</td><td colspan="7">← This report</td></tr>
    <tr><td>0.5</td><td>11</td><td>84</td><td>171,501</td><td>-0.0156</td><td>43,241</td><td>1.7385</td><td>0.0270</td></tr>
    <tr><td>0.8</td><td>18</td><td>20</td><td>89,320</td><td>0.0625</td><td>41,021</td><td>1.8942</td><td>0.0247</td></tr>
    </table>
    """

    # Stats row for r=0.4
    r04 = stats_df.iloc[0]
    stats_block = f"""
    <div class="stats-box">
    <h3>r=0.4 Statistics</h3>
    <table class="styled-table" border="0">
    <tr><th>Metric</th><th>Value</th></tr>
    <tr><td>N clusters</td><td><b>{int(r04['n_clusters'])}</b></td></tr>
    <tr><td>Min cluster size</td><td>{int(r04['min_cluster_size']):,}</td></tr>
    <tr><td>Max cluster size</td><td>{int(r04['max_cluster_size']):,}</td></tr>
    <tr><td>Size entropy</td><td>{r04['size_entropy']:.4f}</td></tr>
    <tr><td>Silhouette score</td><td>{r04['silhouette']:.4f}</td></tr>
    <tr><td>Calinski-Harabasz</td><td>{r04['calinski_harabasz']:.2f}</td></tr>
    <tr><td>Davies-Bouldin</td><td>{r04['davies_bouldin']:.4f}</td></tr>
    <tr><td>Study ARI (batch mixing)</td><td>{r04.get('study_ari', 'N/A')}</td></tr>
    </table>
    </div>
    """

    # Small clusters warning
    min_warn = ""
    if r04["min_cluster_size"] < 200:
        min_warn = f'<p class="warning">Smallest cluster has only {int(r04["min_cluster_size"])} cells</p>'

    parts = [f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Epithelial — Resolution 0.4 Explorer</title>
{css}
</head><body>
<h1>Epithelial — Resolution 0.4 Explorer</h1>

<div class="overview-box">
<p><b>Compartment:</b> Epithelial | <b>Cells for clustering:</b> {n_cells:,} |
<b>Resolution:</b> {RESOLUTION}</p>
{ciliated_note}
<p><b>Source:</b> hgsc_atlas_celltype_level1.h5ad (X_scanvi → neighbors → UMAP)</p>
<p><b>Generated:</b> {timestamp}</p>
</div>

<h2>1. Resolution Statistics</h2>
{stats_block}
{min_warn}
{context_html}

<div class="section-divider"></div>
<h2>2. UMAP — Leiden Clusters</h2>
<div class="resolution-header">
Leiden r={RESOLUTION} → <span class="stats-highlight">{n_cl} clusters</span>
</div>
<div class="fig-container">{_fig_tag(figs["umap_leiden"])}</div>
"""]

    # Study UMAP
    if "umap_study" in figs:
        parts.append('<h3>UMAP — Study</h3>')
        parts.append(f'<div class="fig-container">{_fig_tag(figs["umap_study"])}</div>')

    # Dotplot
    parts.append('<div class="section-divider"></div>')
    parts.append("<h2>3. Dotplot — Literature Markers</h2>")
    parts.append(f'<div class="fig-container">{_fig_tag(figs["dotplot"])}</div>')

    # Gene UMAPs
    if gene_umaps_b64:
        parts.append('<div class="section-divider"></div>')
        parts.append("<h2>4. Gene UMAPs — Key Markers</h2>")
        parts.append(f'<div class="fig-container">{_fig_tag(gene_umaps_b64)}</div>')

    # DEGs
    parts.append('<div class="section-divider"></div>')
    parts.append("<h2>5. Top 10 DEGs per Cluster (Wilcoxon)</h2>")
    parts.append(deg_table_html)
    if not markers_df.empty:
        parts.append(f'<p style="font-size:11px;color:#666;">Full top 50: csvs/{prefix}_markers_top50.csv</p>')

    # Purity
    if purity_html:
        parts.append('<div class="section-divider"></div>')
        parts.append("<h2>6. CellAssign Purity & QC Summary</h2>")
        parts.append(purity_html)

    # QC violins
    parts.append('<div class="section-divider"></div>')
    parts.append("<h2>7. QC Violins</h2>")
    parts.append(f'<div class="fig-container">{_fig_tag(figs["qc_violins"])}</div>')

    # LQ assessment
    parts.append(f"""<div class="lq-box">
<p><b>MALAT1 + NEAT1 Low-Quality Score</b></p>
<p>Median lq_score: <b>{lq_median:.4f}</b></p>
<p>95th percentile: <b>{lq_95:.4f}</b></p>
</div>""")

    # Stacked bars
    if bar_figs:
        parts.append('<div class="section-divider"></div>')
        parts.append("<h2>8. Metadata Composition</h2>")
        parts.append('<div class="grid-2col">')
        for meta_key, b64 in bar_figs.items():
            parts.append(f'<div class="fig-container">{_fig_tag(b64, width="100%")}</div>')
        parts.append("</div>")

    parts.append("</body></html>")

    html_path = os.path.join(OUT_DIR, f"{prefix}_report.html")
    with open(html_path, "w") as f:
        f.write("\n".join(parts))
    print(f"   Saved: {html_path}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  COMPLETE — {elapsed/60:.1f} min")
    print(f"{'='*70}")
    print(f"  Resolution: {RESOLUTION}")
    print(f"  Clusters:   {n_cl}")
    print(f"  Cells:      {n_cells:,} (+ {n_cil:,} ciliated separated)")
    print(f"  Report:     {html_path}")
    print(f"  Stats CSV:  {stats_csv}")
    print()

    del adata
    if ciliated_adata is not None:
        del ciliated_adata
    gc.collect()


if __name__ == "__main__":
    main()
