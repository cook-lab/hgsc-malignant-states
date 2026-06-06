#!/usr/bin/env python3
"""
Epithelial NMF — SecA/SecB polarization programs
================================================
HGSC malignant-states atlas backend.

Non-negative Matrix Factorization on epithelial cells to learn the gene
programs underlying the SecA (progenitor-like) -> SecB (adaptive/differentiated)
polarization axis. Trains on a stratified 50k-cell subsample (top 3000 HVGs),
runs NMF (k=10) on log-normalized expression, then projects all epithelial
cells onto the learned factors. Factor_2 is the canonical SecB-defining factor
used downstream by 02_prepare_nmf_labels.py.

INPUTS:
  - <data_root>/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
    (epithelial subset; layers["counts"], obsm["X_umap_local"])
  - SecA/SecB signatures from shared/signatures.yml (used only to label
    the learned factors; the NMF is unsupervised)

OUTPUTS (output_root/03_epithelial_nmf/):
  - 11d_nmf_usage.csv        normalized per-cell factor usage (KEY cache, Supp Data 4)
  - 11d_nmf_usage_raw.csv    raw per-cell factor usage
  - 11d_nmf_loadings.csv     gene loadings (H matrix)
  - 11d_factor_scores.csv    SecA/SecB enrichment scores per factor
  - 11d_*.svg/pdf            factor UMAPs, polarization scatter, loadings, violins

MANUSCRIPT PANELS: feeds Fig 1G/H, Fig 1I, Fig 2A/C, SF5, SF6, Supp Data 4
  (the 11d_nmf_usage.csv Factor_2 column is the polarization substrate).

RUNTIME TIER: heavy (loads full epithelial counts; NMF fit + projection).

SEEDING: NMF inits, training subsample, and plot subsamples all use SEED
  from config (np.random.seed(SEED) + RandomState(SEED)) for determinism.

Usage:
    python 01_epithelial_nmf.py
    python 01_epithelial_nmf.py --k 8
    python 01_epithelial_nmf.py --figures-only
"""

import argparse
import gc
import os
import sys
import time
import warnings

import anndata as ad
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse
from sklearn.decomposition import NMF

# config is 3 levels up: atlas/03_epithelial_nmf/01_epithelial_nmf.py -> repo root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ============================================================================
# COOK LAB STYLE GUIDE v1.2
# ============================================================================

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       7,
    "axes.titlesize":  8,
    "axes.labelsize":  7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.dpi":      450,
    "savefig.dpi":     450,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
    "savefig.bbox":    "tight",
})

# ============================================================================
# PATHS
# ============================================================================

H5AD    = path("data_root", "2026_final_atlas", "celltype_h5ad",
               "hgsc_atlas_final_epithelial.h5ad")
OUT_DIR = path("output_root", "03_epithelial_nmf")

# ============================================================================
# GENE SIGNATURES (shared/signatures.yml — noBCAM 7-gene set)
# ============================================================================

_SIG_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "shared", "signatures.yml")
with open(os.path.abspath(_SIG_PATH)) as _fh:
    _SIGS = yaml.safe_load(_fh)
SECA_GENES = list(_SIGS["SecA"])
SECB_GENES = list(_SIGS["SecB"])

# ============================================================================
# EPITHELIAL PALETTE (celltype_level2; "Transitioning" -> "Intermediate")
# ============================================================================

EPI_PALETTE = {
    "Adaptive secretory epithelial cell":        "#B8741A",
    "Ciliated epithelial cell":                  "#E05A2C",
    "Cycling secretory epithelial cell":         "#F6D28B",
    "Secretory epithelial cell":                 "#E6A141",
    "Stress-response secretory epithelial cell": "#D9C5A2",
    "Intermediate epithelial cell":              "#7D4E4E",
}

# ============================================================================
# NMF PARAMETERS (validated — do not change)
# ============================================================================

N_TRAIN   = 50_000   # cells for NMF training
N_HVG     = 3000     # highly variable genes
DEFAULT_K = 10       # number of factors


# ============================================================================
# HELPERS
# ============================================================================

def save_fig(fig, name, out_dir=OUT_DIR):
    svg = os.path.join(out_dir, f"{name}.svg")
    pdf = os.path.join(out_dir, f"{name}.pdf")
    fig.savefig(svg, format="svg", bbox_inches="tight")
    fig.savefig(pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {name}.svg / .pdf")


def theme_void(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")


# ============================================================================
# LOAD & PREPROCESS
# ============================================================================

def load_epithelial():
    """Load epithelial h5ad via backed mode, extract needed data."""
    print("  [LOAD] Opening epithelial h5ad (backed)...", flush=True)
    t0 = time.time()

    bdata = ad.read_h5ad(H5AD, backed="r")
    n_cells = bdata.n_obs
    genes = bdata.var_names.tolist()

    level2 = bdata.obs["celltype_level2"].values.astype(str)
    barcodes = bdata.obs.index.values.copy()
    umap_local = np.array(bdata.obsm["X_umap_local"])

    print(f"    Extracting counts matrix ({n_cells:,} x {len(genes):,})...", flush=True)
    counts = bdata.layers["counts"][:, :]

    bdata.file.close()
    del bdata
    gc.collect()

    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"celltype_level2": level2}, index=barcodes),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["X_umap_local"] = umap_local

    print(f"    {n_cells:,} cells x {len(genes):,} genes loaded in {time.time()-t0:.0f}s")
    return adata


def preprocess_for_nmf(adata):
    """Normalize, log-transform, select HVGs. Returns log-normalized adata."""
    print("  [PREPROCESS] Normalize + HVG selection...", flush=True)
    t0 = time.time()

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3",
                                layer="counts")
    n_hvg = adata.var["highly_variable"].sum()
    print(f"    {n_hvg:,} HVGs selected")

    # Ensure SecA/SecB genes are included
    for gene in SECA_GENES + SECB_GENES:
        if gene in adata.var_names:
            adata.var.loc[gene, "highly_variable"] = True

    n_hvg_final = adata.var["highly_variable"].sum()
    if n_hvg_final > n_hvg:
        print(f"    Added {n_hvg_final - n_hvg} SecA/SecB genes -> {n_hvg_final:,} total")

    print(f"    Done in {time.time()-t0:.0f}s")
    return adata


# ============================================================================
# NMF
# ============================================================================

def run_nmf(adata, k=DEFAULT_K):
    """Downsample, train NMF, project all cells. Seeded from config."""
    print(f"  [NMF] k={k}, training on {N_TRAIN:,} cells...", flush=True)
    t0 = time.time()

    np.random.seed(SEED)
    rng = np.random.RandomState(SEED)
    level2 = adata.obs["celltype_level2"].values

    # Stratified downsample for training
    train_idx = []
    subtypes = np.unique(level2)
    for st in subtypes:
        st_idx = np.where(level2 == st)[0]
        n_st = len(st_idx)
        n_sample = max(100, int(N_TRAIN * n_st / len(level2)))
        n_sample = min(n_sample, n_st)
        train_idx.extend(rng.choice(st_idx, size=n_sample, replace=False).tolist())

    train_idx = np.array(sorted(train_idx))
    print(f"    Training set: {len(train_idx):,} cells")

    hvg_mask = adata.var["highly_variable"].values
    hvg_genes = adata.var_names[hvg_mask].tolist()

    X_train = adata[train_idx, :][:, hvg_mask].X
    if sparse.issparse(X_train):
        X_train = X_train.toarray()
    X_train = np.clip(X_train, 0, None)

    print(f"    Training matrix: {X_train.shape[0]:,} x {X_train.shape[1]:,}", flush=True)

    t_fit = time.time()
    model = NMF(
        n_components=k,
        init="nndsvda",
        random_state=SEED,
        max_iter=500,
    )
    model.fit_transform(X_train)
    H = model.components_  # (k, n_hvg)
    recon_err = model.reconstruction_err_
    print(f"    NMF fit in {time.time()-t_fit:.0f}s (recon error: {recon_err:.2f})")

    del X_train
    gc.collect()

    print("    Projecting all cells...", flush=True)
    t_proj = time.time()
    X_all = adata[:, hvg_mask].X
    if sparse.issparse(X_all):
        X_all = X_all.toarray()
    X_all = np.clip(X_all, 0, None)

    W_all = model.transform(X_all)
    print(f"    Projection done in {time.time()-t_proj:.0f}s")

    del X_all
    gc.collect()

    W_norm = W_all / (W_all.sum(axis=1, keepdims=True) + 1e-10)

    factor_names = [f"Factor_{i+1}" for i in range(k)]
    usage_df = pd.DataFrame(W_norm, index=adata.obs.index, columns=factor_names)
    usage_raw = pd.DataFrame(W_all, index=adata.obs.index, columns=factor_names)
    loadings_df = pd.DataFrame(H, index=factor_names, columns=hvg_genes)

    print(f"    NMF complete ({time.time()-t0:.0f}s)")
    return usage_df, usage_raw, loadings_df, model


# ============================================================================
# IDENTIFY SECA / SECB FACTORS
# ============================================================================

def identify_polarization_factors(loadings_df):
    """Score each factor for SecA and SecB gene enrichment."""
    print("  [IDENTIFY] Scoring factors for SecA/SecB enrichment...")

    scores = []
    for factor in loadings_df.index:
        loadings = loadings_df.loc[factor]
        ranked = loadings.sort_values(ascending=False)
        top100 = set(ranked.index[:100])
        top200 = set(ranked.index[:200])

        factor_max = loadings.max()
        secA_mean = loadings[loadings.index.isin(SECA_GENES)].mean() / factor_max
        secB_mean = loadings[loadings.index.isin(SECB_GENES)].mean() / factor_max

        scores.append({
            "factor": factor,
            "secA_top100": len(set(SECA_GENES) & top100),
            "secA_top200": len(set(SECA_GENES) & top200),
            "secB_top100": len(set(SECB_GENES) & top100),
            "secB_top200": len(set(SECB_GENES) & top200),
            "secA_loading": secA_mean,
            "secB_loading": secB_mean,
        })

    scores_df = pd.DataFrame(scores).set_index("factor")
    best_secA = scores_df["secA_loading"].idxmax()
    best_secB = scores_df["secB_loading"].idxmax()

    print(f"    Best SecA factor: {best_secA} "
          f"(loading={scores_df.loc[best_secA, 'secA_loading']:.3f}, "
          f"{scores_df.loc[best_secA, 'secA_top100']} genes in top100)")
    print(f"    Best SecB factor: {best_secB} "
          f"(loading={scores_df.loc[best_secB, 'secB_loading']:.3f}, "
          f"{scores_df.loc[best_secB, 'secB_top100']} genes in top100)")

    return scores_df, best_secA, best_secB


# ============================================================================
# FIGURES
# ============================================================================

def plot_factor_umaps(usage_df, umap_coords, k, out_dir):
    """UMAP colored by each NMF factor usage."""
    print("  [FIG] Factor usage UMAPs...", flush=True)

    ncols = 4
    nrows = (k + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.8 * nrows))
    axes = axes.flatten()

    rng = np.random.RandomState(SEED)
    n = len(usage_df)
    n_plot = min(150_000, n)
    idx = rng.choice(n, size=n_plot, replace=False)

    for i in range(k):
        ax = axes[i]
        factor = f"Factor_{i+1}"
        vals = usage_df[factor].values[idx]
        sc_plot = ax.scatter(
            umap_coords[idx, 0], umap_coords[idx, 1],
            c=vals, cmap="YlOrRd", s=0.3, alpha=0.6,
            edgecolors="none", rasterized=True,
            vmin=np.percentile(vals, 1), vmax=np.percentile(vals, 99),
        )
        theme_void(ax)
        ax.set_title(factor, fontsize=8, fontweight="bold")
        plt.colorbar(sc_plot, ax=ax, shrink=0.6, aspect=20, pad=0.02)

    for j in range(k, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("NMF Factor Usage on Epithelial UMAP", fontsize=10, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "11d_A_nmf_factor_umaps", out_dir)


def plot_polarization(usage_df, level2, best_secA, best_secB, out_dir):
    """SecA vs SecB factor usage scatter — the key polarization plot."""
    print("  [FIG] SecA/SecB polarization scatter...", flush=True)

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    rng = np.random.RandomState(SEED)
    subtypes = level2.values

    counts = pd.Series(subtypes).value_counts()
    for st in counts.index:
        mask = subtypes == st
        idx = np.where(mask)[0]
        n_plot = min(len(idx), 20_000)
        idx_sub = rng.choice(idx, size=n_plot, replace=False)
        ax.scatter(
            usage_df[best_secA].values[idx_sub],
            usage_df[best_secB].values[idx_sub],
            c=[EPI_PALETTE.get(st, "#999999")], s=0.8, alpha=0.5,
            edgecolors="none", rasterized=True, label=st,
        )

    ax.set_xlabel(f"{best_secA} usage (SecA / progenitor-like)", fontsize=8)
    ax.set_ylabel(f"{best_secB} usage (SecB / adaptive)", fontsize=8)
    ax.set_title("Epithelial SecA <-> SecB Polarization", fontsize=9, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [Line2D([0], [0], marker="o", color="none",
                      markerfacecolor=EPI_PALETTE[st], markeredgecolor="none",
                      markersize=5, label=st)
               for st in sorted(EPI_PALETTE.keys())]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, fontsize=5.5, handletextpad=0.3, labelspacing=0.5)

    plt.tight_layout()
    save_fig(fig, "11d_B_secA_secB_polarization", out_dir)


def plot_gene_loadings(loadings_df, best_secA, best_secB, out_dir):
    """Top gene loadings for SecA and SecB factors."""
    print("  [FIG] Gene loading barplots...", flush=True)

    fig, axes = plt.subplots(1, 2, figsize=(9, 5))
    for ax, factor, sig_genes, sig_name in [
        (axes[0], best_secA, SECA_GENES, "SecA (progenitor)"),
        (axes[1], best_secB, SECB_GENES, "SecB (adaptive)"),
    ]:
        loadings = loadings_df.loc[factor].sort_values(ascending=False)
        top30 = loadings.head(30)
        colors = ["#D14E6C" if g in sig_genes else "#888888" for g in top30.index]
        ax.barh(range(len(top30)), top30.values[::-1], color=colors[::-1],
                edgecolor="none", height=0.7)
        ax.set_yticks(range(len(top30)))
        ax.set_yticklabels(top30.index[::-1], fontsize=5.5)
        ax.set_xlabel("Loading weight", fontsize=7)
        ax.set_title(f"{factor} — {sig_name}\nTop 30 genes", fontsize=8, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(0.95, 0.05, f"Red = {sig_name} signature",
                transform=ax.transAxes, fontsize=5.5, ha="right",
                color="#D14E6C", fontweight="bold")

    plt.tight_layout()
    save_fig(fig, "11d_C_gene_loadings", out_dir)


def plot_factor_by_subtype(usage_df, level2, best_secA, best_secB, out_dir):
    """Violin plots of SecA/SecB factor usage by subtype."""
    print("  [FIG] Factor usage by subtype (violin)...", flush=True)
    import seaborn as sns

    rng = np.random.RandomState(SEED)
    n = len(usage_df)
    n_plot = min(100_000, n)
    idx = rng.choice(n, size=n_plot, replace=False)

    plot_df = pd.DataFrame({
        "SecA_usage": usage_df[best_secA].values[idx],
        "SecB_usage": usage_df[best_secB].values[idx],
        "Subtype": level2.values[idx],
    })
    order = (plot_df.groupby("Subtype")["SecA_usage"].mean()
             .sort_values(ascending=False).index.tolist())

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    palette = {st: EPI_PALETTE.get(st, "#999999") for st in order}

    for ax, col, title in [
        (axes[0], "SecA_usage", f"{best_secA} (SecA / progenitor)"),
        (axes[1], "SecB_usage", f"{best_secB} (SecB / adaptive)"),
    ]:
        sns.violinplot(data=plot_df, y="Subtype", x=col, order=order,
                       palette=palette, ax=ax, cut=0,
                       linewidth=0.5, inner="quartile", density_norm="width")
        ax.set_title(title, fontsize=8, fontweight="bold")
        ax.set_xlabel("Factor usage (normalized)", fontsize=7)
        ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=5.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    plt.tight_layout()
    save_fig(fig, "11d_D_factor_by_subtype", out_dir)


def plot_secA_secB_umap(usage_df, umap_coords, best_secA, best_secB, out_dir):
    """Side-by-side UMAP: SecA usage vs SecB usage."""
    print("  [FIG] SecA/SecB UMAP comparison...", flush=True)

    rng = np.random.RandomState(SEED)
    n = len(usage_df)
    n_plot = min(150_000, n)
    idx = rng.choice(n, size=n_plot, replace=False)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, factor, title, cmap in [
        (axes[0], best_secA, "SecA (progenitor-like)", "YlOrBr"),
        (axes[1], best_secB, "SecB (adaptive)", "OrRd"),
    ]:
        vals = usage_df[factor].values[idx]
        sc_plot = ax.scatter(
            umap_coords[idx, 0], umap_coords[idx, 1],
            c=vals, cmap=cmap, s=0.4, alpha=0.6,
            edgecolors="none", rasterized=True,
            vmin=np.percentile(vals, 2), vmax=np.percentile(vals, 98),
        )
        theme_void(ax)
        ax.set_title(title, fontsize=8, fontweight="bold")
        plt.colorbar(sc_plot, ax=ax, shrink=0.6, aspect=20, pad=0.02)

    fig.suptitle("NMF Factor Usage — SecA/SecB Polarization", fontsize=10, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig(fig, "11d_E_secA_secB_umap", out_dir)


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Epithelial NMF — SecA/SecB polarization")
    parser.add_argument("--k", type=int, default=DEFAULT_K,
                        help=f"Number of NMF factors (default: {DEFAULT_K})")
    parser.add_argument("--figures-only", action="store_true",
                        help="Load saved NMF results and regenerate figures")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    t0 = time.time()
    print("=" * 65)
    print("  Epithelial NMF: SecA/SecB Polarization")
    print("=" * 65)

    usage_path = os.path.join(OUT_DIR, "11d_nmf_usage.csv")
    loadings_path = os.path.join(OUT_DIR, "11d_nmf_loadings.csv")

    if args.figures_only and os.path.exists(usage_path) and os.path.exists(loadings_path):
        print("\n  [LOAD] Reading saved NMF results...")
        usage_df = pd.read_csv(usage_path, index_col=0)
        loadings_df = pd.read_csv(loadings_path, index_col=0)
        k = len(loadings_df)

        print("  [LOAD] Extracting UMAP + metadata from h5ad...")
        with h5py.File(H5AD, "r") as f:
            umap_coords = np.array(f["obsm"]["X_umap_local"])
            level2 = pd.Series(
                [x.decode() if isinstance(x, bytes) else x
                 for x in f["obs"]["celltype_level2"][:]],
                index=[x.decode() if isinstance(x, bytes) else x
                       for x in f["obs"]["_index"][:]]
            )
    else:
        adata = load_epithelial()
        adata = preprocess_for_nmf(adata)
        usage_df, usage_raw, loadings_df, model = run_nmf(adata, k=args.k)
        k = args.k

        umap_coords = adata.obsm["X_umap_local"]
        level2 = adata.obs["celltype_level2"]

        print("\n  [SAVE] Writing NMF results...")
        usage_df.to_csv(usage_path)
        usage_raw.to_csv(os.path.join(OUT_DIR, "11d_nmf_usage_raw.csv"))
        loadings_df.to_csv(loadings_path)
        print(f"    Saved to {OUT_DIR}")

        del adata
        gc.collect()

    scores_df, best_secA, best_secB = identify_polarization_factors(loadings_df)
    scores_df.to_csv(os.path.join(OUT_DIR, "11d_factor_scores.csv"))

    print(f"\n  [FIGURES] Generating to {OUT_DIR}")
    plot_factor_umaps(usage_df, umap_coords, k, OUT_DIR)
    plot_secA_secB_umap(usage_df, umap_coords, best_secA, best_secB, OUT_DIR)
    plot_polarization(usage_df, level2, best_secA, best_secB, OUT_DIR)
    plot_gene_loadings(loadings_df, best_secA, best_secB, OUT_DIR)
    plot_factor_by_subtype(usage_df, level2, best_secA, best_secB, OUT_DIR)

    print(f"\n{'=' * 65}")
    print(f"  COMPLETE ({time.time() - t0:.0f}s)")
    print(f"  Output: {OUT_DIR}")
    print(f"  NMF: k={k}, best SecA={best_secA}, best SecB={best_secB}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
