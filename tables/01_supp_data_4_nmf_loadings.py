#!/usr/bin/env python3
"""
Supp Data 4 — Epithelial NMF gene loadings (k=10)
=================================================
HGSC malignant-states atlas · supplemental table generator.

PURPOSE
    Run (or reload) the k=10 NMF decomposition of epithelial cells that defines
    the SecA/SecB polarization programs (canonical analysis step 11d), then emit
    the wide gene-loading matrix AND a tidy long-form table for Supp Data 4.
    Analytical logic is preserved verbatim from the canonical 11d script — this
    migration only re-routes paths/seeds/signatures and adds the tidy reshape.

INPUTS
    - obj("atlas_celltype_dir")/hgsc_atlas_final_epithelial.h5ad   (epithelial subset; counts layer + X_umap_local)
    - shared/signatures.yml  (SecA/SecB noBCAM 7-gene sets, used only to flag/retain signature genes)

OUTPUTS  (under output_root/11d_epithelial_nmf/ and output_root/supplemental/)
    - 11d_nmf_usage.csv          normalized per-cell factor usage (all cells)
    - 11d_nmf_usage_raw.csv      raw per-cell factor usage
    - 11d_nmf_loadings.csv       wide H matrix (factor x gene)            <- source matrix for SD4
    - 11d_factor_scores.csv      SecA/SecB enrichment scores per factor
    - supplemental/supplementary_table_NMF_factor_genes.csv   tidy SD4 (factor, gene, loading, rank, is_SecA_gene, is_SecB_gene)
    Diagnostic figures (11d_A..E svg/pdf) are also written, unchanged.

MANUSCRIPT PANEL(S)
    Supp Data 4 (NMF factor gene loadings). The 11d_nmf_usage.csv cache also feeds
    Fig 1G/1H (SecA Factor_3 / SecB Factor_2 UMAP colormaps) and SF5/SF6.

RUNTIME TIER
    heavy (full NMF over ~575k cells; use --figures-only to reshape from cached CSVs).
"""

import argparse
import gc
import os
import sys
import time
import warnings
from pathlib import Path

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

# --- central config (tables/ is 1 level below repo root) ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
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
# PATHS / SIGNATURES
# ============================================================================
H5AD    = os.path.join(obj("atlas_celltype_dir"), "hgsc_atlas_final_epithelial.h5ad")
OUT_DIR = path("output_root", "11d_epithelial_nmf")

_SIG = yaml.safe_load(open(Path(__file__).resolve().parents[1] / "shared" / "signatures.yml"))
SECA_GENES = list(_SIG["SecA"])   # noBCAM 7-gene set (shared source of truth)
SECB_GENES = list(_SIG["SecB"])   # noBCAM 7-gene set (shared source of truth)

# Epithelial palette ("Intermediate" replaces legacy "Transitioning").
EPI_PALETTE = {
    "Adaptive secretory epithelial cell":        "#B8741A",
    "Ciliated epithelial cell":                  "#E05A2C",
    "Cycling secretory epithelial cell":         "#F6D28B",
    "Secretory epithelial cell":                 "#E6A141",
    "Stress-response secretory epithelial cell": "#D9C5A2",
    "Intermediate epithelial cell":              "#7D4E4E",
}

# ============================================================================
# NMF PARAMETERS  (identical to canonical 11d)
# ============================================================================
N_TRAIN   = 50_000
N_HVG     = 3000
DEFAULT_K = 10
RNG_SEED  = SEED


def save_fig(fig, name, out_dir=OUT_DIR):
    fig.savefig(os.path.join(out_dir, f"{name}.svg"), format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, f"{name}.pdf"), format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"    Saved: {name}.svg / .pdf")


def theme_void(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal")


def _relabel_intermediate(level2):
    """Standardize legacy 'Transitioning' epithelial label to 'Intermediate'."""
    return pd.Series(level2).astype(str).str.replace(
        "Transitioning", "Intermediate", regex=False).values


# ============================================================================
# LOAD & PREPROCESS
# ============================================================================
def load_epithelial():
    print("  [LOAD] Opening epithelial h5ad (backed)...", flush=True)
    t0 = time.time()
    bdata = ad.read_h5ad(H5AD, backed="r")
    n_cells = bdata.n_obs
    genes = bdata.var_names.tolist()
    level2 = _relabel_intermediate(bdata.obs["celltype_level2"].values.astype(str))
    barcodes = bdata.obs.index.values.copy()
    umap_local = np.array(bdata.obsm["X_umap_local"])
    print(f"    Extracting counts matrix ({n_cells:,} x {len(genes):,})...", flush=True)
    counts = bdata.layers["counts"][:, :]
    bdata.file.close(); del bdata; gc.collect()
    adata = ad.AnnData(
        X=counts,
        obs=pd.DataFrame({"celltype_level2": level2}, index=barcodes),
        var=pd.DataFrame(index=genes),
    )
    adata.obsm["X_umap_local"] = umap_local
    print(f"    {n_cells:,} cells x {len(genes):,} genes loaded in {time.time()-t0:.0f}s")
    return adata


def preprocess_for_nmf(adata):
    print("  [PREPROCESS] Normalize + HVG selection...", flush=True)
    t0 = time.time()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat_v3", layer="counts")
    n_hvg = adata.var["highly_variable"].sum()
    print(f"    {n_hvg:,} HVGs selected")
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
    print(f"  [NMF] k={k}, training on {N_TRAIN:,} cells...", flush=True)
    t0 = time.time()
    np.random.seed(SEED)
    rng = np.random.RandomState(RNG_SEED)
    level2 = adata.obs["celltype_level2"].values

    train_idx = []
    for st in np.unique(level2):
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
    model = NMF(n_components=k, init="nndsvda", random_state=RNG_SEED, max_iter=500)
    W_train = model.fit_transform(X_train)
    H = model.components_
    print(f"    NMF fit in {time.time()-t_fit:.0f}s (recon error: {model.reconstruction_err_:.2f})")
    del X_train, W_train; gc.collect()

    print("    Projecting all cells...", flush=True)
    t_proj = time.time()
    X_all = adata[:, hvg_mask].X
    if sparse.issparse(X_all):
        X_all = X_all.toarray()
    X_all = np.clip(X_all, 0, None)
    W_all = model.transform(X_all)
    print(f"    Projection done in {time.time()-t_proj:.0f}s")
    del X_all; gc.collect()

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
    print("  [IDENTIFY] Scoring factors for SecA/SecB enrichment...")
    scores = []
    for factor in loadings_df.index:
        loadings = loadings_df.loc[factor]
        ranked = loadings.sort_values(ascending=False)
        top100 = set(ranked.index[:100]); top200 = set(ranked.index[:200])
        factor_max = loadings.max()
        scores.append({
            "factor": factor,
            "secA_top100": len(set(SECA_GENES) & top100),
            "secA_top200": len(set(SECA_GENES) & top200),
            "secB_top100": len(set(SECB_GENES) & top100),
            "secB_top200": len(set(SECB_GENES) & top200),
            "secA_loading": loadings[loadings.index.isin(SECA_GENES)].mean() / factor_max,
            "secB_loading": loadings[loadings.index.isin(SECB_GENES)].mean() / factor_max,
        })
    scores_df = pd.DataFrame(scores).set_index("factor")
    best_secA = scores_df["secA_loading"].idxmax()
    best_secB = scores_df["secB_loading"].idxmax()
    print(f"    Best SecA factor: {best_secA} ({scores_df.loc[best_secA, 'secA_top100']} genes in top100)")
    print(f"    Best SecB factor: {best_secB} ({scores_df.loc[best_secB, 'secB_top100']} genes in top100)")
    return scores_df, best_secA, best_secB


# ============================================================================
# SUPP DATA 4 — tidy reshape
# ============================================================================
def write_supp_data_4(loadings_df):
    """Wide H matrix (factor x gene) -> tidy long table for Supp Data 4."""
    print("  [SD4] Building tidy NMF loadings table...")
    tidy = (loadings_df.reset_index().rename(columns={"index": "factor"})
            .melt(id_vars="factor", var_name="gene", value_name="loading"))
    # within-factor descending rank of each gene's loading
    tidy["rank_in_factor"] = (tidy.groupby("factor")["loading"]
                              .rank(ascending=False, method="first").astype(int))
    tidy["is_SecA_gene"] = tidy["gene"].isin(SECA_GENES)
    tidy["is_SecB_gene"] = tidy["gene"].isin(SECB_GENES)
    tidy = tidy.sort_values(["factor", "rank_in_factor"]).reset_index(drop=True)
    out = path("output_root", "supplemental", "supplementary_table_NMF_factor_genes.csv")
    tidy.to_csv(out, index=False)
    print(f"    Saved: {out}  ({len(tidy):,} rows)")


# ============================================================================
# FIGURES  (diagnostic; logic unchanged)
# ============================================================================
def plot_factor_umaps(usage_df, umap_coords, k, out_dir):
    print("  [FIG] Factor usage UMAPs...", flush=True)
    ncols = 4; nrows = (k + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.8 * nrows))
    axes = axes.flatten()
    rng = np.random.RandomState(SEED)
    n = len(usage_df); idx = rng.choice(n, size=min(150_000, n), replace=False)
    for i in range(k):
        ax = axes[i]; vals = usage_df[f"Factor_{i+1}"].values[idx]
        sp = ax.scatter(umap_coords[idx, 0], umap_coords[idx, 1], c=vals, cmap="YlOrRd",
                        s=0.3, alpha=0.6, edgecolors="none", rasterized=True,
                        vmin=np.percentile(vals, 1), vmax=np.percentile(vals, 99))
        theme_void(ax); ax.set_title(f"Factor_{i+1}", fontsize=8, fontweight="bold")
        plt.colorbar(sp, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    for j in range(k, len(axes)):
        axes[j].set_visible(False)
    fig.suptitle("NMF Factor Usage on Epithelial UMAP", fontsize=10, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    save_fig(fig, "11d_A_nmf_factor_umaps", out_dir)


def plot_secA_secB_umap(usage_df, umap_coords, best_secA, best_secB, out_dir):
    print("  [FIG] SecA/SecB UMAP comparison...", flush=True)
    rng = np.random.RandomState(SEED)
    n = len(usage_df); idx = rng.choice(n, size=min(150_000, n), replace=False)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, factor, title, cmap in [
        (axes[0], best_secA, "SecA (progenitor-like)", "YlOrBr"),
        (axes[1], best_secB, "SecB (adaptive)", "OrRd"),
    ]:
        vals = usage_df[factor].values[idx]
        sp = ax.scatter(umap_coords[idx, 0], umap_coords[idx, 1], c=vals, cmap=cmap,
                        s=0.4, alpha=0.6, edgecolors="none", rasterized=True,
                        vmin=np.percentile(vals, 2), vmax=np.percentile(vals, 98))
        theme_void(ax); ax.set_title(title, fontsize=8, fontweight="bold")
        plt.colorbar(sp, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    fig.suptitle("NMF Factor Usage - SecA/SecB Polarization", fontsize=10, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    save_fig(fig, "11d_E_secA_secB_umap", out_dir)


def plot_polarization(usage_df, level2, best_secA, best_secB, out_dir):
    print("  [FIG] SecA/SecB polarization scatter...", flush=True)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    rng = np.random.RandomState(SEED)
    subtypes = level2.values
    for st in pd.Series(subtypes).value_counts().index:
        idx = np.where(subtypes == st)[0]
        idx_sub = rng.choice(idx, size=min(len(idx), 20_000), replace=False)
        ax.scatter(usage_df[best_secA].values[idx_sub], usage_df[best_secB].values[idx_sub],
                   c=[EPI_PALETTE.get(st, "#999999")], s=0.8, alpha=0.5,
                   edgecolors="none", rasterized=True, label=st)
    ax.set_xlabel(f"{best_secA} usage (SecA / progenitor-like)", fontsize=8)
    ax.set_ylabel(f"{best_secB} usage (SecB / adaptive)", fontsize=8)
    ax.set_title("Epithelial SecA <-> SecB Polarization", fontsize=9, fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=EPI_PALETTE[st],
                      markeredgecolor="none", markersize=5, label=st)
               for st in sorted(EPI_PALETTE.keys())]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              frameon=False, fontsize=5.5, handletextpad=0.3, labelspacing=0.5)
    plt.tight_layout()
    save_fig(fig, "11d_B_secA_secB_polarization", out_dir)


def plot_gene_loadings(loadings_df, best_secA, best_secB, out_dir):
    print("  [FIG] Gene loading barplots...", flush=True)
    fig, axes = plt.subplots(1, 2, figsize=(9, 5))
    for ax, factor, sig_genes, sig_name in [
        (axes[0], best_secA, SECA_GENES, "SecA (progenitor)"),
        (axes[1], best_secB, SECB_GENES, "SecB (adaptive)"),
    ]:
        top30 = loadings_df.loc[factor].sort_values(ascending=False).head(30)
        colors = ["#D14E6C" if g in sig_genes else "#888888" for g in top30.index]
        ax.barh(range(len(top30)), top30.values[::-1], color=colors[::-1], edgecolor="none", height=0.7)
        ax.set_yticks(range(len(top30))); ax.set_yticklabels(top30.index[::-1], fontsize=5.5)
        ax.set_xlabel("Loading weight", fontsize=7)
        ax.set_title(f"{factor} - {sig_name}\nTop 30 genes", fontsize=8, fontweight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.text(0.95, 0.05, f"Red = {sig_name} signature", transform=ax.transAxes,
                fontsize=5.5, ha="right", color="#D14E6C", fontweight="bold")
    plt.tight_layout()
    save_fig(fig, "11d_C_gene_loadings", out_dir)


def plot_factor_by_subtype(usage_df, level2, best_secA, best_secB, out_dir):
    print("  [FIG] Factor usage by subtype (violin)...", flush=True)
    import seaborn as sns
    rng = np.random.RandomState(SEED)
    n = len(usage_df); idx = rng.choice(n, size=min(100_000, n), replace=False)
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
        sns.violinplot(data=plot_df, y="Subtype", x=col, order=order, palette=palette, ax=ax,
                       cut=0, linewidth=0.5, inner="quartile", density_norm="width")
        ax.set_title(title, fontsize=8, fontweight="bold")
        ax.set_xlabel("Factor usage (normalized)", fontsize=7); ax.set_ylabel("")
        ax.tick_params(axis="y", labelsize=5.5)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.tight_layout()
    save_fig(fig, "11d_D_factor_by_subtype", out_dir)


# ============================================================================
# MAIN
# ============================================================================
def main():
    parser = argparse.ArgumentParser(description="Supp Data 4: Epithelial NMF gene loadings")
    parser.add_argument("--k", type=int, default=DEFAULT_K, help=f"NMF factors (default {DEFAULT_K})")
    parser.add_argument("--figures-only", action="store_true",
                        help="Reload saved NMF CSVs, regenerate SD4 table + figures")
    args = parser.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()
    print("=" * 65)
    print("  Supp Data 4 — Epithelial NMF (SecA/SecB loadings)")
    print("=" * 65)

    usage_path = os.path.join(OUT_DIR, "11d_nmf_usage.csv")
    loadings_path = os.path.join(OUT_DIR, "11d_nmf_loadings.csv")

    if args.figures_only and os.path.exists(usage_path) and os.path.exists(loadings_path):
        print("\n  [LOAD] Reading saved NMF results...")
        usage_df = pd.read_csv(usage_path, index_col=0)
        loadings_df = pd.read_csv(loadings_path, index_col=0)
        k = len(loadings_df)
        with h5py.File(H5AD, "r") as f:
            umap_coords = np.array(f["obsm"]["X_umap_local"])
            raw_level2 = [x.decode() if isinstance(x, bytes) else x
                          for x in f["obs"]["celltype_level2"][:]]
            idx = [x.decode() if isinstance(x, bytes) else x for x in f["obs"]["_index"][:]]
        level2 = pd.Series(_relabel_intermediate(raw_level2), index=idx)
    else:
        adata = load_epithelial()
        adata = preprocess_for_nmf(adata)
        usage_df, usage_raw, loadings_df, _ = run_nmf(adata, k=args.k)
        k = args.k
        umap_coords = adata.obsm["X_umap_local"]
        level2 = adata.obs["celltype_level2"]
        print("\n  [SAVE] Writing NMF results...")
        usage_df.to_csv(usage_path)
        usage_raw.to_csv(os.path.join(OUT_DIR, "11d_nmf_usage_raw.csv"))
        loadings_df.to_csv(loadings_path)
        del adata; gc.collect()

    scores_df, best_secA, best_secB = identify_polarization_factors(loadings_df)
    scores_df.to_csv(os.path.join(OUT_DIR, "11d_factor_scores.csv"))

    # Supp Data 4 tidy export
    write_supp_data_4(loadings_df)

    print(f"\n  [FIGURES] Generating to {OUT_DIR}")
    plot_factor_umaps(usage_df, umap_coords, k, OUT_DIR)
    plot_secA_secB_umap(usage_df, umap_coords, best_secA, best_secB, OUT_DIR)
    plot_polarization(usage_df, level2, best_secA, best_secB, OUT_DIR)
    plot_gene_loadings(loadings_df, best_secA, best_secB, OUT_DIR)
    plot_factor_by_subtype(usage_df, level2, best_secA, best_secB, OUT_DIR)

    print(f"\n{'=' * 65}")
    print(f"  Supp Data 4 COMPLETE ({time.time()-t0:.0f}s)")
    print(f"  NMF: k={k}, best SecA={best_secA}, best SecB={best_secB}")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
