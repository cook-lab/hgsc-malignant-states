#!/usr/bin/env python3
"""
Cell-cell communication with NMF epitype labels (LIANA+)
========================================================
HGSC malignant-states atlas backend.

Runs LIANA+ rank_aggregate on the whole atlas using the NMF-derived epithelial
labels (`celltype_nmf`: SecA / Intermediate / SecB / Ciliated epithelium) so
ligand-receptor interactions are resolved at epitype resolution. Global
analysis only (no treatment stratification). Memory-frugal: row-selective h5py
loading of the subsampled cells.

INPUTS:
  - hgsc_atlas_final.h5ad   [config: atlas_final]   (obs + layers/counts via h5py)
  - output_root/03_epithelial_nmf/celltype_nmf_mapping.csv  (from 02_prepare_nmf_labels.py)

OUTPUTS (output_root/06_cellcomm/):
  - tables/17b_liana_global.csv  (KEY cache — Fig 5F, Supp Data 7)
  - tables/17b_liana_global_top50.csv, 17b_interaction_counts_global.csv,
    17b_subsampling_global.csv
  - figs/17b_*.svg/pdf

MANUSCRIPT PANELS: Fig 5F (autocrine), Supp Data 7 (autocrine L-R pairs).

RUNTIME TIER: heavy (LIANA permutations; hours).

SEEDING: stratified subsample uses RANDOM_SEED = config SEED; the LIANA
  permutation seed is a fixed analytical parameter (preserved as published).

Usage:
    python 01_cellcomm_nmf.py
    python 01_cellcomm_nmf.py --skip-figures
    python 01_cellcomm_nmf.py --n-perms 50    # fast test
"""

import argparse
import gc
import os
import sys
import time
import warnings
from datetime import datetime

import numpy as np
import pandas as pd
import anndata as ad
import h5py
from scipy.sparse import csr_matrix

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path, SEED  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")

try:
    import liana as li
except ImportError:
    sys.exit("ERROR: liana not installed.  Run: pip install liana>=1.6.0")

# ── Cook Lab style v1.2 ─────────────────────────────────────
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.dpi":      300,
    "savefig.dpi":     450,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
    "savefig.bbox":    "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# ── paths ────────────────────────────────────────────────────
ATLAS_H5AD  = obj("atlas_final")
NMF_MAPPING = path("output_root", "03_epithelial_nmf", "celltype_nmf_mapping.csv")
OUT_DIR     = path("output_root", "06_cellcomm")
FIG_DIR     = os.path.join(OUT_DIR, "figs")
TABLE_DIR   = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# ── parameters ───────────────────────────────────────────────
GROUPBY = "celltype_nmf"          # NMF-derived labels
MAX_CELLS_PER_TYPE = 200
MIN_CELLS_PER_TYPE = 10
N_PERMS = 100
N_JOBS  = 1                       # sequential: avoids joblib deadlock on macOS
RANDOM_SEED = SEED                # subsample reproducibility (from config)

# ── palettes ─────────────────────────────────────────────────
CELLTYPE_PALETTE = {
    "Epithelial":    "#E6A141", "Mesothelial":   "#A8A298",
    "Fibroblast":    "#DDD5CA", "Smooth muscle": "#D14E6C",
    "Pericyte":      "#B87A7A", "Endothelial":   "#7D4E4E",
    "T/NK cell":     "#87CEFA", "B cell":        "#5665B6",
    "Plasma cell":   "#8A5DAF", "Macrophage":    "#8FBC8F",
    "DC":            "#2E8B57", "Neutrophil":    "#6B8E23",
    "Mast cell":     "#8B9B6B",
}
CELLTYPE_ORDER = [
    "Epithelial", "Mesothelial", "Fibroblast", "Smooth muscle",
    "Pericyte", "Endothelial", "T/NK cell", "B cell", "Plasma cell",
    "Macrophage", "DC", "Neutrophil", "Mast cell",
]

# NMF epithelial labels ("Transitioning" -> "Intermediate")
NMF_EPI_LABELS = {"SecA epithelium", "Intermediate epithelium",
                  "SecB epithelium", "Ciliated epithelial cell"}


# ============================================================================
# HELPERS
# ============================================================================

def get_subsample_indices(obs_df, groupby, max_per_type, min_per_type, seed):
    rng = np.random.RandomState(seed)
    labels = obs_df[groupby].values
    unique_labels = [l for l in sorted(np.unique(labels)) if pd.notna(l)]
    keep_idx, summary = [], []
    for label in unique_labels:
        idx = np.where(labels == label)[0]
        n_orig = len(idx)
        if n_orig < min_per_type:
            summary.append((label, n_orig, 0, "SKIP (<min)"))
            continue
        n_keep = min(n_orig, max_per_type)
        keep_idx.append(rng.choice(idx, size=n_keep, replace=False))
        summary.append((label, n_orig, n_keep, "OK"))
    if keep_idx:
        keep_idx = np.concatenate(keep_idx)
        rng.shuffle(keep_idx)
    else:
        keep_idx = np.array([], dtype=int)
    summary_df = pd.DataFrame(summary,
                              columns=["celltype", "n_original", "n_subsampled", "status"])
    return keep_idx, summary_df


def load_counts_for_cells(h5_path, row_indices, n_cols):
    """Row-selective CSR loader from h5py. Peak memory ~100-300 MB."""
    row_indices = np.sort(np.asarray(row_indices))
    n_rows = len(row_indices)
    print(f"   Target: {n_rows:,} cells (row-selective h5py read)")
    t0 = time.time()
    with h5py.File(h5_path, "r") as f:
        grp = f["layers/counts"]
        indptr_full = grp["indptr"][:]
        row_starts = indptr_full[row_indices]
        row_ends = indptr_full[row_indices + 1]
        row_lengths = row_ends - row_starts
        total_nnz = int(row_lengths.sum())
        del indptr_full
        print(f"   Need {total_nnz:,} non-zero entries across {n_rows:,} rows")

        new_indptr = np.zeros(n_rows + 1, dtype=np.int64)
        np.cumsum(row_lengths, out=new_indptr[1:])
        data_dtype = grp["data"].dtype
        new_data = np.empty(total_nnz, dtype=data_dtype)
        new_indices = np.empty(total_nnz, dtype=np.int32)
        data_ds = grp["data"]
        indices_ds = grp["indices"]

        out_pos = i = n_batches = 0
        while i < n_rows:
            j = i + 1
            while j < n_rows and row_indices[j] == row_indices[j - 1] + 1:
                j += 1
            read_start = int(row_starts[i])
            read_end = int(row_ends[j - 1])
            read_len = read_end - read_start
            if read_len > 0:
                new_data[out_pos:out_pos + read_len] = data_ds[read_start:read_end]
                chunk = indices_ds[read_start:read_end]
                new_indices[out_pos:out_pos + read_len] = (
                    chunk.astype(np.int32) if chunk.dtype != np.int32 else chunk)
                out_pos += read_len
            n_batches += 1
            i = j
            if n_batches % 500 == 0:
                print(f"      {i/n_rows*100:.0f}% ({i:,}/{n_rows:,} rows)")
    counts = csr_matrix((new_data, new_indices, new_indptr), shape=(n_rows, n_cols))
    print(f"   Done: {counts.shape}, {counts.nnz:,} nnz ({time.time()-t0:.0f}s)")
    return counts


def build_adata(counts_matrix, obs_df, var_names):
    adata = ad.AnnData(X=counts_matrix.copy(), obs=obs_df.reset_index(drop=True),
                       var=pd.DataFrame(index=var_names))
    adata.obs.index = [f"cell_{i}" for i in range(adata.n_obs)]
    return adata


def save_table(df, name):
    p = os.path.join(TABLE_DIR, f"{name}.csv")
    df.to_csv(p, index=False)
    print(f"   Saved: {os.path.basename(p)} ({len(df):,} rows)")
    return p


def save_fig(fig, name):
    for fmt in ("svg", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{fmt}"), format=fmt, dpi=450,
                    bbox_inches="tight")
    plt.close(fig)
    print(f"   Saved: {name}.svg/.pdf")


# ============================================================================
# VISUALIZATION
# ============================================================================

def plot_dotplot(liana_df, top_n=30, title="Top L-R Interactions", name="dotplot"):
    top = liana_df.nsmallest(top_n, "magnitude_rank").copy()
    top["lr_pair"] = top["ligand_complex"].astype(str) + " -> " + top["receptor_complex"].astype(str)
    fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.25)))
    scatter = ax.scatter(x=top["specificity_rank"], y=range(len(top)),
                         s=(1 - top["magnitude_rank"]) * 200 + 10,
                         c=top["magnitude_rank"], cmap="YlOrRd_r",
                         edgecolors="black", linewidths=0.3, alpha=0.8)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(
        [f"{row.lr_pair}  ({row.source} -> {row.target})" for _, row in top.iterrows()],
        fontsize=5)
    ax.set_xlabel("Specificity Rank (lower = more specific)")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.invert_yaxis()
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Magnitude Rank", fontsize=7)
    fig.tight_layout()
    save_fig(fig, name)


def plot_interaction_heatmap(liana_df, level1_map, sig_threshold=0.05,
                             title="Cell-Cell Interaction Counts", name="heatmap"):
    sig = liana_df[liana_df["magnitude_rank"] <= sig_threshold].copy()
    if len(sig) == 0:
        print(f"   WARNING: No significant interactions at threshold {sig_threshold}")
        return
    counts = sig.groupby(["source", "target"]).size().unstack(fill_value=0)
    all_types = sorted(set(counts.index) | set(counts.columns))
    type_order = sorted(all_types, key=lambda x: (
        CELLTYPE_ORDER.index(level1_map.get(x, "Unknown"))
        if level1_map.get(x, "Unknown") in CELLTYPE_ORDER else 99, x))
    counts = counts.reindex(index=[t for t in type_order if t in counts.index],
                            columns=[t for t in type_order if t in counts.columns],
                            fill_value=0)
    fig, ax = plt.subplots(figsize=(16, 14))
    sns.heatmap(counts, cmap="YlOrRd", ax=ax, linewidths=0.3, linecolor="white",
                cbar_kws={"label": f"Significant interactions (mag_rank <= {sig_threshold})",
                          "shrink": 0.5})
    ax.set_title(title, fontsize=10, fontweight="bold", pad=10)
    ax.set_xlabel("Receiver (target)", fontsize=8)
    ax.set_ylabel("Sender (source)", fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=4.5)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=4.5)
    fig.tight_layout()
    save_fig(fig, name)


def plot_level1_chord(liana_df, level1_map, sig_threshold=0.05,
                      title="Level-1 Interaction Summary", name="chord_level1"):
    sig = liana_df[liana_df["magnitude_rank"] <= sig_threshold].copy()
    if len(sig) == 0:
        return
    sig["source_l1"] = sig["source"].map(level1_map).fillna("Other")
    sig["target_l1"] = sig["target"].map(level1_map).fillna("Other")
    counts = sig.groupby(["source_l1", "target_l1"]).size().unstack(fill_value=0)
    order = [c for c in CELLTYPE_ORDER if c in counts.index or c in counts.columns]
    counts = counts.reindex(index=[c for c in order if c in counts.index],
                            columns=[c for c in order if c in counts.columns],
                            fill_value=0)
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(counts, cmap="YlOrRd", ax=ax, annot=True, fmt="d",
                linewidths=0.5, linecolor="white",
                cbar_kws={"label": "Significant L-R pairs", "shrink": 0.6})
    for lbl in ax.get_yticklabels() + ax.get_xticklabels():
        lbl.set_color(CELLTYPE_PALETTE.get(lbl.get_text(), "#333333"))
        lbl.set_fontweight("bold")
    ax.set_title(title, fontsize=10, fontweight="bold", pad=10)
    ax.set_xlabel("Receiver compartment", fontsize=8)
    ax.set_ylabel("Sender compartment", fontsize=8)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=7)
    fig.tight_layout()
    save_fig(fig, name)


# ============================================================================
# MAIN
# ============================================================================

def main():
    global MAX_CELLS_PER_TYPE
    parser = argparse.ArgumentParser(
        description="Cell Communication with NMF labels (LIANA+)")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--max-cells", type=int, default=MAX_CELLS_PER_TYPE)
    parser.add_argument("--n-perms", type=int, default=N_PERMS)
    parser.add_argument("--n-jobs", type=int, default=N_JOBS)
    args = parser.parse_args()

    MAX_CELLS_PER_TYPE = args.max_cells
    n_perms, n_jobs = args.n_perms, args.n_jobs

    print("=" * 70)
    print("  Cell-Cell Communication (LIANA+) with NMF labels")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  LIANA v{li.__version__}")
    print(f"  Groupby: {GROUPBY}  Perms: {n_perms}  Max cells/type: {MAX_CELLS_PER_TYPE}")
    print("=" * 70)

    # ── PHASE 1: metadata ─────────────────────────────────────
    print("\nPHASE 1: Loading atlas metadata (backed mode)")
    adata_backed = ad.read_h5ad(ATLAS_H5AD, backed="r")
    obs_full = adata_backed.obs[["celltype_level1", "celltype_level2"]].copy()
    var_names = adata_backed.var_names.tolist()
    n_cells_total, n_genes = adata_backed.shape
    print(f"   Atlas: {n_cells_total:,} cells x {n_genes:,} genes")
    del adata_backed
    gc.collect()

    if not os.path.exists(NMF_MAPPING):
        sys.exit(f"ERROR: NMF mapping not found at {NMF_MAPPING}.\n"
                 f"Run 02_prepare_nmf_labels.py first.")
    nmf_map = pd.read_csv(NMF_MAPPING, index_col=0)
    # Barcode (index) join, NOT positional .values: atlas_final.obs and the NMF
    # mapping are both barcode-indexed over the full atlas, but a positional assign
    # silently mis-wires cell types if their row order ever diverges. reindex aligns
    # on barcode and the assertion fails loudly on any mismatch (audit follow-up).
    obs_full[GROUPBY] = nmf_map["celltype_nmf"].reindex(obs_full.index)
    n_missing = int(obs_full[GROUPBY].isna().sum())
    assert n_missing == 0, (
        f"{n_missing} atlas cells have no celltype_nmf in {os.path.basename(NMF_MAPPING)} "
        f"— barcode mismatch between atlas_final.obs and the NMF mapping.")
    print(f"   Loaded NMF labels (barcode join); unique {GROUPBY}: {obs_full[GROUPBY].nunique()}")

    level1_map = {}
    for nmf_label, l1 in zip(obs_full[GROUPBY], obs_full["celltype_level1"]):
        if nmf_label not in level1_map:
            level1_map[nmf_label] = "Epithelial" if nmf_label in NMF_EPI_LABELS else l1

    # ── PHASE 2: global LIANA ─────────────────────────────────
    print("\nPHASE 2: Global LIANA analysis (NMF labels)")
    global_indices = np.arange(len(obs_full))
    sub_pos, sub_summary = get_subsample_indices(
        obs_full, GROUPBY, MAX_CELLS_PER_TYPE, MIN_CELLS_PER_TYPE, RANDOM_SEED)
    save_table(sub_summary, "17b_subsampling_global")

    n_types = (sub_summary["status"] == "OK").sum()
    n_cells = int(sub_summary["n_subsampled"].sum())
    print(f"   Subsample: {n_cells:,} cells, {n_types} cell types")
    if n_cells == 0:
        sys.exit("ERROR: No cells after subsampling")

    atlas_sub = global_indices[sub_pos]
    counts_matrix = load_counts_for_cells(ATLAS_H5AD, atlas_sub, n_genes)
    adata = build_adata(counts_matrix, obs_full.iloc[atlas_sub], var_names)
    del counts_matrix
    gc.collect()

    print(f"\n   Running LIANA (global, NMF labels)... "
          f"Cells: {adata.n_obs:,}, Groups: {adata.obs[GROUPBY].nunique()}")
    t0 = time.time()
    li.mt.rank_aggregate(
        adata, groupby=GROUPBY, resource_name="consensus", expr_prop=0.1,
        min_cells=MIN_CELLS_PER_TYPE, use_raw=False, layer=None, n_perms=n_perms,
        seed=1337, n_jobs=n_jobs, key_added="liana_res", verbose=True, inplace=True)
    results = adata.uns["liana_res"].copy()
    print(f"   Done in {time.time()-t0:.0f}s — {len(results):,} interaction rows")

    save_table(results, "17b_liana_global")
    save_table(results.nsmallest(50, "magnitude_rank"), "17b_liana_global_top50")

    sig = results[results["magnitude_rank"] <= 0.05]
    if len(sig) > 0:
        cnt = sig.groupby(["source", "target"]).size().unstack(fill_value=0)
        cnt.to_csv(os.path.join(TABLE_DIR, "17b_interaction_counts_global.csv"))
        print(f"   Saved: 17b_interaction_counts_global.csv ({len(sig):,} sig interactions)")

    del adata
    gc.collect()

    # ── PHASE 3: figures ──────────────────────────────────────
    if not args.skip_figures:
        print("\nPHASE 3: Generating figures")
        plot_dotplot(results, top_n=30,
                     title="HGSC Atlas — Top 30 L-R Interactions (NMF labels)",
                     name="17b_dotplot_global_top30")
        plot_interaction_heatmap(results, level1_map,
                                 title="Cell-Cell Interaction Counts (NMF labels, mag_rank <= 0.05)",
                                 name="17b_heatmap_interaction_counts_global")
        plot_level1_chord(results, level1_map,
                          title="Compartment-Level Interaction Summary (NMF labels)",
                          name="17b_heatmap_level1_global")

    print(f"\n{'='*70}")
    print("  DONE — Cell Communication (LIANA+, NMF labels)")
    print(f"  Output:  {OUT_DIR}")
    print(f"  Global:  {len(results):,} interactions; significant: {len(sig):,}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
