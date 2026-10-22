#!/usr/bin/env python3
"""
Figure S6 — SecA / SecB marker-gene expression on the epithelial UMAP
=====================================================================

Purpose
    NMF factor usage + the 7 SecA and 7 SecB marker-gene expression UMAPs on the
    epithelial-only UMAP (n=575,366), as two 2 x 4 panels:
      a: SecA factor usage (Factor 3) + the 7 SecA genes
      b: SecB factor usage (Factor 2) + the 7 SecB genes
    Each panel is one half of a 4 x 4 portrait-letter-page grid (PAGE_W x PAGE_H).

INPUTS
    figures._prep.atlas_final_tables.load_epithelial()  (575,366 epithelial cells: epithelial
        UMAP1/UMAP2 + Factor_3 = SecA, Factor_2 = SecB usage).
    data_root/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad
        (log-normalised X, chunked h5py scan; gene expression ONLY — its obs metadata are not
        used; cells matched by barcode)
    shared/signatures.yml  (canonical 7-gene noBCAM SecA / SecB lists)

OUTPUTS
    output_root/figures/supplementary/SF6a_atlas_SecA_signature_umaps.{svg,png}    (panel a)
    output_root/figures/supplementary/SF6b_atlas_SecB_signature_umaps.{svg,png}    (panel b)

MANUSCRIPT PANEL(S)
    Figure S6a, S6b.

RUNTIME TIER
    moderate (chunked expression extraction for 14 markers + 16 rasterized scatters).
"""

import os
import sys
import gc

import numpy as np
import pandas as pd
import anndata as ad
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO_ROOT)
from config.config import path, SEED  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

H5AD = path("data_root", "2026_final_atlas", "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
OUT_STANDALONE = path("output_root", "figures", "supplementary", "SF6{part}_atlas_Sec{sig}_signature_umaps.{ext}")

assert os.path.exists(H5AD), f"Missing: {H5AD}"

# ============================================================================
# MARKER GENES — loaded from shared/signatures.yml (canonical noBCAM 7-gene set)
# (Do NOT inline divergent lists; convention #6.)
# ============================================================================

with open(os.path.join(REPO_ROOT, "shared", "signatures.yml")) as fh:
    SIGS = yaml.safe_load(fh)
SECA_GENES = list(SIGS["SecA"])
SECB_GENES = list(SIGS["SecB"])
ALL_GENES = SECA_GENES + SECB_GENES

# ============================================================================
# STYLE
# ============================================================================

FA, FK, FN, FG = 6, 5.5, 5, 6

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        FA,
    "axes.labelsize":   FA,
    "axes.titlesize":   0,
    "xtick.labelsize":  FK,
    "ytick.labelsize":  FK,
    "legend.fontsize":  FN,
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,
    "figure.dpi":       450,
    "savefig.dpi":      450,
})

CMAP = LinearSegmentedColormap.from_list(
    "beige_burgundy",
    [(0.00, "#DDD5CA"), (0.48, "#DDD5CA"), (0.65, "#C8878A"),
     (0.85, "#A04858"), (1.00, "#6B2838")],
    N=256,
)

# ============================================================================
# LOAD DATA
# ============================================================================

print("\nLoading epithelial table (load_epithelial)...", flush=True)
meta = load_epithelial(columns=["UMAP1", "UMAP2", "Factor_2", "Factor_3"])
print(f"  {len(meta):,} cells")

# Chunked h5py scan of the CSR X (log-normalised) for the 14 markers, over ALL
# cells (the full matrix is never read into memory).
import h5py  # noqa: E402
print("Reading X (chunked CSR scan, all cells)...", flush=True)
with h5py.File(H5AD, "r") as fh5:
    var_names = pd.Index(ad.io.read_elem(fh5["var"]).index.astype(str))
    h5_index = pd.Index(ad.io.read_elem(fh5["obs"][fh5["obs"].attrs["_index"]]).astype(str))
    missing = [g for g in ALL_GENES if g not in var_names]
    assert not missing, f"Genes not found in epithelial h5ad: {missing}"
    colmap = np.full(len(var_names), -1, dtype=np.int64)
    colmap[[var_names.get_loc(g) for g in ALL_GENES]] = np.arange(len(ALL_GENES))
    Xg = fh5["X"]
    indptr = Xg["indptr"][:]
    n_obs = len(indptr) - 1
    dense = np.zeros((n_obs, len(ALL_GENES)), dtype=np.float32)
    CH = 200_000_000
    r0 = 0
    while r0 < n_obs:
        s = indptr[r0]
        r1 = min(max(int(np.searchsorted(indptr, s + CH, side="right") - 1), r0 + 1), n_obs)
        e = indptr[r1]
        c = colmap[Xg["indices"][s:e]]
        m = c >= 0
        if m.any():
            dat = Xg["data"][s:e][m]
            rows = np.repeat(np.arange(r0, r1), np.diff(indptr[r0:r1 + 1]))[m]
            dense[rows, c[m]] = dat
        print(f"  rows {r0:,}-{r1:,} / {n_obs:,}", flush=True)
        r0 = r1
print(f"  h5ad shape: ({n_obs}, {len(var_names)})")
meta_idx_in_h5 = h5_index.get_indexer(meta.index)
assert (meta_idx_in_h5 == -1).sum() == 0, "Unmatched cells"
expr = dense[meta_idx_in_h5]   # order = meta rows, columns = ALL_GENES
print(f"  Expression matrix: {expr.shape}")
del dense
gc.collect()

secA_usage = meta["Factor_3"].to_numpy(dtype=np.float32)
secB_usage = meta["Factor_2"].to_numpy(dtype=np.float32)

x = meta["UMAP1"].to_numpy()
y = meta["UMAP2"].to_numpy()
xr = x.max() - x.min()
yr = y.max() - y.min()
pad_x = 0.04 * xr
pad_y = 0.04 * yr

SHUFFLE_ORDER = np.random.default_rng(seed=SEED).permutation(len(meta))
x_shuf = x[SHUFFLE_ORDER]
y_shuf = y[SHUFFLE_ORDER]

# ============================================================================
# GRID DEFINITION — (label, values, is_gene_italic)
# ============================================================================

panels_secA = [("Factor 3 (SecA)", secA_usage, False)]
for g in SECA_GENES:
    panels_secA.append((g, expr[:, ALL_GENES.index(g)], True))

panels_secB = [("Factor 2 (SecB)", secB_usage, False)]
for g in SECB_GENES:
    panels_secB.append((g, expr[:, ALL_GENES.index(g)], True))


def draw_umap(fig, ax, label, vals, is_gene):
    """One all-cell (n = 575,366) rasterized UMAP panel coloured by vals."""
    vo = vals[SHUFFLE_ORDER]
    nz = vals[vals > 0]
    vmax = float(np.quantile(nz, 0.99)) if nz.size else 1.0
    if vmax <= 0:
        vmax = 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)

    ax.scatter(x_shuf, y_shuf, c=vo, cmap=CMAP, norm=norm, s=0.02, marker="o",
               linewidths=0, rasterized=True, alpha=0.95)

    ax.set_xlim(x.min() - pad_x, x.max() + pad_x)
    ax.set_ylim(y.min() - pad_y, y.max() + pad_y)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_bounds(x.min() - pad_x, x.min() - pad_x + 0.18 * xr)
    ax.spines["left"].set_bounds(y.min() - pad_y, y.min() - pad_y + 0.18 * yr)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.spines["left"].set_linewidth(0.4)

    ax.text(0.03, 0.97, label, transform=ax.transAxes, fontsize=FG,
            fontstyle="italic" if is_gene else "normal",
            fontweight="normal" if is_gene else "bold", ha="left", va="top")

    cax = ax.inset_axes([0.88, 0.05, 0.04, 0.3])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="vertical")
    cb.outline.set_linewidth(0.3)
    cb.ax.tick_params(length=1.5, width=0.3, pad=1, labelsize=FN - 1)
    cb.set_ticks([0.0, vmax])
    fmt = "{:.1f}" if is_gene else "{:.2f}"
    cb.set_ticklabels([fmt.format(0.0), fmt.format(vmax)])


# ============================================================================
# PAGE GEOMETRY — a 4 rows × 4 cols portrait letter page; each panel is one half of it
# ============================================================================

print("Plotting...", flush=True)

NROWS, NCOLS = 4, 4
PAGE_W = 7.5
# page height from the UMAP aspect ratio
PAGE_H = NROWS * (PAGE_W / NCOLS) * (yr + 2 * pad_y) / (xr + 2 * pad_x) * 1.08

# ============================================================================
# PANELS — SF6a = SecA block (Figure S6a), SF6b = SecB block (Figure S6b); 2 x 4 each,
# no panel letters / section labels.
# ============================================================================

for part, sig, panels in (("a", "A", panels_secA), ("b", "B", panels_secB)):
    fig, axes = plt.subplots(2, 4, figsize=(PAGE_W, PAGE_H / 2),
                             gridspec_kw={"wspace": 0.08, "hspace": 0.12})
    for i, pnl in enumerate(panels):
        draw_umap(fig, axes[i // 4, i % 4], *pnl)
    for ext in ("svg", "png"):
        out = OUT_STANDALONE.format(part=part, sig=sig, ext=ext)
        fig.savefig(out, format=ext, dpi=450, bbox_inches="tight", pad_inches=0.06)
        print(f"  Saved: {out}")
    plt.close(fig)
