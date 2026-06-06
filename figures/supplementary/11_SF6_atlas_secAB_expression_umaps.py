#!/usr/bin/env python3
"""
SF6 — SecA / SecB marker-gene expression on the epithelial UMAP
===============================================================

Purpose
    Portrait letter page combining NMF factor usage + the 7 SecA and 7 SecB
    marker-gene expression UMAPs on the epithelial-only UMAP (n=575,366).
    Layout: 4 rows x 4 cols
      Row 1: SecA factor usage + SecA genes 1-3
      Row 2: SecA genes 4-7
      Row 3: SecB factor usage + SecB genes 1-3
      Row 4: SecB genes 4-7

INPUTS
    output_root/fig_secretory_polarization/data/meta.parquet  (epithelial UMAP coords)
    obj("atlas_epithelial")  (hgsc_atlas_epithelial.h5ad; log-normalised expression, backed)
    output_root/11d_epithelial_nmf/11d_nmf_usage.csv  (Factor_3 = SecA, Factor_2 = SecB)
    shared/signatures.yml  (canonical 7-gene noBCAM SecA / SecB lists)

OUTPUTS
    output_root/figures/supplementary/SF6_atlas_secAB_expression_umaps.{svg,png}

MANUSCRIPT PANEL(S)
    SF6A-B.

RUNTIME TIER
    moderate (backed expression extraction for 14 markers + 16 rasterized scatters).
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
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (central config)
# ============================================================================

META_PQ = path("output_root", "fig_secretory_polarization", "data", "meta.parquet")
H5AD = obj("atlas_epithelial")
USAGE_CSV = path("output_root", "11d_epithelial_nmf", "11d_nmf_usage.csv")
OUT_SVG = path("output_root", "figures", "supplementary", "SF6_atlas_secAB_expression_umaps.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF6_atlas_secAB_expression_umaps.png")

assert os.path.exists(META_PQ), f"Missing: {META_PQ}"
assert os.path.exists(H5AD), f"Missing: {H5AD}"
assert os.path.exists(USAGE_CSV), f"Missing: {USAGE_CSV}"

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

print("\nLoading epithelial metadata...", flush=True)
meta = pd.read_parquet(META_PQ, columns=["UMAP1", "UMAP2"])
print(f"  {len(meta):,} cells")

print("Loading h5ad (backed)...", flush=True)
adata = ad.read_h5ad(H5AD, backed="r")
print(f"  h5ad shape: {adata.shape}")

missing = [g for g in ALL_GENES if g not in adata.var_names]
assert not missing, f"Genes not found in epithelial h5ad: {missing}"

gene_to_col = {g: int(adata.var_names.get_loc(g)) for g in ALL_GENES}
h5_index = pd.Index(adata.obs_names)
meta_idx_in_h5 = h5_index.get_indexer(meta.index)
assert (meta_idx_in_h5 == -1).sum() == 0, "Unmatched cells"

print(f"Extracting expression for {len(ALL_GENES)} markers...", flush=True)
col_indices = [gene_to_col[g] for g in ALL_GENES]
expr = adata.X[meta_idx_in_h5, :][:, col_indices]
if hasattr(expr, "toarray"):
    expr = expr.toarray()
expr = np.asarray(expr, dtype=np.float32)
print(f"  Expression matrix: {expr.shape}")

adata.file.close()
del adata
gc.collect()

print("Loading NMF usage...", flush=True)
usage = pd.read_csv(USAGE_CSV, index_col=0)
usage = usage.loc[meta.index]
secA_usage = usage["Factor_3"].to_numpy(dtype=np.float32)
secB_usage = usage["Factor_2"].to_numpy(dtype=np.float32)

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

panels_secA = [("SecA NMF factor", secA_usage, False)]
for g in SECA_GENES:
    panels_secA.append((g, expr[:, ALL_GENES.index(g)], True))

panels_secB = [("SecB NMF factor", secB_usage, False)]
for g in SECB_GENES:
    panels_secB.append((g, expr[:, ALL_GENES.index(g)], True))

grid = [
    panels_secA[:4],
    panels_secA[4:],
    panels_secB[:4],
    panels_secB[4:],
]
SECTION_LABELS = {0: "SecA markers", 2: "SecB markers"}

# ============================================================================
# PLOT — 4 rows × 4 cols on portrait letter page
# ============================================================================

print("Plotting...", flush=True)

NROWS, NCOLS = 4, 4
PAGE_W, PAGE_H = 7.5, 10.0

fig, axes = plt.subplots(NROWS, NCOLS, figsize=(PAGE_W, PAGE_H),
                         gridspec_kw={"wspace": 0.08, "hspace": 0.12})

for row_idx, row_panels in enumerate(grid):
    for col_idx in range(NCOLS):
        ax = axes[row_idx, col_idx]
        if col_idx >= len(row_panels):
            ax.set_visible(False)
            continue

        label, vals, is_gene = row_panels[col_idx]
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
        cb.set_ticklabels([f"{0.0:.1f}", f"{vmax:.1f}"])

    if row_idx in SECTION_LABELS:
        axes[row_idx, 0].text(-0.08, 0.5, SECTION_LABELS[row_idx],
                              transform=axes[row_idx, 0].transAxes, fontsize=FA,
                              fontweight="bold", ha="right", va="center", rotation=90)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.06)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.06)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
