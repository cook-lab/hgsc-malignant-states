#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202608_kpca_secb_mouse/output/41_fig3_umaps/41b_export_1p8in.py, the script that ran.
# Kept: 41_mouse_umap_compartment_1p8in.{svg,png}, which is Fig 3i. Removed: the epithelial SecB, polarization and
# Krt7 panels (figures/figure3/04_... and 05_... draw Fig 3j, k) and an unused palette. The untrimmed script, which
# the docstring below describes, is in the lab project folder.

"""
Figure 3 — KPCA.B mouse UMAP panels re-exported as EQUAL 1.8x1.8" SQUARES
=========================================================================
Same data/columns as 41_build_mouse_umaps.py; layout only. Intended for a
4-across letter-page row (4th square = KRT7 IHC added physically by the user).

Changes vs the 3.0"/2.4" build:
  * All three panels forced to a true 1.8 x 1.8 inch footprint (equal squares).
  * Panel i (compartments): NO in-plot legend (won't fit at 1.8") — clean plot
    for an EXTERNAL legend. Compartment colours made all-distinct using atlas
    cell-type hexes (two too-close pairs reassigned; see COMP_PALETTE notes).
  * Panels j/k: same brown SecB colourmap for SecB<->KRT7 comparability; small
    inset colourbar retained.

Outputs (svg + png, svg.fonttype:none):
  41_mouse_umap_compartment_1p8in.{svg,png}
  41_mouse_epi_umap_secB_1p8in.{svg,png}
  41_mouse_epi_umap_KRT7_1p8in.{svg,png}
"""
from __future__ import annotations
import os
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SEED = 0
H5AD = "/Volumes/CookLab/Sarah/202608_kpca_secb_mouse/output/33_handoff/kpca_final_annotated.h5ad"
OUTDIR = "/Volumes/CookLab/Sarah/202608_kpca_secb_mouse/output/41_fig3_umaps"
os.makedirs(OUTDIR, exist_ok=True)
SIDE = 1.8  # inches, equal squares

FA, FK, FN = 8, 6, 6
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       FA, "axes.labelsize": FA, "axes.titlesize": 0,
    "xtick.labelsize": FK, "ytick.labelsize": FK, "legend.fontsize": FN,
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "figure.dpi": 450, "savefig.dpi": 450,
})

# ---------- Panel i is coloured by the CORRECTED `celltype` column ----------
# (`compartment` is the old coarse call that mis-split the tumour manifold into
#  Fibroblast/Other via a Col1a1 artifact. `celltype` Tumour = the true
#  epithelial/malignant manifold = the whole bottom-left blob = the 12,429-cell
#  set used by panels j/k.) Tumour is relabelled "Epithelial" for the legend key.
# All-distinct atlas cell-type hexes (min pairwise CIE76 dE = 26).
CELLTYPE_RELABEL = {"Tumour": "Epithelial"}
CELLTYPE_PALETTE = {
    "Epithelial":  "#E6A141",   # atlas Epithelial (orange)  <- celltype 'Tumour'
    "Macrophage":  "#8FBC8F",   # atlas Macrophage (green)
    "T/NK":        "#87CEFA",   # atlas T/NK cell (light blue)
    "Endothelial": "#7D4E4E",   # atlas Endothelial (dark brown)
    "Neutrophil":  "#6B8E23",   # atlas Neutrophil (olive)
    "Fibroblast":  "#DDD5CA",   # atlas Fibroblast (beige)
    "Monocyte":    "#D14E6C",   # atlas Smooth muscle (red) — distinct from greens
    "Dendritic":   "#2E8B57",   # atlas DC (dark sea green)
    "B":           "#5665B6",   # atlas B cell (indigo)
}
CELLTYPE_ORDER = ["Epithelial", "Macrophage", "T/NK", "Endothelial",
                  "Neutrophil", "Fibroblast", "Monocyte", "Dendritic", "B"]

a = ad.read_h5ad(H5AD)
print(f"[load] total cells: {a.n_obs:,}")


def stubs(ax, x, y, xr, yr, px, py):
    ax.set_xlim(x.min() - px, x.max() + px)
    ax.set_ylim(y.min() - py, y.max() + py)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_bounds(x.min() - px, x.min() - px + 0.18 * xr)
    ax.spines["left"].set_bounds(y.min() - py, y.min() - py + 0.18 * yr)
    ax.spines["bottom"].set_linewidth(0.4); ax.spines["left"].set_linewidth(0.4)


# ============ PANEL i — full UMAP by CELLTYPE, external legend ============
# Canvas widened to 2.7 x 1.8"; the UMAP axes stay a true 1.8 x 1.8" square on
# the left, the cell-type key sits in the 0.9" strip to its right.
xy = np.asarray(a.obsm["X_umap"]); x, y = xy[:, 0], xy[:, 1]
ct_raw = a.obs["celltype"].astype(str).values
ct = np.array([CELLTYPE_RELABEL.get(v, v) for v in ct_raw])   # Tumour -> Epithelial
o = np.random.default_rng(SEED).permutation(a.n_obs)
xr, yr = x.max() - x.min(), y.max() - y.min(); px, py = 0.04 * xr, 0.04 * yr

CANVAS_W, CANVAS_H = 2.7, 1.8   # UMAP square preserved at 1.8 x 1.8"
fig = plt.figure(figsize=(CANVAS_W, CANVAS_H))
ax = fig.add_axes([0.0, 0.0, SIDE / CANVAS_W, 1.0])   # 1.8" wide x 1.8" tall square
fig.subplots_adjust()
for cat in CELLTYPE_ORDER:
    m = ct[o] == cat
    if not m.any():
        continue
    ax.scatter(x[o][m], y[o][m], c=CELLTYPE_PALETTE[cat], s=1.0, marker="o",
               linewidths=0, rasterized=True, alpha=0.9)
stubs(ax, x, y, xr, yr, px, py)
handles = [Line2D([0], [0], marker="o", linestyle="", markerfacecolor=CELLTYPE_PALETTE[c],
                  markeredgecolor="none", markersize=4, label=c)
           for c in CELLTYPE_ORDER if (ct == c).any()]
leg = ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.00, 0.5),
                frameon=False, fontsize=7, handletextpad=0.4, labelspacing=0.42,
                borderaxespad=0.0)
for h in leg.legend_handles:
    h.set_markersize(4)
pi = os.path.join(OUTDIR, "41_mouse_umap_compartment_1p8in")
fig.savefig(pi + ".svg", format="svg", bbox_inches="tight", pad_inches=0.02)
fig.savefig(pi + ".png", format="png", bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print(f"[save] {pi}.svg/.png  canvas={CANVAS_W}x{CANVAS_H}in, UMAP square={SIDE}in")

print("Done.")
