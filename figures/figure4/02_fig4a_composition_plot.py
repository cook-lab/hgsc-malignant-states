#!/usr/bin/env python3
"""
Figure 4a — Per-patient cell-type composition with hierarchical-clustering dendrogram
=====================================================================================
PURPOSE
    Per-patient stacked-bar cell-type composition (SingleR labels) with a ward.D2
    hclust dendrogram on top and a tissue-type strip at the bottom (FTE cores,
    whole tissues, TMA patients); tissue and cell-type keys on the right.

INPUTS
    - COMP_CSV (env; default output_root/singler_recompute/fig4a_composition_by_patient_SINGLER.csv): per-patient composition table (patient_id, tissue_type,
      cell_label, n) from figure4/01_fig4a_composition_table.py (TMA arm = the final
      TMA core set)
    - PANEL_A_ONLY (env, default 1): the Fig 4a layout (bars + dendrogram + keys)

OUTPUTS
    - figures_dir/figure4/xenium_composition_by_tissue_dendro_panelA.{png,svg}
      (SVG: the bars and tissue strip are one 600-dpi image; dendrogram, axes and keys stay vector)

MANUSCRIPT PANEL(S): Fig 4a.

RUNTIME TIER: fast.

NOTE: an epithelial label "Transitioning epithelium" in the input is read as "Intermediate epithelium".
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import gridspec
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram as scipy_dendro
from scipy.spatial.distance import pdist

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, CFG  # noqa: E402

# ---------- Paths ----------
STEM = "xenium_composition_by_tissue_dendro"

# ---------- Palette & order ----------
CELL_ORDER = [
    "SecA epithelium", "Intermediate epithelium", "SecB epithelium",
    "Ciliated epithelium", "Mesothelial", "Fibroblast", "Smooth muscle",
    "Pericyte", "Endothelial", "T cell", "NK cell", "B cell", "Plasma cell",
    "Macrophage", "Conventional dendritic cell", "Neutrophil", "Mast cell",
    "Plasmacytoid dendritic cell",
]
REF_PALETTE = {
    "SecA epithelium": "#E6A141", "Intermediate epithelium": "#C08E48",
    "SecB epithelium": "#9A7D55", "Ciliated epithelium": "#E07850",
    "Mesothelial": "#D4A574", "Fibroblast": "#DDD5CA", "Smooth muscle": "#D14E6C",
    "Pericyte": "#B87A7A", "Endothelial": "#7D4E4E", "T cell": "#87CEFA",
    "NK cell": "#56AFC4", "B cell": "#5665B6", "Plasma cell": "#8A5DAF",
    "Macrophage": "#8FBC8F", "Conventional dendritic cell": "#2E8B57",
    "Neutrophil": "#6B8E23", "Mast cell": "#8B9B6B",
    "Plasmacytoid dendritic cell": "#A0A0A0",
}

TISSUE_ORDER = ["Fallopian tube", "HGSC whole tissue", "HGSC TMA"]
TISSUE_COLORS = {"Fallopian tube": "#E07850", "HGSC TMA": "#D4A574",
                 "HGSC whole tissue": "#8A5DAF"}

# ---------- Style ----------
FA, FK, FN, FL = 7, 6, 6, 6
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FK,
    "axes.titlesize":    FA,
    "axes.labelsize":    FA,
    "xtick.labelsize":   FK,
    "ytick.labelsize":   FK,
    "legend.fontsize":   FL,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
    "axes.linewidth":    0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})


def main():
    t0 = time.time()
    print(f"  {STEM}\n")

    # COMP_CSV: the per-patient composition table (patient_id, tissue_type, cell_label, n) from
    # figure4/01_fig4a_composition_table.py, whose TMA arm is the final TMA core set (default: where that
    # script writes it). SPLIT_COMP_CSV is an older name of the same variable.
    _comp = os.environ.get("COMP_CSV", os.environ.get("SPLIT_COMP_CSV", str(path(
        "output_root", "singler_recompute", "fig4a_composition_by_patient_SINGLER.csv"))))
    if _comp:
        full = pd.read_csv(_comp, dtype={"patient_id": str})
        full["cell_label"] = full["cell_label"].replace({"Transitioning epithelium": "Intermediate epithelium"})
        full = full[["patient_id", "tissue_type", "cell_label", "n"]]
        print(f"  composition table: {_comp}")

    pt_tis = full.groupby(["patient_id", "tissue_type"], as_index=False)["n"].sum()
    grid = pt_tis[["patient_id", "tissue_type"]].merge(
        pd.DataFrame({"cell_label": CELL_ORDER}), how="cross")
    full = grid.merge(full, on=["patient_id", "tissue_type", "cell_label"],
                      how="left").fillna({"n": 0})
    full["n"] = full["n"].astype(int)
    full["total"] = full.groupby(["patient_id", "tissue_type"])["n"].transform("sum")
    full = full[full["total"] > 0].copy()
    full["proportion"] = full["n"] / full["total"]
    full["bar"] = full["patient_id"].astype(str) + "|" + full["tissue_type"]
    n_bars = full["bar"].nunique()
    print(f"  bars: {n_bars}")

    mat = (full.pivot_table(index="bar", columns="cell_label",
                            values="proportion", fill_value=0)
               .reindex(columns=CELL_ORDER, fill_value=0))
    from scipy.cluster.hierarchy import optimal_leaf_ordering, to_tree
    d = pdist(mat.values, metric="euclidean")
    Z = linkage(d, method="ward")
    Z = optimal_leaf_ordering(Z, d)   # group similar leaves (un-split the FTE)

    bar_to_tissue = dict(zip(full["bar"], full["tissue_type"]))
    # Rotate every node so the FTE-rich subtree is its RIGHT child -> the FTE
    # clade is pushed to the right edge. Only branch orientation changes; the
    # clustering (merges/heights) is untouched. Dendrogram is drawn from this
    # rotated traversal so bars and tree stay aligned.
    ft_leaf = np.array([bar_to_tissue[b] == "Fallopian tube" for b in mat.index], dtype=int)
    _tree = to_tree(Z)
    _ftc = {}

    def _cnt(nd):
        v = int(ft_leaf[nd.id]) if nd.is_leaf() else _cnt(nd.left) + _cnt(nd.right)
        _ftc[nd.id] = v
        return v

    _cnt(_tree)
    dend_lines = []
    leaf_x = {}
    _ctr = [0]

    def _place(nd):
        if nd.is_leaf():
            x = _ctr[0]; leaf_x[nd.id] = x; _ctr[0] += 1
            return x, 0.0
        l, r = nd.left, nd.right
        if _ftc[l.id] > _ftc[r.id]:
            l, r = r, l                      # FTE-rich subtree -> right
        xl, hl = _place(l); xr, hr = _place(r); h = nd.dist
        dend_lines.extend([("v", xl, hl, h), ("v", xr, hr, h), ("h", xl, xr, h)])
        return (xl + xr) / 2.0, h

    _place(_tree)
    ordered_bars = [mat.index[i] for i in sorted(leaf_x, key=leaf_x.get)]
    tissue_counts = {t: sum(1 for b in ordered_bars if bar_to_tissue[b] == t)
                     for t in TISSUE_ORDER}
    print(f"  tissue breakdown: {tissue_counts}")

    props = mat.reindex(index=ordered_bars)

    print("\n  Rendering ...")
    # PANEL_A_ONLY=1 (the default): dendrogram + composition bars, the Fig 4a layout.
    PANEL_A = os.environ.get("PANEL_A_ONLY", "1") != "0"
    if PANEL_A:
        # Short panel for the Fig 4 layout (~0.99 in tall). Font sizes for the y-axis
        # label/title/ticks are intentionally NOT reduced (already small); legend
        # retained so the repo reproduces the full panel.
        fig = plt.figure(figsize=(7.6, 0.99))
        gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[0.6, 3.0],
                               hspace=0.02, left=0.06, right=0.80, top=0.96, bottom=0.05)
        ax_dend = fig.add_subplot(gs[0]); ax_bar = fig.add_subplot(gs[1]); ax_don = None

    # Draw the (rotated) dendrogram from the collected segments; leaf x-positions
    # are 0..n-1, matching the composition bars below.
    _hmax = max(h for _, _, _, h in dend_lines)
    for kind, a, b, h in dend_lines:
        if kind == "v":
            ax_dend.plot([a, a], [b, h], color="black", linewidth=0.4)
        else:
            ax_dend.plot([a, b], [h, h], color="black", linewidth=0.4)
    ax_dend.set_xlim(-0.5, n_bars - 0.5)
    ax_dend.set_ylim(0, _hmax * 1.02)
    ax_dend.set_xticks([]); ax_dend.set_yticks([])
    for spine in ax_dend.spines.values():
        spine.set_visible(False)

    # bars and tissue strip are rasterized in the SVG (2026-10-21, author): ~7,600 thin
    # rectangles become one image; consecutive rasterized artists are merged by matplotlib
    x = np.arange(n_bars)
    bottoms = np.zeros(n_bars)
    for ct in CELL_ORDER:
        y = props[ct].values
        ax_bar.bar(x, y, bottom=bottoms, width=1.0, color=REF_PALETTE[ct],
                   edgecolor="none", linewidth=0, label=ct, rasterized=True)
        bottoms += y

    strip_y, strip_h = -0.035, 0.03
    for i, b in enumerate(ordered_bars):
        t = bar_to_tissue[b]
        ax_bar.add_patch(plt.Rectangle((i - 0.5, strip_y), 1.0, strip_h,
                                       facecolor=TISSUE_COLORS[t], edgecolor="none", rasterized=True))
    ax_bar.set_xlim(-0.5, n_bars - 0.5)
    ax_bar.set_ylim(strip_y - 0.01, 1.002)
    ax_bar.set_xticks([])
    ax_bar.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax_bar.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax_bar.set_ylabel("Cell type proportion")
    for s in ("top", "right", "bottom"):
        ax_bar.spines[s].set_visible(False)
    ax_bar.tick_params(axis="y", length=2, pad=1.5)
    ax_dend.set_xlim(ax_bar.get_xlim())

    tissue_handles = [Patch(facecolor=TISSUE_COLORS[t], edgecolor="none",
                            label=f"{t} (n={tissue_counts[t]} patients)") for t in TISSUE_ORDER]
    leg_tissue = ax_bar.legend(handles=tissue_handles, title="Tissue",
                               bbox_to_anchor=(1.01, 1.0), loc="upper left",
                               frameon=False, fontsize=FL, title_fontsize=FL,
                               handlelength=1.2, handleheight=0.8, borderpad=0.2,
                               labelspacing=0.3)
    ax_bar.add_artist(leg_tissue)
    # The cell-type key sits beside the tissue key, not under it.
    ax_bar.legend(handles=[Patch(facecolor=REF_PALETTE[c], edgecolor="none", label=c)
                           for c in CELL_ORDER], title="Cell type",
                  bbox_to_anchor=(1.29, 1.0) if PANEL_A else (1.01, 0.72), loc="upper left", frameon=False,
                  fontsize=FL, title_fontsize=FL, handlelength=1.2, handleheight=0.8,
                  borderpad=0.2, labelspacing=0.25)

    _suffix = "_panelA" if PANEL_A else ""
    out_png = path("figures_dir", "figure4", f"{STEM}{_suffix}.png")
    out_svg = path("figures_dir", "figure4", f"{STEM}{_suffix}.svg")
    fig.savefig(out_png, dpi=450, bbox_inches="tight")
    fig.savefig(out_svg, dpi=600, bbox_inches="tight")   # dpi = the rasterized bars only
    plt.close(fig)
    print(f"\n  -> {out_png}\n  -> {out_svg}\n  Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
