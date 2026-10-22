#!/usr/bin/env python3
"""
Figure 1j — Orthogonal NMF program UMAPs (G2M / AP-1 / IFN / EMT)
=================================================================
PURPOSE
    Four single-panel epithelial UMAPs coloured by continuous NMF factor usage
    for the biological programs that vary ORTHOGONALLY to the SecA<->SecB
    polarization axis (Fig 1h,i). Same embedding + styling as the SecA/SecB
    factor UMAPs of Fig 1h,i; a shared magenta ramp distinguishes these "other
    programs" from the gold(SecA)/brown(SecB) poles.

    Factor -> program identity is taken from the atlas NMF characterization
    (11d_epithelial_nmf). The four sub-panels of Fig 1j, in order:
        G2M / cell cycle       = Factor_5  (TUBA1B, HMGB2, CENPF, STMN1, UBE2S)
        AP-1 / immediate-early = Factor_6  (FOS, JUN, ATF3, KLF6, JUNB)
        IFN / interferon       = Factor_4  (IFI27, ISG15, IFI6, LY6E, WFDC2)
        EMT-like               = Factor_9  (VIM, FN1, LGALS1, S100A6, S100A10)

    tables/05_tableS5_epithelial_nmf_programs.py reads PANELS from this file by
    path (Table S5 factor labels); keep the file name and PANELS as they are.

INPUTS
    - figures/_prep/atlas_final_tables.py :: load_epithelial()  (575,366
      epithelial cells: epithelial-embedding UMAP1/UMAP2 + NMF usage
      Factor_5/6/4/9 from 11d_epithelial_nmf/11d_nmf_usage.csv; cache
      output_root/_prep_caches/atlas_final_epithelial.parquet).
OUTPUTS
    - figures_dir/figure1/atlas_program_{g2m,ap1,ifn,emt}_umap.{svg,png}
MANUSCRIPT PANEL(S): Fig 1j (four sub-panels: G2M, AP-1, IFN, EMT).
RUNTIME TIER: fast (reads one cached parquet; renders ~575k points x4).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

# ---------- Panels: (out-stem, factor, label, top genes) --------------------
PANELS = [
    ("atlas_program_g2m_umap", "Factor_5", "G2M",  "(TUBA1B, HMGB2, CENPF)"),
    ("atlas_program_ap1_umap", "Factor_6", "AP-1", "(FOS, JUN, ATF3)"),
    ("atlas_program_ifn_umap", "Factor_4", "IFN",  "(IFI27, ISG15, IFI6)"),
    ("atlas_program_emt_umap", "Factor_9", "EMT",  "(VIM, FN1, LGALS1)"),
]

# Neutral beige -> magenta ramp (distinct from SecA gold / SecB brown).
CMAP = LinearSegmentedColormap.from_list(
    "beige_program_magenta",
    [(0.00, "#DDD5CA"), (0.40, "#DDD5CA"), (0.55, "#D08FB4"),
     (0.70, "#B0508A"), (0.85, "#8E2C6B"), (1.00, "#5E1747")], N=256,
)
LABEL_COLOR = "#8E2C6B"

# ---------- Style (identical to the Fig 1h,i UMAPs) -------------------------
FA, FK, FN, FG = 8, 6, 8, 8
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "axes.titlesize": 0,
    "xtick.labelsize": FK, "ytick.labelsize": FK, "legend.fontsize": FN,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 450, "savefig.dpi": 450,
})

# ---------- Load (once) -----------------------------------------------------
meta = load_epithelial(columns=["UMAP1", "UMAP2"] + [p[1] for p in PANELS])
assert len(meta) == 575_366, len(meta)
usage = meta

x = meta["UMAP1"].to_numpy(); y = meta["UMAP2"].to_numpy()
xr, yr = x.max() - x.min(), y.max() - y.min()
pad_x, pad_y = 0.04 * xr, 0.04 * yr
SHUF = np.random.default_rng(seed=SEED).permutation(len(meta))
x_shuf, y_shuf = x[SHUF], y[SHUF]

# ---------- Render each panel ----------------------------------------------
for stem, factor, label, genes in PANELS:
    vals = usage[factor].to_numpy()
    vo = vals[SHUF]
    vmax = float(np.quantile(vals, 0.99)) or 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)

    fig, ax = plt.subplots(1, 1, figsize=(3.0, 2.30))
    ax.scatter(x_shuf, y_shuf, c=vo, cmap=CMAP, norm=norm, s=0.08, marker="o",
               linewidths=0, rasterized=True, alpha=0.95)
    ax.set_xlim(x.min() - pad_x, x.max() + pad_x)
    ax.set_ylim(y.min() - pad_y, y.max() + pad_y)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_bounds(x.min() - pad_x, x.min() - pad_x + 0.18 * xr)
    ax.spines["left"].set_bounds(y.min() - pad_y, y.min() - pad_y + 0.18 * yr)
    ax.spines["bottom"].set_linewidth(0.4); ax.spines["left"].set_linewidth(0.4)

    ax.text(0.03, 0.97, label, transform=ax.transAxes, fontsize=FG,
            fontweight="bold", ha="left", va="top", color=LABEL_COLOR)
    ax.text(0.03, 0.88, genes, transform=ax.transAxes, fontsize=FG,
            fontstyle="italic", ha="left", va="top", color="#444444")

    cax = ax.inset_axes([0.88, 0.05, 0.04, 0.3])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax, orientation="vertical")
    cb.outline.set_linewidth(0.3)
    cb.ax.tick_params(length=1.5, width=0.3, pad=1, labelsize=FN - 1)
    cb.set_ticks([0.0, vmax]); cb.set_ticklabels([f"{0.0:.2f}", f"{vmax:.2f}"])

    out_svg = path("figures_dir", "figure1", f"{stem}.svg")
    out_png = path("figures_dir", "figure1", f"{stem}.png")
    fig.savefig(out_svg, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(out_png, format="png", dpi=450, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print(f"  {label:5s} ({factor}) vmax={vmax:.3f} -> {stem}.svg")

print("Done: Fig 1J-M program UMAPs.")
