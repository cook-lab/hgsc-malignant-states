#!/usr/bin/env python3
"""
SF9 — Epithelial composition by clinical / genomic metadata
===========================================================

Purpose
    Six panels of mean epitype composition (SecA / Intermediate / SecB /
    Ciliated) stratified by: treatment status, anatomic site, metastatic site,
    TP53 status, HRD status, BRCA status. Portrait letter page.

INPUTS
    output_root/fig_secretory_polarization/data/meta.parquet
        (575,366 epithelial cells; schema_nmf + treatment_status, anatomic_site,
         metastatic_site)
    obj("atlas_epithelial")  (hgsc_atlas_epithelial.h5ad; TP53_status, HRD_status,
        BRCA_status, joined via cell barcode index)

OUTPUTS
    output_root/figures/supplementary/SF9_atlas_epitype_by_metadata.{svg,png}

MANUSCRIPT PANEL(S)
    SF9A-E.

RUNTIME TIER
    moderate (backed read of epithelial obs for genomic columns).
"""

import os
import gc
import sys

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

META_PQ = path("output_root", "fig_secretory_polarization", "data", "meta.parquet")
EPI_H5AD = obj("atlas_epithelial")
OUT_SVG = path("output_root", "figures", "supplementary", "SF9_atlas_epitype_by_metadata.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF9_atlas_epitype_by_metadata.png")

assert os.path.exists(META_PQ), f"Missing: {META_PQ}"

# ============================================================================
# STYLE
# ============================================================================

FA, FK, FN = 6, 5.5, 5

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       FA,
    "axes.labelsize":  FA,
    "axes.titlesize":  0,
    "xtick.labelsize": FK,
    "ytick.labelsize": FK,
    "legend.fontsize": FN,
    "svg.fonttype":    "none",
    "pdf.fonttype":    42,
    "figure.dpi":      450,
    "savefig.dpi":     450,
})

# Epitype labels: "Intermediate" (was "Transitioning").
EPI_PALETTE = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55", "Ciliated": "#E07850"}
EPI_ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]

# ============================================================================
# LOAD & MERGE DATA
# ============================================================================

print("Loading epithelial metadata...", flush=True)
meta = pd.read_parquet(META_PQ, columns=[
    "schema_nmf", "treatment_status", "anatomic_site", "metastatic_site",
])
print(f"  {len(meta):,} cells")

# Standardize label "Transitioning" -> "Intermediate" in the schema column.
meta["schema_nmf"] = meta["schema_nmf"].replace({"Transitioning": "Intermediate"})

print("Loading epithelial h5ad for genomic annotations...", flush=True)
adata = ad.read_h5ad(EPI_H5AD, backed="r")
genomic = adata.obs[["TP53_status", "HRD_status", "BRCA_status"]].copy()
adata.file.close()
del adata
gc.collect()

common = meta.index.intersection(genomic.index)
print(f"  Matched {len(common):,} / {len(meta):,} cells for genomic annotations")
meta = meta.loc[common].copy()
for col in ["TP53_status", "HRD_status", "BRCA_status"]:
    meta[col] = genomic.loc[common, col].values
    meta[col] = meta[col].fillna("NA")
    meta.loc[meta[col].isin(["False", "nan"]), col] = "NA"
del genomic
gc.collect()

# ============================================================================
# PANEL DEFINITIONS
# ============================================================================

PANELS = [
    {
        "col": "treatment_status", "title": "Treatment status",
        "order": ["pre-treatment", "post-chemotherapy", "post-niraparib",
                  "post-chemotherapy_niraparib", "post-chemotherapy_olaparib",
                  "post-chemotherapy_pembro"],
        "labels": {
            "pre-treatment": "Pre-treatment", "post-chemotherapy": "Post-chemo",
            "post-niraparib": "Post-niraparib",
            "post-chemotherapy_niraparib": "Post-chemo\n+ niraparib",
            "post-chemotherapy_olaparib": "Post-chemo\n+ olaparib",
            "post-chemotherapy_pembro": "Post-chemo\n+ pembro",
        },
    },
    {
        "col": "anatomic_site", "title": "Anatomic site",
        "order": ["adnexa", "omentum", "ascites", "peritoneum", "bowel",
                  "upper_quadrant", "pelvic_cavity", "lymph_node"],
        "labels": {
            "adnexa": "Adnexa", "omentum": "Omentum", "ascites": "Ascites",
            "peritoneum": "Peritoneum", "bowel": "Bowel",
            "upper_quadrant": "Upper\nquadrant", "pelvic_cavity": "Pelvic\ncavity",
            "lymph_node": "Lymph\nnode",
        },
    },
    {
        "col": "metastatic_site", "title": "Metastatic site",
        "order": ["primary", "metastasis", "ascites"],
        "labels": {"primary": "Primary", "metastasis": "Metastasis", "ascites": "Ascites"},
    },
    {
        "col": "TP53_status", "title": "TP53 status",
        "order": ["mutated", "wildtype"],
        "labels": {"mutated": "Mutated", "wildtype": "Wildtype"},
    },
    {
        "col": "HRD_status", "title": "HRD status",
        "order": ["HRD", "HRP"],
        "labels": {"HRD": "HRD", "HRP": "HR Proficient"},
    },
    {
        "col": "BRCA_status", "title": "BRCA status",
        "order": ["mutated", "wildtype"],
        "labels": {"mutated": "Mutated", "wildtype": "Wildtype"},
    },
]

# ============================================================================
# COMPUTE COMPOSITIONS
# ============================================================================

print("\nComputing compositions...", flush=True)
for panel in PANELS:
    col = panel["col"]
    order = [g for g in panel["order"] if g in meta[col].values]
    panel["order"] = order
    records = []
    for grp in order:
        mask = meta[col] == grp
        n = mask.sum()
        if n == 0:
            continue
        row = {"group": grp, "n": n}
        for ep in EPI_ORDER:
            row[ep] = (meta.loc[mask, "schema_nmf"] == ep).sum() / n * 100
        records.append(row)
    panel["data"] = pd.DataFrame(records)
    print(f"  {panel['title']}: {len(panel['data'])} groups")

# ============================================================================
# PLOT — Row 1: 2 wide panels; Row 2: 4 narrow panels
# ============================================================================

print("\nPlotting...", flush=True)

PAGE_W, PAGE_H = 7.5, 6.0
fig = plt.figure(figsize=(PAGE_W, PAGE_H))
gs = GridSpec(2, 4, figure=fig, hspace=0.60, wspace=0.12,
              top=0.94, bottom=0.09, left=0.07, right=0.97, height_ratios=[1, 1])

panel_specs = [
    gs[0, 0:2], gs[0, 2:4], gs[1, 0], gs[1, 1], gs[1, 2], gs[1, 3],
]

for idx, panel in enumerate(PANELS):
    ax = fig.add_subplot(panel_specs[idx])
    df = panel["data"]
    labels_map = panel["labels"]
    n_bars = len(df)
    x_pos = np.arange(n_bars)
    BAR_W = 0.65

    bottom = np.zeros(n_bars)
    for ep in EPI_ORDER:
        vals = df[ep].values
        ax.bar(x_pos, vals, width=BAR_W, bottom=bottom, color=EPI_PALETTE[ep],
               edgecolor="white", linewidth=0.3)
        bottom += vals

    bottom = np.zeros(n_bars)
    for ep in EPI_ORDER:
        vals = df[ep].values
        for i, v in enumerate(vals):
            if v >= 8:
                ax.text(x_pos[i], bottom[i] + v / 2, f"{v:.0f}%", fontsize=FN - 1,
                        ha="center", va="center", color="white", fontweight="bold", alpha=0.9)
        bottom += vals

    x_labels = []
    for grp in df["group"]:
        lbl = labels_map.get(grp, grp)
        n = int(df.loc[df["group"] == grp, "n"].values[0])
        x_labels.append(f"{lbl}\n(n={n:,})")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, fontsize=FN - 0.5, linespacing=0.9)
    ax.tick_params(axis="x", length=0, pad=2)

    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    if idx in (0, 2):
        ax.set_ylabel("Epithelial composition (%)", fontsize=FA)
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])

    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.spines["left"].set_linewidth(0.4)
    ax.tick_params(axis="y", length=2, width=0.4)
    ax.text(0.5, 1.06, panel["title"], transform=ax.transAxes, fontsize=FA + 0.5,
            fontweight="bold", ha="center", va="bottom")

handles = [Patch(facecolor=EPI_PALETTE[ep], edgecolor="white", linewidth=0.3, label=ep)
           for ep in EPI_ORDER]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=FN + 0.5,
           handlelength=1.2, handleheight=0.9, handletextpad=0.4, columnspacing=1.5,
           bbox_to_anchor=(0.5, 0.01))

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.06)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.06)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
