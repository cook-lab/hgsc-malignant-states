#!/usr/bin/env python3
"""
Figure S9 — Epithelial composition by clinical / genomic metadata / published annotations
=========================================================================================

Purpose
    Seven panels of epitype composition (SecA / Intermediate / SecB / Ciliated; pooled
    % of epithelial cells in each group, all 575,366 epithelial cells, no subsampling)
    stratified by:
      a treatment status, b anatomic site, c TP53 status, d HRD status, e BRCA status,
      f Vázquez-García et al. 2022 (SPECTRUM) epithelial cancer-cell cluster labels,
      g Geistlinger et al. 2020 per-cell consensusOV subtype.
    One standalone file per panel (no panel letters).

INPUTS
    figures._prep.atlas_final_tables.load_epithelial()
        (575,366 epithelial cells of the final atlas: `state` = SecA / Intermediate / SecB /
         Ciliated — Factor_2 p50/p75 split of non-ciliated cells — plus final-atlas obs study,
         treatment_status, anatomic_site, TP53_status, HRD_status, BRCA_status).
        TP53/HRD/BRCA are plotted as held in the final atlas (no recode;
        "False"/"nan" count as NA and are not plotted).
    data_root/2026_final_atlas/output/spectrum_epithelial_crosstab/GSE180661_GEO_cells.tsv.gz
        (SPECTRUM cell table; cluster_label for cell_type == Ovarian.cancer.cell; join on
         cell_id == atlas obs index)
    data_root/2026_final_atlas/output/spectrum_epithelial_crosstab/geistlinger_cell_annotations.csv
        (Geistlinger SCE colData subtype; join_key == atlas obs index)

OUTPUTS
    output_root/figures/supplementary/SF9{a..g}_atlas_epitype_by_<metadata>.{svg,png}
    output_root/figures/supplementary/SF9_atlas_epitype_by_metadata_counts.csv (n and % per bar)

MANUSCRIPT PANEL(S)
    Figure S9a-g.

RUNTIME TIER
    fast (epithelial cache + two external annotation tables).
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

XTAB_DIR = ("data_root", "2026_final_atlas", "output", "spectrum_epithelial_crosstab")
SPECTRUM_TSV = path(*XTAB_DIR, "GSE180661_GEO_cells.tsv.gz")
GEIST_CSV = path(*XTAB_DIR, "geistlinger_cell_annotations.csv")
OUT_DIR = path("output_root", "figures", "supplementary")
OUT_COUNTS = os.path.join(OUT_DIR, "SF9_atlas_epitype_by_metadata_counts.csv")

for p in (SPECTRUM_TSV, GEIST_CSV):
    assert os.path.exists(p), f"Missing: {p}"

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

# Epitype labels and colours.
EPI_PALETTE = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55", "Ciliated": "#E07850"}
EPI_ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]

# ============================================================================
# LOAD & MERGE DATA
# ============================================================================

OBS_COLS = ["study", "treatment_status", "anatomic_site", "TP53_status", "HRD_status", "BRCA_status"]
print("Loading epithelial table (load_epithelial)...", flush=True)
meta = load_epithelial(columns=["state"] + OBS_COLS)
meta["schema_nmf"] = meta.pop("state").astype(str)
for col in OBS_COLS:
    meta[col] = meta[col].astype(str)
print(f"  {len(meta):,} epithelial cells")
for col in ["TP53_status", "HRD_status", "BRCA_status"]:
    print(f"  {col} categories (epithelial, as held in the final atlas): {meta[col].value_counts().to_dict()}")

for col in ["TP53_status", "HRD_status", "BRCA_status", "treatment_status", "anatomic_site"]:
    meta.loc[meta[col].isin(["False", "nan", "None", ""]), col] = "NA"

# --- F: Vázquez-García (SPECTRUM) epithelial cancer-cell cluster labels ---
print("Joining SPECTRUM cluster_label...", flush=True)
spec = pd.read_csv(SPECTRUM_TSV, sep="\t", usecols=["cell_id", "cell_type", "cluster_label"])
spec = spec[spec["cell_type"] == "Ovarian.cancer.cell"].set_index("cell_id")["cluster_label"]
meta["vazquez_garcia_label"] = "NA"
is_vg = meta["study"] == "vazquez_garcia_2022"
lab = spec.reindex(meta.index[is_vg])
meta.loc[is_vg, "vazquez_garcia_label"] = lab.fillna("Unassigned").values
print(f"  vazquez_garcia_2022 epithelial cells: {is_vg.sum():,}; "
      f"with a SPECTRUM cancer-cell label: {lab.notna().sum():,}")

# --- G: Geistlinger per-cell consensusOV subtype ---
print("Joining Geistlinger subtype...", flush=True)
gann = pd.read_csv(GEIST_CSV).set_index("join_key")
meta["geistlinger_label"] = "NA"
is_g = meta["study"] == "geistlinger_2020"
gl = gann["subtype"].reindex(meta.index[is_g])
meta.loc[is_g, "geistlinger_label"] = gl.fillna("Unassigned").values
print(f"  geistlinger_2020 epithelial cells: {is_g.sum():,}; with a subtype: {gl.notna().sum():,}")
print("  Geistlinger celltype of those cells:",
      gann["celltype"].reindex(meta.index[is_g]).value_counts(dropna=False).to_dict())

# ============================================================================
# PANEL DEFINITIONS
# ============================================================================

PANELS = [
    {
        "key": "a", "file": "treatment_status",
        "col": "treatment_status", "title": "Treatment status",
        "order": ["pre-treatment", "post-chemotherapy", "post-niraparib",
                  "post-chemotherapy_niraparib", "post-chemotherapy_olaparib",
                  "post-chemotherapy_pembro", "post-chemotherapy_niraparib_pembro"],   # last: Nath 2021 fix, 2026-10-20
        "labels": {
            "pre-treatment": "Pre-\ntreatment", "post-chemotherapy": "Post-\nchemo",
            "post-niraparib": "Post-\nniraparib",
            "post-chemotherapy_niraparib": "Post-chemo\n+ niraparib",
            "post-chemotherapy_olaparib": "Post-chemo\n+ olaparib",
            "post-chemotherapy_pembro": "Post-chemo\n+ pembro",
            "post-chemotherapy_niraparib_pembro": "Post-chemo\n+ niraparib\n+ pembro",
        },
    },
    {
        # all 12 atlas sites
        "key": "b", "file": "anatomic_site",
        "col": "anatomic_site", "title": "Anatomic site",
        "order": ["adnexa", "omentum", "ascites", "peritoneum", "bowel",
                  "upper_quadrant", "pelvic_cavity", "diaphragm", "bladder", "colon",
                  "lymph_node", "liver"],
        "labels": {
            "adnexa": "Adnexa", "omentum": "Omentum", "ascites": "Ascites",
            "peritoneum": "Peritoneum", "bowel": "Bowel",
            "upper_quadrant": "Upper\nquadrant", "pelvic_cavity": "Pelvic\ncavity",
            "diaphragm": "Diaphragm", "bladder": "Bladder", "colon": "Colon",
            "lymph_node": "Lymph\nnode", "liver": "Liver",
        },
    },
    # No metastatic-site panel (Figure S9 has none).
    {
        "key": "c", "file": "TP53_status",
        "col": "TP53_status", "title": "TP53 status",
        "order": ["mutated", "wildtype"],
        "labels": {"mutated": "Mutated", "wildtype": "Wildtype"},
    },
    {
        "key": "d", "file": "HRD_status",
        "col": "HRD_status", "title": "HRD status",
        "order": ["HRD", "HRP"],
        "labels": {"HRD": "HRD", "HRP": "HR proficient"},
    },
    {
        "key": "e", "file": "BRCA_status",
        "col": "BRCA_status", "title": "BRCA status",
        "order": ["mutated", "wildtype"],
        "labels": {"mutated": "Mutated", "wildtype": "Wildtype"},
    },
    {
        "key": "f", "file": "vazquez_garcia_annotation",
        "col": "vazquez_garcia_label", "title": "Vázquez-García et al. epithelial annotation",
        "order": [f"Cancer.cell.{i}" for i in range(1, 7)]
                 + ["Cycling.cancer.cell.1", "Cycling.cancer.cell.2",
                    "Ciliated.cell.1", "Ciliated.cell.2"],
        "labels": {**{f"Cancer.cell.{i}": f"Cancer\ncell {i}" for i in range(1, 7)},
                   "Cycling.cancer.cell.1": "Cycling\ncancer 1",
                   "Cycling.cancer.cell.2": "Cycling\ncancer 2",
                   "Ciliated.cell.1": "Ciliated\ncell 1", "Ciliated.cell.2": "Ciliated\ncell 2"},
    },
    {
        "key": "g", "file": "geistlinger_annotation",
        "col": "geistlinger_label", "title": "Geistlinger et al. subtype",
        "order": ["DIF", "PRO", "MES", "IMR"],
        "labels": {"DIF": "Differen-\ntiated", "PRO": "Prolif-\nerative",
                   "MES": "Mesen-\nchymal", "IMR": "Immuno-\nreactive"},
    },
]

# ============================================================================
# COMPUTE COMPOSITIONS (pooled % of cells per group)
# ============================================================================

print("\nComputing compositions...", flush=True)
count_rows = []
for panel in PANELS:
    col = panel["col"]
    present = set(meta[col].unique())
    extra = sorted(present - set(panel["order"]) - {"NA", "Unassigned"})
    assert not extra, f"{col}: unplotted groups {extra}"
    order = [g_ for g_ in panel["order"] if g_ in present]
    panel["order"] = order
    records = []
    for grp in order:
        mask = meta[col] == grp
        n = int(mask.sum())
        row = {"group": grp, "n": n}
        vc = meta.loc[mask, "schema_nmf"].value_counts()
        for ep in EPI_ORDER:
            row[ep] = vc.get(ep, 0) / n * 100
            count_rows.append({"panel": panel["key"], "metadata": col, "group": grp,
                               "n_cells": n, "epitype": ep, "n_epitype": int(vc.get(ep, 0)),
                               "pct": row[ep]})
        records.append(row)
    panel["data"] = pd.DataFrame(records)
    n_excl = int((~meta[col].isin(order)).sum())
    panel["n_excluded"] = n_excl
    print(f"  {panel['key']} {panel['title']}: {len(order)} groups, "
          f"{panel['data']['n'].sum():,} cells plotted ({n_excl:,} NA/unassigned/other-study)")
pd.DataFrame(count_rows).to_csv(OUT_COUNTS, index=False)

# ============================================================================
# PLOTTING
# ============================================================================


def draw_panel(ax, panel, show_ylabel=True, show_title=True):
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
            if v > 8:
                ax.text(x_pos[i], bottom[i] + v / 2, f"{v:.0f}%", fontsize=FN - 1,
                        ha="center", va="center", color="white", fontweight="bold", alpha=0.9)
        bottom += vals

    x_labels = [f"{labels_map.get(g_, g_)}\n(n={int(n):,})" for g_, n in zip(df["group"], df["n"])]
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, fontsize=(FN - 0.5) if n_bars <= 10 else (FN - 1.3), linespacing=0.9)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_xlim(-0.6, n_bars - 0.4)

    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    if show_ylabel:
        ax.set_ylabel("Epithelial composition (%)", fontsize=FA)
    else:
        ax.set_yticklabels([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.spines["left"].set_linewidth(0.4)
    ax.tick_params(axis="y", length=2, width=0.4)
    if show_title:
        ax.text(0.5, 1.04, panel["title"], transform=ax.transAxes, fontsize=FA + 0.5,
                fontweight="bold", ha="center", va="bottom")


def legend_handles():
    return [Patch(facecolor=EPI_PALETTE[ep], edgecolor="white", linewidth=0.3, label=ep)
            for ep in EPI_ORDER]


print("Plotting standalone panels...", flush=True)
for panel in PANELS:
    n_bars = len(panel["data"])
    w = 0.9 + 0.36 * n_bars
    fig, ax = plt.subplots(figsize=(w, 2.3))
    draw_panel(ax, panel, show_ylabel=True, show_title=True)
    ax.legend(handles=legend_handles(), loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, fontsize=FN, handlelength=1.0, handleheight=0.8,
              handletextpad=0.4, labelspacing=0.3, borderaxespad=0.0)
    for ext in ("svg", "png"):
        out = os.path.join(OUT_DIR, f"SF9{panel['key']}_atlas_epitype_by_{panel['file']}.{ext}")
        fig.savefig(out, format=ext, dpi=450, bbox_inches="tight", pad_inches=0.04)
    print(f"  Saved: SF9{panel['key']}_atlas_epitype_by_{panel['file']}.{{svg,png}}")
    plt.close(fig)
