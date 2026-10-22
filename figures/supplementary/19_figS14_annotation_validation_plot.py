#!/usr/bin/env python3
"""
Figure S14 — Cell type annotation validation within the Xenium cohort
=====================================================================

PURPOSE
    Canonical-marker dot plot across the 18 Xenium cell types (colour = per-gene
    min-max scaled mean log-expression, size = % expressing), in the style of the
    Figure S4b atlas dot plot. Genes excluded from SingleR by probe QC are marked
    (dagger); Xenium-only genes (absent from the atlas reference, hence also never
    seen by SingleR) are marked with a double dagger.

INPUTS
    output_root/_prep_caches/sf12/<sample>_agg_{groups,sumlog,nexpr}   (from
        18_figS14_annotation_validation_extract.R: per-group logcounts sums /
        expressing counts over ALL cells)
    obj("xenium_ref") = xenium_celltype_downsampled.h5ad (SingleR reference, 16,000
        cells / 16 types; its gene list defines the Xenium-only genes)
    output_root/05_probe_qc/genes_exclude_singler.txt (28 genes)
    the final TMA core set (figures/_prep/final_tma_cores.py)

COHORT
    8 whole tissues (core 8) + TMA A (HGSC + FTE cores) + TMA B (blocks G/H).
    HGSC TMA: the final TMA core set (spatial/02_qc/06b_final_tma_cores.py), as every
    TMA panel -> 378 patients / 503 cores. FTE: patient x tissue-type with >= 500 cells (all 15 cores).
    -> 8 WT / 378 HGSC TMA / 15 FTE. All cells of those units used; no subsampling.

OUTPUTS
    output_root/figures/supplementary/SF12B_xenium_canonical_markers_dotplot.{svg,png}
        (Figure S14; promoted as SuppFig14b_xenium_canonical_markers_dotplot)
    output_root/_prep_caches/sf12/sf12b_dotplot_stats.csv (source data)
    Output and cache names carry the earlier figure number (SF12 = Figure S14).

MANUSCRIPT PANEL(S)
    Figure S14.

RUNTIME TIER
    fast (reads the per-group aggregates and the 16,000-cell reference h5ad).
"""

import os
import sys

import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, REPO)
from config.config import obj, path, SEED  # noqa: E402
from figures._prep.final_tma_cores import final_core_ids  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS
# ============================================================================

CACHE = path("output_root", "_prep_caches", "sf12")
EXCL = path("output_root", "05_probe_qc", "genes_exclude_singler.txt")
REF = obj("xenium_ref")
if not os.path.exists(REF):  # DATA_ROOT not pointed at the bundle -> canonical copy
    REF = os.path.join(os.environ.get("DATA_ROOT", "/path/to/deposit"), "2026_final_atlas/output/xenium_celltype/xenium_celltype_downsampled.h5ad")
OUT_DIR = path("output_root", "figures", "supplementary")
os.makedirs(OUT_DIR, exist_ok=True)
STEMS = {
    # only "b" is drawn (Figure S14, promoted as SuppFig14b_*); the "a" and "c" names are not written
    "a": "SF12A_xenium_singler_confidence",
    "b": "SF12B_xenium_canonical_markers_dotplot",
    "c": "SF12C_xenium_atlas_agreement",
}

WT8 = ["WT01", "WT02", "WT03", "WT04", "WT05", "WT06",
       "WT07", "WT08"]
TMA = ["tma", "COEUR_TMA_G", "COEUR_TMA_H"]
SAMPLES = TMA + WT8
MIN_CELLS = 500
FINAL = {s: final_core_ids(s) for s in TMA}   # final HGSC TMA cores per array file (raw core ids)

# ============================================================================
# STYLE (Figure S4 / Fig 4 conventions)
# ============================================================================

FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "xtick.labelsize": FK, "ytick.labelsize": FK,
    "legend.fontsize": FN, "figure.dpi": 450, "savefig.dpi": 450,
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.bbox": "tight",
})

# Fig 4 / Fig 5 palette (figures/figure5/11_fig5_compose_panels.py PAL == ref_palette in 00_setup.R)
PAL = {"SecA epithelium": "#E6A141", "Intermediate epithelium": "#C08E48", "SecB epithelium": "#9A7D55",
       "Ciliated epithelium": "#E07850", "Mesothelial": "#D4A574", "Fibroblast": "#DDD5CA",
       "Smooth muscle": "#D14E6C", "Pericyte": "#B87A7A", "Endothelial": "#7D4E4E",
       "Macrophage": "#8FBC8F", "Mast cell": "#8B9B6B", "Neutrophil": "#6B8E23",
       "Conventional dendritic cell": "#2E8B57", "Plasmacytoid dendritic cell": "#A0A0A0",
       "NK cell": "#56AFC4", "T cell": "#87CEFA", "Plasma cell": "#8A5DAF", "B cell": "#5665B6"}
ORDER = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium", "Ciliated epithelium",
         "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
         "T cell", "NK cell", "B cell", "Plasma cell",
         "Macrophage", "Conventional dendritic cell", "Plasmacytoid dendritic cell",
         "Neutrophil", "Mast cell"]
DISPLAY = {"SecA epithelium": "SecA epithelium", "Intermediate epithelium": "Intermediate epithelium",
           "SecB epithelium": "SecB epithelium", "Ciliated epithelium": "Ciliated epithelium",
           "Conventional dendritic cell": "Conventional DC", "Plasmacytoid dendritic cell": "Plasmacytoid DC"}


def disp(ct):
    return DISPLAY.get(ct, ct)


def save(fig, key, pad=0.04):
    for ext in ("svg", "png"):
        fp = os.path.join(OUT_DIR, f"{STEMS[key]}.{ext}")
        fig.savefig(fp, format=ext, dpi=450, bbox_inches="tight", pad_inches=pad)
        print(f"  Saved: {fp}")
    plt.close(fig)


def style_axes(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="both", which="both", length=2, width=0.4)


def tissue_of(sample, sample_type):
    if sample in WT8:
        return "HGSC whole tissue"
    return "Fallopian tube" if sample_type == "fallopian" else "HGSC TMA"


def pid_of(sample, patient_id):
    if sample in WT8:
        return sample
    if pd.isna(patient_id):
        return np.nan
    p = str(patient_id).replace(".0", "") if str(patient_id).endswith(".0") else str(patient_id)
    return ("COEUR-" + p) if sample.startswith("COEUR") else p


# ============================================================================
# 1. COHORT (group aggregates; >=500 cells per patient x tissue type)
# ============================================================================

print("[1] Loading per-group aggregates ...", flush=True)
G, SUM, NEX = [], [], []
for s in SAMPLES:
    g = pd.read_csv(os.path.join(CACHE, f"{s}_agg_groups.csv"),
                    dtype={"patient_id": str, "sample_type": str, "core_id": str})
    G.append(g)
    SUM.append(pd.read_csv(os.path.join(CACHE, f"{s}_agg_sumlog.csv.gz")).values)
    NEX.append(pd.read_csv(os.path.join(CACHE, f"{s}_agg_nexpr.csv.gz")).values)
GENES = list(pd.read_csv(os.path.join(CACHE, f"{SAMPLES[0]}_agg_sumlog.csv.gz"), nrows=1).columns)
grp = pd.concat(G, ignore_index=True)
SUM, NEX = np.vstack(SUM), np.vstack(NEX)
grp["tissue"] = [tissue_of(s, t) for s, t in zip(grp["sample"], grp["sample_type"])]
grp["pid"] = [pid_of(s, p) for s, p in zip(grp["sample"], grp["patient_id"])]
has_pid = grp["pid"].notna()
tot = grp[has_pid].groupby(["pid", "tissue"]).n_cells.sum()
ok_pt = set(tot[tot >= MIN_CELLS].index)


def in_cohort(sample, tissue, pid, core):
    """HGSC TMA: final TMA core set; WT and FTE: patient x tissue with >= 500 cells"""
    if tissue == "HGSC TMA":
        return isinstance(core, str) and core in FINAL[sample]
    return isinstance(pid, str) and (pid, tissue) in ok_pt


grp["keep"] = [in_cohort(s, t, p, c) for s, t, p, c in zip(grp["sample"], grp["tissue"], grp["pid"], grp["core_id"])]
n_all = int(grp.n_cells.sum())
n_keep = int(grp.loc[grp.keep, "n_cells"].sum())
pts = grp[grp.keep].groupby("tissue").pid.nunique()
cells_by_tissue = grp[grp.keep].groupby("tissue").n_cells.sum()
print(f"  cells: {n_keep:,} kept of {n_all:,} annotated; patients: {pts.to_dict()}")
print(f"  cells by tissue: {cells_by_tissue.to_dict()}")

k = grp.keep.values
gk = grp[k].reset_index(drop=True)
Sk, Nk = SUM[k], NEX[k]
idx = gk.groupby("cell_label").indices
n_ct = pd.Series({ct: int(gk.n_cells.values[idx[ct]].sum()) for ct in ORDER})
x_mean = pd.DataFrame([Sk[idx[ct]].sum(0) / n_ct[ct] for ct in ORDER], index=ORDER, columns=GENES)
x_frac = pd.DataFrame([Nk[idx[ct]].sum(0) / n_ct[ct] for ct in ORDER], index=ORDER, columns=GENES)
assert n_ct.sum() == n_keep, "cell_label not covering all cohort cells"

# ============================================================================
# 4. FIGURE S14 — canonical-marker dot plot (Figure S4b styling; drawn as panel "b")
# ============================================================================

print("[4] Panel B ...", flush=True)
excluded = set(g.strip() for g in open(EXCL) if g.strip())
ref = ad.read_h5ad(REF)
ref_genes = set(ref.var_names)
xen_only = set(g for g in GENES if g not in ref_genes)

# Two canonical markers per cell type (Figure S4b logic: grouped and ordered to match the
# y-axis), restricted to the Xenium panel. Figure S4b genes are kept wherever the panel has them.
# Epithelium: pan-epithelial EPCAM/PAX8 (Figure S4b), then the SecA / SecB signature genes
# (shared/signatures.yml); Intermediate cells are a continuum between SecA and SecB and
# have no discrete marker pair. Where Figure S4b genes are absent from the panel the replacement is
# the top panel gene for that type in the atlas reference (e.g. Mesothelial C3: UPK3B/CALB2 are
# not on the panel; Neutrophil MXD1: S100A8 absent and the FCGR3B probe is macrophage-dominated).
MARKERS = [
    ("Epithelial",                   ["EPCAM", "PAX8"]),
    ("SecA epithelium",              ["MECOM", "SOX17"]),
    ("SecB epithelium",              ["KRT7", "SLPI"]),
    ("Ciliated epithelium",          ["FOXJ1", "TPPP3"]),
    ("Mesothelial",                  ["MSLN", "C3"]),
    ("Fibroblast",                   ["DCN", "PDGFRA"]),
    ("Smooth muscle",                ["MYH11", "ACTA2"]),
    ("Pericyte",                     ["NOTCH3", "TPM2"]),
    ("Endothelial",                  ["VWF", "PECAM1"]),
    ("T cell",                       ["CD3E", "CD3D"]),
    ("NK cell",                      ["KLRD1", "PRF1"]),
    ("B cell",                       ["MS4A1", "CD79A"]),
    ("Plasma cell",                  ["MZB1", "XBP1"]),
    ("Macrophage",                   ["CD14", "C1QC"]),
    ("Conventional dendritic cell",  ["FLT3", "IDO1"]),
    ("Plasmacytoid dendritic cell",  ["PLD4", "GZMB"]),
    ("Neutrophil",                   ["CSF3R", "MXD1"]),
    ("Mast cell",                    ["KIT", "CPA3"]),
]
GENE_ORDER = [g for _, gs in MARKERS for g in gs]
missing = [g for g in GENE_ORDER if g not in GENES]
assert not missing, missing
bounds = np.cumsum([len(gs) for _, gs in MARKERS])

mu = x_mean[GENE_ORDER]
norm = (mu - mu.min()) / (mu.max() - mu.min())
fr = x_frac[GENE_ORDER]
long = []
for ct in ORDER:
    for gname in GENE_ORDER:
        long.append({"celltype": ct, "gene": gname, "frac_expressing": fr.loc[ct, gname],
                     "mean_expression": mu.loc[ct, gname], "norm_expression": norm.loc[ct, gname],
                     "n_cells": int(n_ct[ct]), "excluded_from_singler_probe_qc": gname in excluded,
                     "xenium_only": gname in xen_only})
pd.DataFrame(long).to_csv(os.path.join(CACHE, "sf12b_dotplot_stats.csv"), index=False)

LAPAZ_CMAP = LinearSegmentedColormap.from_list(
    "lapaz_r", [(0.00, "#F5F5F5"), (0.20, "#F0D0B0"), (0.40, "#D4916B"),
                (0.60, "#B05A5A"), (0.80, "#7A2040"), (1.00, "#4A0C2A")])
n_r, n_g = len(ORDER), len(GENE_ORDER)
# keep Figure S4b's per-dot spacing (7.5 in for 26 genes, 2.2 in for 13 rows)
fig, ax = plt.subplots(figsize=(7.5 * (n_g + 6) / (26 + 6), 2.2 * (n_r + 1.5) / (13 + 1.5)))
DOT_SCALE = 120
X, Y, Sz, Cc = [], [], [], []
for i, ct in enumerate(ORDER):
    for j, gname in enumerate(GENE_ORDER):
        X.append(j); Y.append(n_r - 1 - i)
        Sz.append(fr.loc[ct, gname] * DOT_SCALE); Cc.append(norm.loc[ct, gname])
ax.scatter(X, Y, s=Sz, c=Cc, cmap=LAPAZ_CMAP, vmin=0, vmax=1, edgecolors="black",
           linewidths=0.3, zorder=3)
for b in bounds[:-1]:
    ax.axvline(b - 0.5, color="#DDDDDD", linewidth=0.4, zorder=1)
# Epithelial+mesothelial | stromal | lymphoid | myeloid
for b in (5, 9, 13):
    ax.axhline(n_r - b - 0.5, color="#DDDDDD", linewidth=0.4, zorder=1)


def gl(gname):
    return gname + ("†" if gname in excluded else "") + ("‡" if gname in xen_only else "")


ax.set_xticks(range(n_g))
ax.set_xticklabels([gl(g) for g in GENE_ORDER], rotation=45, ha="right", fontsize=FA, style="italic")
ax.set_yticks(range(n_r))
ax.set_yticklabels([disp(ct) for ct in reversed(ORDER)], fontsize=FA)
for i, ct in enumerate(reversed(ORDER)):
    lbl = ax.get_yticklabels()[i]
    lbl.set_color(PAL[ct]); lbl.set_fontweight("bold")
    ax.add_patch(plt.Rectangle((-1.1, i - 0.3), 0.5, 0.6, facecolor=PAL[ct], edgecolor="none",
                               clip_on=False, zorder=5))
ax.set_xlim(-0.5, n_g - 0.5)
ax.set_ylim(-0.5, n_r - 0.5)
style_axes(ax)
ax.tick_params(axis="y", pad=10)
size_vals = [0.2, 0.4, 0.6, 0.8]
size_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor="#888888",
                       markeredgecolor="black", markeredgewidth=0.3,
                       markersize=np.sqrt(v * DOT_SCALE), linewidth=0, label=f"{int(v * 100)}%")
                for v in size_vals]
leg = ax.legend(handles=size_handles, title="% expressing", loc="upper left",
                bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=FN, title_fontsize=FA,
                handletextpad=0.4, labelspacing=0.8, borderpad=0.2)
ax.add_artist(leg)
sm = plt.cm.ScalarMappable(cmap=LAPAZ_CMAP, norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.25, pad=0.10, aspect=10, anchor=(0.0, 0.0))
cbar.ax.tick_params(labelsize=FN)
cbar.set_label("Norm. expression", fontsize=FN)
ax.text(1.02, 0.40, "† excluded from SingleR\n   (probe QC)\n‡ Xenium-only (not in\n   atlas reference)",
        transform=ax.transAxes, ha="left", va="top", fontsize=FN - 0.5, color="#555555")
save(fig, "b", pad=0.06)

print("SF12_DONE")
