#!/usr/bin/env python3
"""
Atlas 01 — Step 04: post-integration QC review

PURPOSE
    Assess data integrity and integration quality of the post-filter atlas object
    before downstream annotation: per-study QC distributions (UMI / genes / Scrublet),
    and integration-quality UMAPs. Reads obs + UMAP only (never the full matrix).

INPUTS
    obj("atlas_scanvi")  = hgsc_atlas_scanvi.h5ad  (written by 03b).

OUTPUTS
    output_root/01_preprocess_qc/04_qc_review/atlas_qc_summary.csv, sample_qc_summary.csv
    output_root/01_preprocess_qc/04_qc_review/04*.svg  (QC histograms + UMAPs)

MANUSCRIPT PANEL(S)
    QC backend; substrate for SF1 (QC violins) is extracted downstream in figures/.

RUNTIME TIER
    moderate (backed read of obs + UMAP; no matrix materialisation).
"""

import gc
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

warnings.filterwarnings("ignore")
np.random.seed(SEED)

# ============================================================================
# PATHS (resolved via central config)
# ============================================================================

H5AD_PATH = obj("atlas_scanvi")
FIG_DIR   = path("output_root", "01_preprocess_qc", "04_qc_review")
OUT_DIR   = FIG_DIR
os.makedirs(FIG_DIR, exist_ok=True)

# ============================================================================
# STYLE — Cook Lab v1.2
# ============================================================================

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "axes.linewidth":  0.6,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.dpi":      450,
    "savefig.dpi":     450,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
})

KELLY_22 = [
    "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280",
    "#848482", "#008856", "#E68FAC", "#0067A5", "#F99379",
    "#604E97", "#F6A600", "#B3446C", "#882D17", "#8DB600",
    "#654522", "#E25822", "#2B3D26", "#CC79A7", "#56B4E9",
    "#009E73", "#D55E00",
]

STUDY_ORDER = [
    "denisenko_2024", "geistlinger_2020", "hornburg_2021",
    "loret_2022",     "luo_2024",         "nath_2021",
    "olalekan_2021",  "olbrecht_2021",    "regner_2021",
    "vazquez_garcia_2022", "xu_2022",     "zhang_2022",
    "zheng_2023",
]


def make_study_palette(studies):
    return {s: KELLY_22[i % len(KELLY_22)] for i, s in enumerate(studies)}


# ============================================================================
# LOAD DATA (backed → lightweight extraction → close)
# ============================================================================

print("=" * 60)
print("Step 4 — Atlas QC Review")
print("=" * 60)
print(f"\nLoading (backed): {H5AD_PATH}", flush=True)

adata = sc.read_h5ad(H5AD_PATH, backed="r")
print(f"  Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

OBS_COLS = [
    "study", "sample_id", "patient_id",
    "n_genes_by_counts", "total_counts", "doublet_score_scrublet",
    "celltype_pred",
]
obs = adata.obs[OBS_COLS].copy()
umap_coords = np.array(adata.obsm["X_umap"])
obs["UMAP1"] = umap_coords[:, 0]
obs["UMAP2"] = umap_coords[:, 1]

adata.file.close()
del adata, umap_coords
gc.collect()

for col in ["n_genes_by_counts", "total_counts", "doublet_score_scrublet"]:
    obs[col] = pd.to_numeric(obs[col], errors="coerce")

studies_present = [s for s in STUDY_ORDER if s in obs["study"].unique()]
extras = sorted(set(obs["study"].unique()) - set(STUDY_ORDER))
studies_present += extras
STUDY_PALETTE = make_study_palette(studies_present)
print(f"Studies ({len(studies_present)}): {studies_present}\n")


# ============================================================================
# SUMMARY STATISTICS → CSV
# ============================================================================

summary_rows = []
for study in studies_present:
    sub = obs[obs["study"] == study]
    summary_rows.append({
        "study":              study,
        "n_cells":            len(sub),
        "n_samples":          sub["sample_id"].nunique(),
        "n_patients":         sub["patient_id"].nunique(),
        "median_total_counts":    sub["total_counts"].median(),
        "median_n_genes":         sub["n_genes_by_counts"].median(),
        "median_doublet_score":   sub["doublet_score_scrublet"].median(),
        "max_doublet_score":      sub["doublet_score_scrublet"].max(),
        "pct_doublet_score_gt02": (sub["doublet_score_scrublet"] > 0.2).mean() * 100,
    })
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUT_DIR, "atlas_qc_summary.csv"), index=False)

sample_counts = (
    obs.groupby(["study", "sample_id"])
    .size()
    .reset_index(name="n_cells")
    .sort_values(["study", "n_cells"], ascending=[True, False])
)
sample_counts.to_csv(os.path.join(OUT_DIR, "sample_qc_summary.csv"), index=False)


# ============================================================================
# HELPER: histogram grid (one panel per study)
# ============================================================================

def plot_hist_grid(obs_df, col, xlabel, studies, palette,
                   bins=50, xlog=False, xlim=None, fname=None, title=None):
    n = len(studies)
    ncols = 4
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.8, nrows * 2.0),
                             constrained_layout=True)
    axes = np.array(axes).flatten()
    for i, study in enumerate(studies):
        ax = axes[i]
        vals = obs_df.loc[obs_df["study"] == study, col].dropna()
        color = palette[study]
        if xlog:
            vals = vals[vals > 0]
            ax.set_xscale("log")
            bins_use = np.logspace(np.log10(vals.min()), np.log10(vals.max()), bins)
        else:
            bins_use = bins
        ax.hist(vals, bins=bins_use, color=color, alpha=0.85, linewidth=0)
        med = vals.median()
        ax.axvline(med, color="#333333", linewidth=0.8, linestyle="--")
        ax.text(0.97, 0.95, f"med={med:,.0f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=5.5, color="#333333")
        ax.text(0.97, 0.80, f"n={len(vals):,}", transform=ax.transAxes,
                ha="right", va="top", fontsize=5.5, color="#555555")
        ax.set_title(study.replace("_", "\n"), fontsize=6.5, pad=3)
        ax.set_xlabel(xlabel, fontsize=6)
        ax.set_ylabel("Cells", fontsize=6)
        ax.tick_params(labelsize=5.5)
        ax.spines[["top", "right"]].set_visible(False)
        if xlim is not None:
            ax.set_xlim(xlim)
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda x, _: f"{int(x):,}" if x >= 1000 else f"{int(x)}")
        )
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)
    if title:
        fig.suptitle(title, fontsize=10, fontweight="bold", y=1.01)
    fig.savefig(os.path.join(FIG_DIR, fname), format="svg", dpi=450, bbox_inches="tight")
    plt.close(fig)


plot_hist_grid(obs, "total_counts", xlabel="Total UMI counts",
               studies=studies_present, palette=STUDY_PALETTE, bins=60, xlog=True,
               fname="04a_hist_total_counts.svg",
               title="Total UMI counts per cell, by study")
plot_hist_grid(obs, "n_genes_by_counts", xlabel="Genes detected",
               studies=studies_present, palette=STUDY_PALETTE, bins=60, xlog=False,
               fname="04b_hist_n_genes.svg",
               title="Genes detected per cell, by study")
plot_hist_grid(obs, "doublet_score_scrublet", xlabel="Scrublet doublet score",
               studies=studies_present, palette=STUDY_PALETTE, bins=50,
               xlim=(0, 0.26), fname="04c_hist_doublet_score.svg",
               title="Scrublet doublet scores per cell, by study (all cells < 0.25)")


# ============================================================================
# UMAP HELPER
# ============================================================================

def save_umap_svg(obs_df, color_col, palette, fname, title,
                  continuous=False, cmap="viridis", pt_size=0.5):
    fig, ax = plt.subplots(figsize=(5, 4))
    data = obs_df.dropna(subset=[color_col]).copy()
    data = data.sample(frac=1, random_state=SEED)
    if continuous:
        data[color_col] = pd.to_numeric(data[color_col], errors="coerce")
        data = data.dropna(subset=[color_col]).sort_values(color_col)
        sc_plot = ax.scatter(data["UMAP1"], data["UMAP2"], c=data[color_col],
                             cmap=cmap, s=pt_size, alpha=0.6, linewidths=0, rasterized=True)
        cb = fig.colorbar(sc_plot, ax=ax, shrink=0.6, pad=0.02)
        cb.ax.tick_params(labelsize=6)
        cb.set_label(color_col.replace("_", " ").title(), fontsize=6)
    else:
        categories = sorted(data[color_col].astype(str).unique())
        for cat in categories:
            mask = data[color_col].astype(str) == cat
            ax.scatter(data.loc[mask, "UMAP1"], data.loc[mask, "UMAP2"],
                       c=palette.get(cat, "#A0A0A0"), s=pt_size, alpha=0.6,
                       linewidths=0, rasterized=True, label=cat)
        handles = [
            Line2D([0], [0], marker="o", color="w",
                   markerfacecolor=palette.get(cat, "#A0A0A0"),
                   markersize=4, linewidth=0, label=cat)
            for cat in categories
        ]
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1),
                  ncol=1, fontsize=5.5, frameon=False, markerscale=1.2, handletextpad=0.3)
    ax.set_title(title, fontsize=9, pad=6)
    ax.set_axis_off()
    fig.savefig(os.path.join(FIG_DIR, fname), format="svg", dpi=450, bbox_inches="tight")
    plt.close(fig)


save_umap_svg(obs, "study", palette=STUDY_PALETTE, fname="04d_umap_study.svg",
              title="UMAP — coloured by study", pt_size=0.5)
save_umap_svg(obs, "doublet_score_scrublet", palette=None,
              fname="04e_umap_doublet_score.svg", title="UMAP — Scrublet doublet score",
              continuous=True, cmap="viridis", pt_size=0.5)

SMALL_PT = 0.02
save_umap_svg(obs, "study", palette=STUDY_PALETTE, fname="04f_umap_study_small.svg",
              title="UMAP — coloured by study (small points)", pt_size=SMALL_PT)
save_umap_svg(obs, "doublet_score_scrublet", palette=None,
              fname="04g_umap_doublet_score_small.svg",
              title="UMAP — Scrublet doublet score (small points)",
              continuous=True, cmap="viridis", pt_size=SMALL_PT)

CELLTYPE_PALETTE = {
    "Epithelial":    "#E6A141",
    "Mesothelial":   "#A8A298",
    "Fibroblast":    "#DDD5CA",
    "Smooth_Muscle": "#D14E6C",
    "Pericyte":      "#B87A7A",
    "Endothelial":   "#7D4E4E",
    "T_cell":        "#87CEFA",
    "NK_cell":       "#56AFC4",
    "B_cell":        "#5665B6",
    "Plasma_cell":   "#8A5DAF",
    "Macrophage":    "#8FBC8F",
    "DC":            "#2E8B57",
    "Neutrophil":    "#6B8E23",
    "Mast":          "#8B9B6B",
    "Erythrocyte":   "#CD5C5C",
}
save_umap_svg(obs, "celltype_pred", palette=CELLTYPE_PALETTE,
              fname="04h_umap_celltype_pred_small.svg",
              title="UMAP — celltype_pred (small points)", pt_size=SMALL_PT)

print("Done — QC review outputs written to output_root/01_preprocess_qc/04_qc_review/")
