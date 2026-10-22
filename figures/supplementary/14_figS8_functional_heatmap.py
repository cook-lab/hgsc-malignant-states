#!/usr/bin/env python3
"""
Figure S8 — Functional characterization heatmaps (4 blocks)
===========================================================

Purpose
    Four z-scored activity heatmaps across the four epitypes:
      PROGENy signalling pathways (14 features)
      MSigDB Hallmark gene sets (50 features)
      Metabolic flux modules (scFEA on all 575,366 epithelial cells; 165 features = 168 minus the
         3 constant modules M_20, M_83, M_128)
      DoRoThEA transcription-factor activity (294 features)
    Features (hierarchically clustered) on y-axis, epitypes on x-axis, colour =
    row-normalised z-score, flowing across 4 equal columns on a letter page.

INPUTS
    data_root/2026_final_atlas/output/21_epitype_functional_characterization/
        21_{progeny,hallmark,flux,dorothea}_per_cell.parquet   (per-cell activities / scFEA flux)
    figures._prep.atlas_final_tables.load_epithelial()      (epithelial `state`, final atlas)
        -> per-state mean of each feature -> z-score across the 4 states (pandas std, ddof=1;
           the upstream 21_epitype_functional_characterization.py zscore_df), rebuilt here into
           output_root/_prep_caches/sf8_functional/21_{axis}_{means,zscored}.csv.
    output_root/figure2/fig2a_functional/21_{axis}_zscored.csv  (the same tables rebuilt by the
        Figure 2 pipeline for Fig 2a). When present they are used (one source for Fig 2a and Figure S8)
        after checking they agree with the rebuild above (max |diff| <= 1e-6); otherwise the
        rebuild is used.
    scFEA modules whose per-state means are exactly constant (zero variance; M_20, M_83 and
    M_128 on all cells) are dropped (their z-scores are 0 by construction).
    data_root/2026_final_atlas/tools/scFEA/data/Human_M168_information.symbols.csv
        (scFEA module-id -> human-readable Compound_IN -> Compound_OUT label map;
         resolves the int-Module_id vs M_x label-map quirk — convention #8)

OUTPUTS
    output_root/figures/supplementary/SF8_atlas_functional_heatmap.{svg,png}
    output_root/_prep_caches/sf8_functional/21_{axis}_{means,zscored}.csv (the rebuilt per-state values)

MANUSCRIPT PANEL(S)
    Figure S8.

RUNTIME TIER
    fast (renders from pre-computed z-score CSVs).
"""

import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, OUTPUT_ROOT  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

PER_CELL_DIR = path("data_root", "2026_final_atlas", "output", "21_epitype_functional_characterization")
SF8_CACHE = path("output_root", "_prep_caches", "sf8_functional")
FIG2_DIR = os.path.join(str(OUTPUT_ROOT), "figure2", "fig2a_functional")   # not created here
MODULE_INFO = path("data_root", "2026_final_atlas", "tools", "scFEA", "data",
                   "Human_M168_information.symbols.csv")
OUT_SVG = path("output_root", "figures", "supplementary", "SF8_atlas_functional_heatmap.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF8_atlas_functional_heatmap.png")

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

CMAP = LinearSegmentedColormap.from_list(
    "bwr_custom",
    ["#2166AC", "#67A9CF", "#D1E5F0", "#FAFAFA", "#FDDBC7", "#EF8A62", "#B2182B"],
    N=256,
)

# Epitype display order & colours.
EPI_ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]
EPI_DISPLAY = ["SecA", "Int.", "SecB", "Ciliated"]
EPI_PALETTE = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55", "Ciliated": "#E07850"}

# ============================================================================
# LOAD DATA
# ============================================================================

print("Rebuilding per-state z-scores from the per-cell outputs...", flush=True)
STATE = load_epithelial(columns=["state"])["state"].astype(str)
CONSTANT_FLUX = []


def rebuild(axis):
    src = path("output_root", "figure2", "fig2a_scfea_all") if axis == "flux" else PER_CELL_DIR
    pc = pd.read_parquet(os.path.join(src, f"21_{axis}_per_cell.parquet"))  # flux: all cells (2026-10-05)
    pc.index = pc.index.astype(str)
    lab = STATE.reindex(pc.index)
    assert lab.notna().all(), f"{axis}: per-cell barcodes missing from the epithelial table"
    means = pc.groupby(lab.values).mean().reindex(EPI_ORDER)
    means.index.name = "group"
    z = (means - means.mean(axis=0)) / means.std(axis=0).replace(0, 1)
    means.to_csv(os.path.join(SF8_CACHE, f"21_{axis}_means.csv"))
    z.to_csv(os.path.join(SF8_CACHE, f"21_{axis}_zscored.csv"))
    print(f"  {axis}: {len(pc):,} cells x {pc.shape[1]} features; per state "
          f"{lab.value_counts().reindex(EPI_ORDER).to_dict()}")
    if axis == "flux":
        CONSTANT_FLUX.extend(means.columns[means.std(axis=0) == 0].tolist())
    return z


def load_zscored(fname, label_clean_fn=None):
    """Rebuilt z-scores (or the identical Fig 2 tables) -> (feature_names, matrix [n_features × 4])."""
    axis = fname.split("_")[1]
    df = rebuild(axis)
    f2 = os.path.join(FIG2_DIR, fname)
    if os.path.exists(f2):
        d2 = pd.read_csv(f2).set_index("group").reindex(EPI_ORDER)
        diff = float(np.nanmax(np.abs(d2[df.columns].values - df.values))) if set(df.columns) <= set(d2.columns) else np.inf
        if diff <= 1e-6:
            print(f"  {fname}: using Figure 2 table {f2} (max |diff| vs rebuild {diff:.1e})")
            df = d2[df.columns]
        else:
            print(f"  [WARN] {fname}: Figure 2 table differs from the rebuild (max |diff| {diff}); using the rebuild")
    if axis == "flux":
        df = df.drop(columns=CONSTANT_FLUX)
        print(f"  flux: dropped {len(CONSTANT_FLUX)} constant scFEA modules {CONSTANT_FLUX}")
    mat = df.values.T
    features = list(df.columns)
    if label_clean_fn:
        features = [label_clean_fn(f) for f in features]
    print(f"  {fname}: {mat.shape[0]} features × {mat.shape[1]} epitypes")
    return features, mat


def clean_hallmark(name):
    name = name.replace("HALLMARK_", "").replace("_", " ").title()
    if len(name) > 30:
        name = name[:28] + "…"
    return name


# --- Metabolic module name mapping (scFEA M168) ---
assert os.path.exists(MODULE_INFO), f"Missing: {MODULE_INFO}"
mod_df = pd.read_csv(MODULE_INFO, index_col=0)
module_name_map = {}
for mid, row in mod_df.iterrows():
    cin = str(row["Compound_IN_name"]).strip()
    cout = str(row["Compound_OUT_name"]).strip()
    label = f"{cin} → {cout}"
    if len(label) > 40:
        label = label[:38] + "…"
    module_name_map[mid] = label
print(f"  Loaded {len(module_name_map)} metabolic module names")


def clean_flux(name):
    return module_name_map.get(name, name)


progeny_feat, progeny_mat = load_zscored("21_progeny_zscored.csv")
hallmark_feat, hallmark_mat = load_zscored("21_hallmark_zscored.csv", clean_hallmark)
flux_feat, flux_mat = load_zscored("21_flux_zscored.csv", clean_flux)
dorothea_feat, dorothea_mat = load_zscored("21_dorothea_zscored.csv")

# ============================================================================
# HIERARCHICAL CLUSTERING (reorder rows)
# ============================================================================


def cluster_rows(mat):
    if mat.shape[0] <= 1:
        return np.arange(mat.shape[0])
    Z = linkage(pdist(mat, metric="euclidean"), method="ward")
    return leaves_list(Z)


print("\nClustering features...", flush=True)
progeny_order = cluster_rows(progeny_mat)
hallmark_order = cluster_rows(hallmark_mat)
flux_order = cluster_rows(flux_mat)
dorothea_order = cluster_rows(dorothea_mat)

progeny_mat = progeny_mat[progeny_order]
progeny_feat = [progeny_feat[i] for i in progeny_order]
hallmark_mat = hallmark_mat[hallmark_order]
hallmark_feat = [hallmark_feat[i] for i in hallmark_order]
flux_mat = flux_mat[flux_order]
flux_feat = [flux_feat[i] for i in flux_order]
dorothea_mat = dorothea_mat[dorothea_order]
dorothea_feat = [dorothea_feat[i] for i in dorothea_order]

# ============================================================================
# PLOT — continuous heatmap flowing across 4 equal columns
# ============================================================================

print("\nPlotting...", flush=True)

PAGE_W, PAGE_H = 7.5, 10.0
section_names = ["PROGENy", "Hallmark", "Metabolic flux", "DoRoThEA"]
SECTION_COLORS = {
    "PROGENy": "#5B8FA8", "Hallmark": "#8B6DAF",
    "Metabolic flux": "#4A967A", "DoRoThEA": "#C47A3D",
}
panel_data = [
    (progeny_feat, progeny_mat),
    (hallmark_feat, hallmark_mat),
    (flux_feat, flux_mat),
    (dorothea_feat, dorothea_mat),
]

all_feats = []
all_mat_rows = []
all_section_ids = []
section_starts = []

for si, (feats, mat) in enumerate(panel_data):
    section_starts.append(len(all_feats))
    all_feats.extend(feats)
    all_mat_rows.append(mat)
    all_section_ids.extend([si] * len(feats))

all_mat = np.vstack(all_mat_rows)
total = len(all_feats)
print(f"  Total features: {total}")

N_COLS = 4
rows_per_col = int(np.ceil(total / N_COLS))
col_slices = []
for c in range(N_COLS):
    i0 = c * rows_per_col
    i1 = min(i0 + rows_per_col, total)
    col_slices.append((i0, i1))

VMIN, VMAX = -1.5, 1.5
fig = plt.figure(figsize=(PAGE_W, PAGE_H))
gs_outer = GridSpec(1, N_COLS, figure=fig, wspace=0.50,
                    left=0.10, right=0.98, top=0.96, bottom=0.05)

im = None
for ci, (i0, i1) in enumerate(col_slices):
    gs_inner = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[0, ci],
                                       width_ratios=[0.035, 1.0], wspace=0.03)

    ax_side = fig.add_subplot(gs_inner[0, 0])
    chunk_sids = all_section_ids[i0:i1]
    n_rows = len(chunk_sids)
    sidebar_colors = np.zeros((n_rows, 1, 3))
    for ri, sid in enumerate(chunk_sids):
        sidebar_colors[ri, 0, :] = to_rgb(SECTION_COLORS[section_names[sid]])
    ax_side.imshow(sidebar_colors, aspect="auto", interpolation="nearest")
    ax_side.set_xticks([]); ax_side.set_yticks([])
    for sp in ax_side.spines.values():
        sp.set_linewidth(0.3); sp.set_color("#888888")

    ax = fig.add_subplot(gs_inner[0, 1])
    chunk_mat = all_mat[i0:i1]
    chunk_feats = all_feats[i0:i1]
    im = ax.imshow(chunk_mat, aspect="auto", cmap=CMAP, vmin=VMIN, vmax=VMAX,
                   interpolation="nearest")

    ax.set_xticks(range(4))
    ax.set_xticklabels(EPI_DISPLAY, fontsize=FN, rotation=45, ha="right")
    ax.tick_params(axis="x", length=0, pad=5)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    for xi, ep in enumerate(EPI_ORDER):
        ax.plot(xi, 1.006, "s", color=EPI_PALETTE[ep], markersize=3,
                clip_on=False, transform=ax.get_xaxis_transform())
    # pin the y-range to the heatmap rows; the colour squares above are drawn off-axis, in axes
    # coordinates just above the heatmap.
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xlim(-0.5, 3.5)

    ax.yaxis.tick_right()
    ax.set_yticks(range(n_rows))
    label_fs = max(0.9, min(3.0, 180 / n_rows))
    ax.set_yticklabels(chunk_feats, fontsize=label_fs)
    ax.tick_params(axis="y", length=0, pad=1)

    for sec_start in section_starts:
        if i0 < sec_start < i1:
            local_y = sec_start - i0
            ax.axhline(local_y - 0.5, color="white", linewidth=1.2, zorder=5)
            ax_side.axhline(local_y - 0.5, color="white", linewidth=1.2, zorder=5)

    for sp in ax.spines.values():
        sp.set_linewidth(0.3); sp.set_color("#888888")

sec_handles = [Patch(facecolor=SECTION_COLORS[s], edgecolor="white", linewidth=0.3, label=s)
               for s in section_names]
fig.legend(handles=sec_handles, loc="lower left", ncol=4, frameon=False,
           fontsize=FN + 0.5, handlelength=1.0, handleheight=0.8, handletextpad=0.3,
           columnspacing=1.2, bbox_to_anchor=(0.10, 0.005))

cbar_ax = fig.add_axes([0.55, 0.018, 0.30, 0.007])
cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
cbar.ax.tick_params(labelsize=FN, length=2, width=0.4)
cbar.set_ticks([-1.5, -0.75, 0, 0.75, 1.5])
cbar.set_label("Z-score", fontsize=FN + 0.5, labelpad=2)
cbar.outline.set_linewidth(0.3)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.06)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.06)
plt.close(fig)

print(f"\n  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
