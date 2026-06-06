#!/usr/bin/env python3
"""
SF11 — Atlas vs Xenium whole-tissue polarization violins
========================================================

Purpose
    Violin strip: atlas epitypes (SecA / Intermediate / SecB) followed by each
    Xenium whole-tissue sample, y-axis = polarization UCell (SecB_UCell -
    SecA_UCell). TMA is intentionally excluded (plotted separately). Dashed
    reference lines at the atlas SecB and Intermediate p25 of polarization UCell.

INPUTS
    output_root/18_ucell_atlas/atlas_ucell_scores_fulllabels.csv
        (atlas secretory cells; UCell on the xenium-shared gene panel; labels =
         current SecB-NMF schema_nmf. Uses the canonical 18b scoring that matches
         the xenium noBCAM signature — convention #6, Q5 — NOT the 18b_v2 output.)
    output_root/06f_reclassification_polarization/reclassified_xenium_scores.csv
        (pooled 06f reclassification; whole-tissue samples only, sfe_tma dropped)

OUTPUTS
    output_root/figures/supplementary/SF11_atlas_xenium_wt_polarization_violins.{png,svg}

MANUSCRIPT PANEL(S)
    SF11.

RUNTIME TIER
    fast (reads two score CSVs; pooled violins).
"""

import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, CFG  # noqa: E402

# ---------- Paths (central config) ----------
ATLAS_CSV = path("output_root", "18_ucell_atlas", "atlas_ucell_scores_fulllabels.csv")
XEN_CSV = path("output_root", "06f_reclassification_polarization", "reclassified_xenium_scores.csv")

STEM = "SF11_atlas_xenium_wt_polarization_violins"
OUT_PNG = path("output_root", "figures", "supplementary", f"{STEM}.png")
OUT_SVG = path("output_root", "figures", "supplementary", f"{STEM}.svg")

# ---------- Epitype + xenium sample colours ----------
# "Intermediate" (was "Transitioning").
EPI_ORDER = ["SecA", "Intermediate", "SecB"]
EPI_COLORS = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#6B5530"}

# Published whole-tissue cohort (config.cohort.whole_tissue; TMA + FTE excluded here).
XEN_ORDER = list(CFG["cohort"]["whole_tissue"])
XEN_COLORS = dict(zip(XEN_ORDER, [
    "#87CEFA", "#56AFC4", "#5665B6", "#8A5DAF",
    "#8FBC8F", "#2E8B57", "#D14E6C", "#B87A7A",
]))

# ---------- Font sizes ----------
FONT_TITLE = 8
FONT_TICK = 6
FONT_LABEL = 7
FONT_ANNO = 6

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FONT_TICK,
    "axes.titlesize":    FONT_TITLE,
    "axes.labelsize":    FONT_LABEL,
    "xtick.labelsize":   FONT_TICK,
    "ytick.labelsize":   FONT_TICK,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
    "axes.linewidth":    0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})


def main():
    t0 = time.time()
    print("=" * 72)
    print(f"  {STEM}")
    print("=" * 72)

    # ---- Atlas polarization (xenium-panel UCell, current NMF labels) ----
    atl = pd.read_csv(ATLAS_CSV, usecols=["SecA_UCell", "SecB_UCell", "schema_nmf"])
    # Standardize label "Transitioning" -> "Intermediate".
    atl["schema_nmf"] = atl["schema_nmf"].replace({"Transitioning": "Intermediate"})
    atl = atl[atl["schema_nmf"].isin(EPI_ORDER)].copy()
    atl["polarization"] = atl["SecB_UCell"] - atl["SecA_UCell"]
    atl["group"] = atl["schema_nmf"]
    print(f"   atlas : {len(atl):,} cells, {atl['group'].value_counts().to_dict()}")

    # ---- Xenium polarization, whole-tissue only (drop sfe_tma) ----
    xen = pd.read_csv(XEN_CSV, usecols=["sample", "SecA_UCell", "SecB_UCell", "polarization_UCell"])
    xen = xen[xen["sample"] != "sfe_tma"].copy()
    xen["group"] = xen["sample"].str.replace("sfe_", "", regex=False)
    reconstructed = xen["SecB_UCell"] - xen["SecA_UCell"]
    assert np.allclose(reconstructed, xen["polarization_UCell"], atol=1e-6), \
        "polarization_UCell column disagrees with SecB - SecA"
    xen["polarization"] = xen["polarization_UCell"]
    missing = set(xen["group"].unique()) - set(XEN_ORDER)
    if missing:
        raise RuntimeError(f"Unexpected xenium samples: {missing}")
    print(f"   xenium WT: {len(xen):,} cells across {xen['group'].nunique()} samples")

    # ---- p25 reference lines (polarization UCell) from atlas ----
    secb_p25 = atl.loc[atl["group"] == "SecB", "polarization"].quantile(0.25)
    inter_p25 = atl.loc[atl["group"] == "Intermediate", "polarization"].quantile(0.25)
    print(f"   Atlas polarization p25:  SecB={secb_p25:.3f}  Intermediate={inter_p25:.3f}")

    # ---- Assemble violins ----
    all_groups = EPI_ORDER + XEN_ORDER
    all_colors = {**EPI_COLORS, **XEN_COLORS}

    data = []
    for g in all_groups:
        if g in EPI_ORDER:
            vals = atl.loc[atl["group"] == g, "polarization"].values
        else:
            vals = xen.loc[xen["group"] == g, "polarization"].values
        data.append(vals)

    fig_w = 0.55 * len(all_groups) + 1.2
    fig, ax = plt.subplots(figsize=(fig_w, 3.6))

    positions = np.arange(1, len(all_groups) + 1)
    parts = ax.violinplot(data, positions=positions, widths=0.82,
                          showmeans=False, showmedians=False, showextrema=False)
    for body, g in zip(parts["bodies"], all_groups):
        body.set_facecolor(all_colors[g])
        body.set_edgecolor("#4d4d4d")
        body.set_linewidth(0.4)
        body.set_alpha(0.95)

    ax.axhline(secb_p25, linestyle="--", linewidth=0.5, color=EPI_COLORS["SecB"])
    ax.axhline(inter_p25, linestyle="--", linewidth=0.5, color=EPI_COLORS["Intermediate"])
    xright = len(all_groups) + 0.35
    ax.text(xright, secb_p25, " SecB p25", ha="left", va="center",
            fontsize=FONT_ANNO, color=EPI_COLORS["SecB"])
    ax.text(xright, inter_p25, " Int p25", ha="left", va="center",
            fontsize=FONT_ANNO, color=EPI_COLORS["Intermediate"])

    ax.axvline(len(EPI_ORDER) + 0.5, linestyle="-", linewidth=0.4, color="#b3b3b3")

    ax.set_xticks(positions)
    ax.set_xticklabels(all_groups, rotation=45, ha="right")
    ax.set_ylabel("Polarization UCell  (SecB − SecA)", fontsize=FONT_LABEL)
    ax.set_xlabel("")
    ax.set_ylim(-0.6, 0.6)
    ax.axhline(0, color="0.7", linewidth=0.3, zorder=0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="both", length=2, pad=1.5)

    ax.set_title("Atlas vs Xenium whole-tissue — polarization UCell",
                 fontsize=FONT_TITLE, pad=4, loc="left")

    ymax = ax.get_ylim()[1]
    ax.text(2, ymax + 0.02, "Atlas (current SecB NMF)", ha="center", va="bottom",
            fontsize=FONT_ANNO, color="0.25")
    ax.text(len(EPI_ORDER) + (len(XEN_ORDER) + 1) / 2, ymax + 0.02,
            "Xenium whole tissue (per sample)", ha="center", va="bottom",
            fontsize=FONT_ANNO, color="0.25")

    fig.subplots_adjust(left=0.08, right=0.92, top=0.88, bottom=0.22)

    fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight")
    fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight")
    plt.close(fig)

    print(f"   -> {OUT_PNG}")
    print(f"   -> {OUT_SVG}")
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
