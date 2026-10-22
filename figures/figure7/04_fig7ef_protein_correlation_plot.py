#!/usr/bin/env python3
# ============================================================================
# Figure 7e,7f — Transcript vs protein and polarization vs protein scatter row
# ----------------------------------------------------------------------------
# PURPOSE
#   6-panel row across letter-page width:
#     Panel 1: Xenium transcript vs IF protein MFI (KRT7 + KRT19 combined) (7e)
#     Panels 2-6: Polarization score vs each protein marker
#                 (KRT7, KRT19, KRT18, VIM, E-cadherin) (7f)
#   Per-core Pearson r + Spearman ρ.
#
# INPUTS
#   $PROTEIN_CSV (this is the default): output_root/41_xenium_protein_correlation/per_core_xenium_protein_final.csv from
#     03_fig7ef_protein_correlation.R (final TMA A cores, SingleR annotation).
#
# OUTPUTS
#   figures_dir/figure7/xenium_protein_correlation_row.{svg,png}
#
# MANUSCRIPT PANEL(S): Fig 7e (transcript vs protein), Fig 7f (polarization vs protein)
# RUNTIME TIER: fast (reads one per-core CSV)
# ============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# --- Config (script is 2 levels deep: figures/figure7/) ---------------------
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED

np.random.seed(SEED)

DATA_CSV = os.environ.get("PROTEIN_CSV", str(path(   # default: where 03_fig7ef_protein_correlation.R writes it
    "output_root", "41_xenium_protein_correlation", "per_core_xenium_protein_final.csv")))
FIG_DIR  = path("figures_dir", "figure7")

OUT_STEM = "xenium_protein_correlation_row"
OUT_SVG  = os.path.join(FIG_DIR, f"{OUT_STEM}.svg")
OUT_PNG  = os.path.join(FIG_DIR, f"{OUT_STEM}.png")

# ============================================================================
# STYLE
# ============================================================================
# Figure 7 house typography: 8pt titles/labels, 6pt ticks, 7pt in-panel r/rho.
FA, FK, FN = 8, 6, 7

plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FA,
    "axes.labelsize":    FA,
    "axes.titlesize":    FA,
    "xtick.labelsize":   FK,
    "ytick.labelsize":   FK,
    "legend.fontsize":   FN,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

MARKER_PALETTE = {
    "KRT7":        "#5665B6",
    "KRT19":       "#D14E6C",
    "KRT18":       "#8A5DAF",
    "VIM":         "#2E8B57",
    "E-cadherin":  "#E6A141",
}

DOT_COLOR = "#555555"
LINE_COLOR = "#D14E6C"

# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading per-core data...", flush=True)
df = pd.read_csv(DATA_CSV)
print(f"  Cores: {len(df)}")

# ============================================================================
# PLOT — 1×6 across letter width (8.5 in)
# ============================================================================
print("Plotting...", flush=True)

fig, axes = plt.subplots(1, 6, figsize=(8.5, 2.1))

# --- Panel 1 (Fig 7e): Transcript vs protein (KRT7 + KRT19 combined) --------
ax = axes[0]

for xcol, ycol, label in [("KRT7_xenium", "KRT7", "KRT7"),
                           ("KRT19_xenium", "KRT19", "KRT19")]:
    mask = df[[xcol, ycol]].notna().all(axis=1)
    x = df.loc[mask, xcol].values
    y = df.loc[mask, ycol].values
    color = MARKER_PALETTE[label]

    ax.scatter(x, y, s=5, color=color, alpha=0.45,
               edgecolors="none", rasterized=False, zorder=2)

    slope, intercept, r, p, se = stats.linregress(x, y)
    xfit = np.linspace(x.min(), x.max(), 100)
    rho = stats.spearmanr(x, y)[0]
    ax.plot(xfit, slope * xfit + intercept,
            color=color, linewidth=1.0, zorder=3,
            label=f"{label}  r={r:.2f}")

ax.set_xlabel("Xenium transcript\n(mean logcounts)", fontsize=FA)
ax.set_ylabel("Protein MFI", fontsize=FA)
ax.legend(fontsize=FN - 0.5, frameon=False, loc="upper left",
          handlelength=1.2, handletextpad=0.3, labelspacing=0.3)

# --- Panels 2–6 (Fig 7f): Polarization vs each protein marker ---------------
POLAR_MARKERS = [
    ("KRT7",  "KRT7"),
    ("KRT19", "KRT19"),
    ("KRT18", "KRT18"),
    ("VIM",   "VIM"),
    ("ECAD",  "E-cadherin"),
]

for ax, (ycol, label) in zip(axes[1:], POLAR_MARKERS):
    mask = df[["polarization_mean", ycol]].notna().all(axis=1)
    x = df.loc[mask, "polarization_mean"].values
    y = df.loc[mask, ycol].values
    color = MARKER_PALETTE[label]

    ax.scatter(x, y, s=5, color=color, alpha=0.45,
               edgecolors="none", rasterized=False, zorder=2)

    slope, intercept, r, p, se = stats.linregress(x, y)
    xfit = np.linspace(x.min(), x.max(), 100)
    rho = stats.spearmanr(x, y)[0]
    ax.plot(xfit, slope * xfit + intercept,
            color=color, linewidth=1.0, zorder=3)

    ax.text(0.03, 0.97,
            f"r = {r:.2f}\nρ = {rho:.2f}",
            transform=ax.transAxes, fontsize=FN, va="top", ha="left",
            color=color, linespacing=1.2)

    # Marker identity lives on the y-axis ("KRT7 MFI" etc.); the redundant bold
    # per-panel title and the "(mean UCell)" x-label second line are dropped.
    ax.set_xlabel("Polarization", fontsize=FA)
    ax.set_ylabel(f"{label} MFI", fontsize=FA)

# --- Shared formatting -------------------------------------------------------
for ax in axes:
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=2)

fig.tight_layout(w_pad=0.8)
fig.savefig(OUT_SVG, format="svg", bbox_inches="tight", pad_inches=0.04)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.04)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
