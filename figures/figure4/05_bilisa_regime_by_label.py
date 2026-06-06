#!/usr/bin/env python3
"""
Figure 4F — BiLISA regime composition by discrete secretory label
=================================================================
PURPOSE
    Horizontal stacked bars: among significant cells, the fraction of each
    epitype (SecA / Intermediate / SecB) in each BiLISA regime (HH / HL / LH /
    LL) on the SecA_UCell vs neighbour-SecB_UCell axes. Key result: most SecB
    cells reside in SecB-dominant (LH) neighbourhoods.

INPUTS
    - output_root/44_spatial_autocorrelation/bilisa_vs_label_crosstab.csv
      (columns: cell_label, regime_sig, pct, total)

OUTPUTS
    - figures_dir/bilisa_regime_by_label.{png,svg}

MANUSCRIPT PANEL(S): Fig 4F.

RUNTIME TIER: fast.

NOTE: epitype label standardized "Transitioning" -> "Intermediate".
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
CROSSTAB = path("output_root", "44_spatial_autocorrelation", "bilisa_vs_label_crosstab.csv")
OUT_PNG = path("figures_dir", "bilisa_regime_by_label.png")
OUT_SVG = path("figures_dir", "bilisa_regime_by_label.svg")

# ---------- Style ----------
FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       FA,
    "axes.labelsize":  FA,
    "axes.titlesize":  FA,
    "xtick.labelsize": FK,
    "ytick.labelsize": FK,
    "legend.fontsize": FN,
    "pdf.fonttype":    42,
    "svg.fonttype":    "none",
    "savefig.dpi":     450,
    "figure.dpi":      150,
})

REGIME_PAL = {"HL": "#E6A141", "HH": "#C08E48", "LH": "#9A7D55", "LL": "#B0B0B0"}
REGIME_ORDER = ["HL", "HH", "LH", "LL"]
REGIME_LABELS = {"HL": "SecA neighborhood", "HH": "Co-localized",
                 "LH": "SecB neighborhood", "LL": "Both low"}

CELL_ORDER = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
CELL_LABELS = {"SecA epithelium": "SecA", "Intermediate epithelium": "Intermediate",
               "SecB epithelium": "SecB"}

# ---------- Load ----------
ct = pd.read_csv(CROSSTAB)
piv = ct.pivot(index="cell_label", columns="regime_sig", values="pct").fillna(0)
piv = piv.reindex(index=CELL_ORDER, columns=REGIME_ORDER).fillna(0)

for cl in CELL_ORDER:
    row = piv.loc[cl]
    n = ct.loc[ct["cell_label"] == cl, "total"].values[0]
    print(f"  {CELL_LABELS[cl]:15s} (n={n:,}): "
          + "  ".join(f"{r}={row[r]:.1f}%" for r in REGIME_ORDER))

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(4.25, 1.65))
y_pos = np.arange(len(CELL_ORDER))
bar_h = 0.55
lefts = np.zeros(len(CELL_ORDER))
for regime in REGIME_ORDER:
    widths = piv[regime].values
    ax.barh(y_pos, widths, height=bar_h, left=lefts, color=REGIME_PAL[regime],
            edgecolor="white", linewidth=0.3, label=REGIME_LABELS[regime])
    for i, (w, l) in enumerate(zip(widths, lefts)):
        if w >= 8:
            ax.text(l + w / 2, y_pos[i], f"{w:.0f}%", ha="center", va="center",
                    fontsize=FN, color="white")
    lefts += widths

ax.set_yticks(y_pos)
ax.set_yticklabels([CELL_LABELS[c] for c in CELL_ORDER])
ax.invert_yaxis()
ax.set_xlim(0, 100)
ax.set_xlabel("% of significant cells")
ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, frameon=False,
          fontsize=FN, handlelength=1.0, handleheight=0.7, labelspacing=0.4)

for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(left=False, bottom=True, width=0.4, length=2)
ax.spines["bottom"].set_visible(True)
ax.spines["bottom"].set_linewidth(0.4)

fig.tight_layout(pad=0.3)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {os.path.basename(OUT_PNG)}\n  Saved: {os.path.basename(OUT_SVG)}\nDone.")
