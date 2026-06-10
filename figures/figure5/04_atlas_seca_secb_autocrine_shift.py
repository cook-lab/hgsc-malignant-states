#!/usr/bin/env python3
# ============================================================================
# Figure 5F — SecA -> SecB autocrine L-R shift dumbbell
# ----------------------------------------------------------------------------
# PURPOSE
#   Dumbbell plot of autocrine ligand-receptor signalling categories. Each of
#   7 functional categories is a horizontal dumbbell from the SecA autocrine
#   pair count (gold) to the SecB count (brown), ordered by Δ (SecB − SecA).
#   Autocrine = self-loop pairs (source == target == pole) with lrscore > 0.5;
#   categorisation follows fig2k_autocrine_budget.py (2026-04-11).
#
# INPUTS
#   data_root/2026_final_atlas/output/17_cellcomm_nmf/tables/17b_liana_global.csv
#     (LIANA NMF; Fig 5F / Supp Data 7 upstream)
#
# OUTPUTS
#   figures_dir/figure5/atlas_seca_secb_autocrine_shift.{png,svg}
#
# MANUSCRIPT PANEL(S): Fig 5F
# RUNTIME TIER: fast (reads one CSV)
#
# NOTE: Fig 5F LIANA autocrine is flagged "under investigation";
#       code migrated faithfully, no logic change.
# ============================================================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Config (script is 2 levels deep: figures/figure5/) ---------------------
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED

np.random.seed(SEED)

LIANA = path("data_root", "2026_final_atlas", "output",
             "17_cellcomm_nmf", "tables", "17b_liana_global.csv")
FIG_DIR = path("figures_dir", "figure5")
OUT_PNG = os.path.join(FIG_DIR, "atlas_seca_secb_autocrine_shift.png")
OUT_SVG = os.path.join(FIG_DIR, "atlas_seca_secb_autocrine_shift.svg")

# ── Style ─────────────────────────────────────────────────────────────────
FA, FK, FN = 8, 7, 6.5
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       FA,
    "axes.labelsize":  FA,
    "axes.titlesize":  0,
    "xtick.labelsize": FK,
    "ytick.labelsize": FK,
    "legend.fontsize": FN,
    "pdf.fonttype":    42,
    "svg.fonttype":    "none",
    "savefig.dpi":     450,
    "figure.dpi":      150,
})

SECA_COL  = "#E6A141"
SECB_COL  = "#6B4D2E"
GAIN_COL  = "#7A2B35"   # burgundy — SecB > SecA
LOSS_COL  = "#3D6B85"   # cool slate — SecB < SecA
NULL_COL  = "#999999"
SIG_THR   = 0.5

CATEGORISE = {
    "ECM / Plasminogen":   ["PLAU","SERPIN","FN1","COL","POSTN","LAMA","LAMB",
                            "PTN","MXRA","AGRN"],
    "Integrin / Receptor": ["ITGB","ITGA","IGF2R","LRP","DAG1"],
    "DAMP / TLR":          ["HMGB","S100A","ANXA","DCN → TLR","TLR","HSPA"],
    "Chemokine":           ["CCL","CCR","CXCL","CXCR"],
    "Activin / BMP":       ["INHBA","ACTR","GDF","BMPR"],
    "Lipid / GPCR":        ["PSAP","SORT","GNAI","S1PR","APOE"],
    "Semaphorin":          ["SEMA","PLXN","NRP","FARP"],
}

def categorise(lr):
    s = lr.upper()
    for cat, pats in CATEGORISE.items():
        if any(p in s for p in pats):
            return cat
    return "Other"

# ── Load + tally ──────────────────────────────────────────────────────────
print(f"Loading: {LIANA}")
li = pd.read_csv(LIANA)

def autocrine_pairs(pole):
    sub = li[(li["source"] == pole) & (li["target"] == pole)].copy()
    sub["lr"] = sub["ligand_complex"] + " → " + sub["receptor_complex"]
    return sub.set_index("lr")["lrscore"]

sA = autocrine_pairs("SecA epithelium")
sB = autocrine_pairs("SecB epithelium")

def category_counts(scores, thr=SIG_THR):
    sig = scores[scores > thr]
    cats = sig.index.map(categorise)
    return pd.Series(cats).value_counts()

cA = category_counts(sA)
cB = category_counts(sB)

CATS = list(CATEGORISE.keys())
df = pd.DataFrame({"SecA": [int(cA.get(c, 0)) for c in CATS],
                   "SecB": [int(cB.get(c, 0)) for c in CATS]},
                  index=CATS)
df["delta"] = df["SecB"] - df["SecA"]
df = df.sort_values("delta", ascending=False)
print(df)

# Pull SecB-exclusive HMGB1 / S100A → TLR receptor pairs (narrative anchors)
sA_sig = sA[sA > SIG_THR]
key_pairs_idx = [
    p for p in sB.index
    if p.startswith(("HMGB1 ", "S100A8 ", "S100A9 "))
    and "TLR" in p
    and sB[p] > SIG_THR
    and p not in sA_sig.index
]
key_pairs = sB.loc[key_pairs_idx].sort_values(ascending=False)
print(f"\nSecB-exclusive DAMP→TLR pairs ({len(key_pairs)}):")
for k, v in key_pairs.items():
    print(f"  {k:<28s}  lrscore={v:.3f}")

# ── Plot ──────────────────────────────────────────────────────────────────
fig_w_mm, fig_h_mm = 65, 35
fig, ax = plt.subplots(figsize=(fig_w_mm / 25.4, fig_h_mm / 25.4))

cats   = df.index.tolist()
y_pos  = np.arange(len(cats))[::-1]
xmax   = max(df["SecA"].max(), df["SecB"].max()) * 1.12

for x in (0, 20, 40, 60):
    ax.axvline(x, color="#E0E0E0", lw=0.3, zorder=0)

for y, cat in zip(y_pos, cats):
    a = int(df.at[cat, "SecA"])
    b = int(df.at[cat, "SecB"])
    delta = b - a
    if delta > 0:    shaft = GAIN_COL
    elif delta < 0:  shaft = LOSS_COL
    else:            shaft = NULL_COL
    lw = 0.3 + min(abs(delta), 16) * 0.08

    if a != b:
        x0, x1 = (a, b) if a <= b else (b, a)
        ax.plot([x0, x1], [y, y], color=shaft, lw=lw, alpha=0.85,
                solid_capstyle="round", zorder=2)

    if delta != 0:
        ax.annotate("", xy=(b, y), xytext=(a + (b - a) * 0.85, y),
                    arrowprops=dict(arrowstyle="->", color=shaft,
                                    lw=lw, shrinkA=0, shrinkB=0),
                    zorder=3)

    ax.scatter(a, y, s=18, color=SECA_COL, edgecolor="white",
               linewidth=0.4, zorder=4)
    ax.scatter(b, y, s=18, color=SECB_COL, edgecolor="white",
               linewidth=0.4, zorder=4)

    pad = xmax * 0.02
    if a < b:
        ax.text(a - pad, y, str(a), ha="right", va="center",
                fontsize=FN, color=SECA_COL, fontweight="bold")
        ax.text(b + pad, y, str(b), ha="left", va="center",
                fontsize=FN, color=SECB_COL, fontweight="bold")
    elif a > b:
        ax.text(a + pad, y, str(a), ha="left", va="center",
                fontsize=FN, color=SECA_COL, fontweight="bold")
        ax.text(b - pad, y, str(b), ha="right", va="center",
                fontsize=FN, color=SECB_COL, fontweight="bold")
    else:
        ax.text(a + pad, y + 0.18, str(a), ha="left", va="bottom",
                fontsize=FN, color=SECA_COL, fontweight="bold")
        ax.text(b + pad, y - 0.18, str(b), ha="left", va="top",
                fontsize=FN, color=SECB_COL, fontweight="bold")

    sign = "+" if delta > 0 else ("" if delta < 0 else "±")
    ax.text(xmax * 1.05, y, f"{sign}{delta}",
            ha="left", va="center", fontsize=FN,
            color=shaft, fontweight="bold")

ax.set_yticks(y_pos)
ax.set_yticklabels(cats)
ax.set_xlim(-xmax * 0.08, xmax * 1.16)
ax.set_ylim(-0.6, len(cats) - 0.4)
ax.set_xlabel("Autocrine L–R pairs (lrscore > 0.5)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_linewidth(0.4); ax.spines[s].set_color("black")
ax.tick_params(width=0.4, length=2, colors="black")

# SecB-exclusive DAMP/TLR pairs (HMGB1→TLR4/2, S100A8/9→TLR4) described in caption.

fig.tight_layout(pad=0.2)
fig.savefig(OUT_PNG, bbox_inches="tight", dpi=450)
fig.savefig(OUT_SVG, bbox_inches="tight")
plt.close(fig)
print(f"\nSaved:\n  {os.path.basename(OUT_PNG)}\n  {os.path.basename(OUT_SVG)}")
