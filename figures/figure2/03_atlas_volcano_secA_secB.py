#!/usr/bin/env python3
"""
Figure 2C — Volcano plot, SecA vs SecB (wide horizontal, signature markers labelled)
====================================================================================
PURPOSE
    Wide horizontal volcano of the SecA-vs-SecB Wilcoxon DE. Top-10 DEGs per side
    highlighted; SecA/SecB signature markers drawn large/outlined and labelled
    with adjustText repulsion.

INPUTS
    - panel_i_deg_results.csv :
        output_root/fig_secretory_polarization/data/panel_i_deg_results.csv
        (Wilcoxon SecA vs SecB on schema_nmf; produced by
         figures/_prep/fig_secretory_polarization_00_prepare_data.py)
    - SecA/SecB 7-gene signatures from shared/signatures.yml.

OUTPUTS
    - figures_dir/atlas_volcano_secA_secB.{svg,png}

MANUSCRIPT PANEL(S): Fig 2C.

RUNTIME TIER: fast.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import yaml
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
DEG_CSV = path("data_root", "2026_final_atlas", "output", "fig_secretory_polarization", "data", "panel_i_deg_results.csv")
SIG_YML = Path(__file__).resolve().parents[2] / "shared" / "signatures.yml"
OUT_SVG = path("figures_dir", "atlas_volcano_secA_secB.svg")
OUT_PNG = path("figures_dir", "atlas_volcano_secA_secB.png")

# ---------- Signatures (shared source of truth) ----------
SIG = yaml.safe_load(open(SIG_YML))
SECA_MARKERS = SIG["SecA"]
SECB_MARKERS = SIG["SecB"]
ALL_MARKERS = set(SECA_MARKERS + SECB_MARKERS)

# ---------- Style ----------
FA, FK, FN = 7, 6.5, 6
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         FA,
    "axes.labelsize":    FA,
    "xtick.labelsize":   FK,
    "ytick.labelsize":   FK,
    "legend.fontsize":   FN,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

SECA_COLOR = "#E6A141"
SECB_COLOR = "#9A7D55"
NS_COLOR = "#D9D9D9"
N_TOP = 10

# ---------- Load + prepare ----------
df = pd.read_csv(DEG_CSV)
df["pval_adj"] = df["pval_adj"].clip(lower=1e-300)
df["neg_log10p"] = -np.log10(df["pval_adj"])
YMAX = 80
df["neg_log10p_plot"] = df["neg_log10p"].clip(upper=YMAX)
XMIN, XMAX = -5, 6
df["log2fc_plot"] = df["log2fc"].clip(lower=XMIN + 0.1, upper=XMAX - 0.1)
LFC_THRESH, PADJ_THRESH = 0.5, 0.05
df["category"] = "ns"
df.loc[(df["log2fc"] < -LFC_THRESH) & (df["pval_adj"] < PADJ_THRESH), "category"] = "SecA"
df.loc[(df["log2fc"] > LFC_THRESH) & (df["pval_adj"] < PADJ_THRESH), "category"] = "SecB"
df["is_marker"] = df["gene"].isin(ALL_MARKERS)

seca_sig = df[(df["category"] == "SecA") & ~df["is_marker"]].copy()
secb_sig = df[(df["category"] == "SecB") & ~df["is_marker"]].copy()
top10_seca = seca_sig.nlargest(N_TOP, "neg_log10p")["gene"].tolist()
top10_secb = secb_sig.nlargest(N_TOP, "neg_log10p")["gene"].tolist()
top10_all = set(top10_seca + top10_secb)

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(6.0, 2.0))
ns = df[df["category"] == "ns"]
ax.scatter(ns["log2fc_plot"], ns["neg_log10p_plot"], s=1.5, c=NS_COLOR, alpha=0.3,
           linewidths=0, rasterized=True, zorder=1)
seca_bg = df[(df["category"] == "SecA") & ~df["is_marker"] & ~df["gene"].isin(top10_all)]
ax.scatter(seca_bg["log2fc_plot"], seca_bg["neg_log10p_plot"], s=3, c=SECA_COLOR,
           alpha=0.4, linewidths=0, rasterized=True, zorder=2)
secb_bg = df[(df["category"] == "SecB") & ~df["is_marker"] & ~df["gene"].isin(top10_all)]
ax.scatter(secb_bg["log2fc_plot"], secb_bg["neg_log10p_plot"], s=3, c=SECB_COLOR,
           alpha=0.4, linewidths=0, rasterized=True, zorder=2)
t10a = df[df["gene"].isin(top10_seca)]
ax.scatter(t10a["log2fc_plot"], t10a["neg_log10p_plot"], s=8, c=SECA_COLOR,
           alpha=0.6, linewidths=0, rasterized=True, zorder=3)
t10b = df[df["gene"].isin(top10_secb)]
ax.scatter(t10b["log2fc_plot"], t10b["neg_log10p_plot"], s=8, c=SECB_COLOR,
           alpha=0.6, linewidths=0, rasterized=True, zorder=3)
seca_m = df[df["gene"].isin(SECA_MARKERS)]
ax.scatter(seca_m["log2fc_plot"], seca_m["neg_log10p_plot"], s=22, c=SECA_COLOR,
           edgecolors="black", linewidths=0.4, zorder=5)
secb_m = df[df["gene"].isin(SECB_MARKERS)]
ax.scatter(secb_m["log2fc_plot"], secb_m["neg_log10p_plot"], s=22, c=SECB_COLOR,
           edgecolors="black", linewidths=0.4, zorder=5)

ax.axhline(-np.log10(PADJ_THRESH), color="grey", linewidth=0.4, linestyle="--", zorder=0)
ax.axvline(-LFC_THRESH, color="grey", linewidth=0.4, linestyle="--", zorder=0)
ax.axvline(LFC_THRESH, color="grey", linewidth=0.4, linestyle="--", zorder=0)
ax.axhline(YMAX, color="grey", linewidth=0.3, linestyle=":", alpha=0.4, zorder=0)

ax.set_xlim(XMIN, XMAX)
ax.set_ylim(-3, YMAX + 5)
LABEL_SIZE = 5.5
texts = []
for gene in SECA_MARKERS:
    row = df[df["gene"] == gene]
    if len(row) == 0:
        continue
    row = row.iloc[0]
    texts.append(ax.text(row["log2fc_plot"], row["neg_log10p_plot"], gene,
                         fontsize=LABEL_SIZE + 0.5, fontweight="bold", fontstyle="italic",
                         color="#8B6914", ha="center", va="center", zorder=7))
for gene in SECB_MARKERS:
    row = df[df["gene"] == gene]
    if len(row) == 0:
        continue
    row = row.iloc[0]
    texts.append(ax.text(row["log2fc_plot"], row["neg_log10p_plot"], gene,
                         fontsize=LABEL_SIZE + 0.5, fontweight="bold", fontstyle="italic",
                         color="#5C4A2A", ha="center", va="center", zorder=7))
adjust_text(texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="#CCCCCC", linewidth=0.3, shrinkA=0, shrinkB=3),
            expand=(2.0, 2.0), force_text=(2.0, 2.0), force_points=(1.0, 1.0),
            only_move="xy", lim=3000, ensure_inside_axes=True)

ax.set_xlabel("log$_2$ fold change (SecA ← → SecB)", fontsize=FA)
ax.set_ylabel("$-$log$_{10}$ adj. p-value", fontsize=FA)
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)
ax.tick_params(width=0.5, length=2)

legend_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=SECA_COLOR,
           markeredgecolor="black", markeredgewidth=0.4, markersize=4, label="SecA markers"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor=SECB_COLOR,
           markeredgecolor="black", markeredgewidth=0.4, markersize=4, label="SecB markers"),
]
ax.legend(handles=legend_handles, loc="upper right", frameon=False, fontsize=FN, handletextpad=0.3)

fig.tight_layout()
fig.savefig(OUT_SVG, format="svg")
fig.savefig(OUT_PNG, format="png")
plt.close(fig)
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}")
