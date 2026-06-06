#!/usr/bin/env python3
"""
SF7 — CNV validation: within-clone epitype composition
======================================================

Purpose
    Within-clone epitype composition across all samples with CopyKAT subclone
    inference. Each CNV subclone (>= 20 cells) is a stacked bar split into
    SecA / Intermediate / SecB %. Bars grouped by sample, sorted by sample SecB
    fraction descending; a CNV-verdict ribbon sits above each sample group.

INPUTS
    output_root/19_cnv/tables/within_clone_coexistence.csv
    output_root/19_cnv/tables/per_sample_verdict.csv
    output_root/19_cnv/tables/sample_manifest.csv
        (atlas step 19_cnv chain; annotation rev = NMF schema_nmf)

OUTPUTS
    output_root/figures/supplementary/SF7_atlas_cnv_validation.{svg,png}

MANUSCRIPT PANEL(S)
    SF7.

RUNTIME TIER
    fast (renders from pre-computed tables).
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

# ============================================================================
# PATHS (central config)
# ============================================================================

COEX_CSV = path("output_root", "19_cnv", "tables", "within_clone_coexistence.csv")
VERDICT_CSV = path("output_root", "19_cnv", "tables", "per_sample_verdict.csv")
MANIFEST_CSV = path("output_root", "19_cnv", "tables", "sample_manifest.csv")
OUT_SVG = path("output_root", "figures", "supplementary", "SF7_atlas_cnv_validation.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF7_atlas_cnv_validation.png")

for p in (COEX_CSV, VERDICT_CSV, MANIFEST_CSV):
    assert os.path.exists(p), f"Missing: {p}"

# ============================================================================
# STYLE
# ============================================================================

# Epitype display labels: "Intermediate" (was "Transitioning"). The CSV column
# names (frac_transitioning) are the upstream 19_cnv schema and stay as-is.
EPI_PALETTE = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55"}
EPI_ORDER = ["SecA", "Intermediate", "SecB"]
EPI_COLUMN = {"SecA": "frac_secA", "Intermediate": "frac_transitioning", "SecB": "frac_secB"}

VERDICT_PALETTE = {"monoclonal": "#D9D9D9", "clonally_driven": "#5C4A33", "mixed": "#BFA17A"}
VERDICT_DISPLAY = {"monoclonal": "Monoclonal", "clonally_driven": "Clonally driven", "mixed": "Mixed"}

FA, FK, FN = 6, 5.5, 5

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":        FA,
    "axes.labelsize":   FA,
    "axes.titlesize":   0,
    "xtick.labelsize":  FK,
    "ytick.labelsize":  FK,
    "legend.fontsize":  FN,
    "svg.fonttype":     "none",
    "pdf.fonttype":     42,
})

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading tables...", flush=True)
coex = pd.read_csv(COEX_CSV)
verdict = pd.read_csv(VERDICT_CSV)
manifest = pd.read_csv(MANIFEST_CSV)

coex["sample_id"] = coex["sample_id"].astype(int)
verdict["sample_id"] = verdict["sample_id"].astype(int)
manifest["sample_id"] = manifest["sample_id"].astype(int)
print(f"  {len(coex)} clones across {coex['sample_id'].nunique()} samples")

totals = coex[[EPI_COLUMN[c] for c in EPI_ORDER]].sum(axis=1)
assert np.allclose(totals.values, 1.0, atol=1e-6)

# ============================================================================
# SAMPLE ORDER — descending sample_secB_frac
# ============================================================================

sample_order = (
    coex.groupby("sample_id")["sample_secB_frac"].first()
    .sort_values(ascending=False).index.tolist()
)

clone_records = []
for sid in sample_order:
    sub = coex[coex["sample_id"] == sid].sort_values("frac_secB", ascending=False)
    for _, row in sub.iterrows():
        clone_records.append({
            "sample_id": sid,
            "subclone": row["subclone"],
            "n_cells": int(row["n_cells"]),
            "frac_secA": row["frac_secA"],
            "frac_transitioning": row["frac_transitioning"],
            "frac_secB": row["frac_secB"],
        })

n_bars = len(clone_records)
print(f"  Total bars (clones): {n_bars}")

# ============================================================================
# SPLIT SAMPLES INTO ROWS (letter page)
# ============================================================================

sample_groups = []
current_sid = None
current_group = []
for rec in clone_records:
    if rec["sample_id"] != current_sid:
        if current_group:
            sample_groups.append((current_sid, current_group))
        current_sid = rec["sample_id"]
        current_group = [rec]
    else:
        current_group.append(rec)
if current_group:
    sample_groups.append((current_sid, current_group))

BAR_W = 0.78
INTER_GAP = 0.5
PAGE_W_IN = 7.5
PAGE_H_IN = 10.0

# Balance rows by total bar count (clones), each sample kept intact (~63 bars/row).
NROWS = max(4, int(np.ceil(n_bars / 63)))
bars_per_row = int(np.ceil(n_bars / NROWS))

rows_of_samples = []
current_row = []
current_bars = 0
for sid, recs in sample_groups:
    cost = len(recs)
    if current_row and current_bars + cost > bars_per_row and current_bars >= bars_per_row * 0.5:
        rows_of_samples.append(current_row)
        current_row = [(sid, recs)]
        current_bars = cost
    else:
        current_row.append((sid, recs))
        current_bars += cost
if current_row:
    rows_of_samples.append(current_row)

NROWS = len(rows_of_samples)
print(f"  Layout: {NROWS} rows across letter page")

# ============================================================================
# FIGURE — multi-row letter page
# ============================================================================

print("Plotting...", flush=True)

fig, axes = plt.subplots(NROWS, 1, figsize=(PAGE_W_IN, PAGE_H_IN),
                         gridspec_kw={"hspace": 0.55})
if NROWS == 1:
    axes = [axes]

for ri, (ax, row_samples) in enumerate(zip(axes, rows_of_samples)):
    row_records = []
    for sid, recs in row_samples:
        row_records.extend(recs)

    x_pos = []
    sc = {}  # sample_centres for this row
    cur = 0.0
    prev = None
    for rec in row_records:
        sid = rec["sample_id"]
        if prev is not None and sid != prev:
            cur += INTER_GAP
            sc[prev] = (sc[prev][0], x_pos[-1])
        if sid not in sc:
            sc[sid] = (cur, cur)
        x_pos.append(cur)
        cur += 1.0
        prev = sid
    sc[prev] = (sc[prev][0], x_pos[-1])
    x_pos = np.array(x_pos)

    n_row_bars = len(row_records)
    heights = np.array([[rec[EPI_COLUMN[c]] * 100 for c in EPI_ORDER] for rec in row_records])
    bot = np.zeros(n_row_bars)
    for j, cls in enumerate(EPI_ORDER):
        vals = heights[:, j]
        ax.bar(x_pos, vals, width=BAR_W, bottom=bot, color=EPI_PALETTE[cls],
               edgecolor="white", linewidth=0.3)
        bot = bot + vals

    RIB_Y = 104
    RIB_H = 3.0
    for sid, (xl, xr) in sc.items():
        v = verdict.loc[verdict["sample_id"] == sid, "verdict"].iloc[0]
        ax.add_patch(plt.Rectangle((xl - BAR_W / 2, RIB_Y), (xr - xl) + BAR_W, RIB_H,
                                   facecolor=VERDICT_PALETTE[v], edgecolor="white",
                                   linewidth=0.3, clip_on=False))

    for sid, (xl, xr) in sc.items():
        centre = (xl + xr) / 2
        ax.text(centre, -5.0, f"S{sid}", ha="center", va="top",
                fontsize=FN - 1, color="#333333", rotation=90)

    ax.set_xlim(x_pos.min() - 0.6, x_pos.max() + 0.6)
    ax.set_ylim(-0.5, 108)
    ax.set_xticks([])
    ax.set_yticks([0, 50, 100])
    ax.tick_params(axis="y", which="both", length=2, width=0.4, pad=1, labelsize=FK)
    if ri == NROWS // 2:
        ax.set_ylabel("Within-clone epitype composition (%)", fontsize=FA, labelpad=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_linewidth(0.5)

# Shared legend on the last row
ax_last = axes[-1]
epi_handles = [Patch(facecolor=EPI_PALETTE[c], edgecolor="white", linewidth=0.4, label=c)
               for c in EPI_ORDER]
leg1 = ax_last.legend(handles=epi_handles, title="Epitype", loc="upper left",
                      bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=FN, title_fontsize=FN,
                      handlelength=1.0, handleheight=0.9, handletextpad=0.4,
                      labelspacing=0.35, borderaxespad=0.0)
ax_last.add_artist(leg1)

ver_handles = [Patch(facecolor=VERDICT_PALETTE[v], edgecolor="white", linewidth=0.4,
                     label=VERDICT_DISPLAY[v]) for v in ["monoclonal", "clonally_driven", "mixed"]]
ax_last.legend(handles=ver_handles, title="CNV verdict", loc="lower left",
               bbox_to_anchor=(1.01, 0.0), frameon=False, fontsize=FN, title_fontsize=FN,
               handlelength=1.0, handleheight=0.9, handletextpad=0.4,
               labelspacing=0.35, borderaxespad=0.0)

n_mono = int((verdict["verdict"] == "monoclonal").sum())
n_driv = int((verdict["verdict"] == "clonally_driven").sum())
n_mix = int((verdict["verdict"] == "mixed").sum())
n_tot = len(verdict)
axes[0].text(1.01, 0.95,
             f"n = {n_tot} samples, {n_bars} clones\n"
             f"Monoclonal: {n_mono} ({n_mono/n_tot*100:.0f}%)\n"
             f"Clonally driven: {n_driv} ({n_driv/n_tot*100:.0f}%)\n"
             f"Mixed: {n_mix} ({n_mix/n_tot*100:.0f}%)",
             transform=axes[0].transAxes, fontsize=FN, va="top", ha="left", linespacing=1.3)

fig.savefig(OUT_SVG, format="svg", bbox_inches="tight", pad_inches=0.08)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.08)
plt.close(fig)

print(f"  Saved: {OUT_SVG}")
print(f"  Saved: {OUT_PNG}")
