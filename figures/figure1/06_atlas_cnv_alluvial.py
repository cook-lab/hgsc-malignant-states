#!/usr/bin/env python3
"""
Figure 1J — CNV subclone -> epitype alluvial (3 patients)
=========================================================
PURPOSE
    Single alluvial: 3 patients' CopyKAT CNV subclones on the left fan into the
    shared SecA / Intermediate / SecB epitype nodes on the right. Each clone
    flows into all three epitypes — polarization is transcriptional, not clonal.
    Node/column labels added in Illustrator at layout time.

INPUTS
    - 19_cnv within_clone_coexistence.csv :
        output_root/05_cnv/tables/within_clone_coexistence.csv
        (columns: sample_id, n_cells, n_secA, n_transitioning, n_secB)

OUTPUTS
    - figures_dir/atlas_cnv_alluvial.{svg,png}

MANUSCRIPT PANEL(S): Fig 1J.

RUNTIME TIER: fast.

NOTE: epitype node label standardized "Transitioning" -> "Intermediate". The
upstream CNV cache column name (n_transitioning) is read as-is.
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
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, FancyBboxPatch
from matplotlib.colors import to_rgb, to_hex

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ---------- Paths ----------
COEX_CSV = path("data_root", "2026_final_atlas", "output", "19_cnv", "tables", "within_clone_coexistence.csv")
OUT_SVG = path("figures_dir", "atlas_cnv_alluvial.svg")
OUT_PNG = path("figures_dir", "atlas_cnv_alluvial.png")

assert os.path.exists(COEX_CSV), f"Missing: {COEX_CSV}"

# ---------- Example samples ----------
SAMPLE_IDS = [16, 100, 62]
SAMPLE_NAMES = ["Patient A\n(monoclonal)", "Patient B", "Patient C"]

# ---------- Style ----------
FA, FK, FN = 10, 9, 8
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

EPITYPE_COLORS = {"SecA": "#E6A141", "Intermediate": "#CF8C2E", "SecB": "#B8741A"}
EPITYPE_ORDER = ["SecA", "Intermediate", "SecB"]
PATIENT_BASE = ["#8C8C8C", "#4A7FA5", "#4A967A"]


def shade_variants(base_hex, n, spread=0.28):
    r, g, b = to_rgb(base_hex)
    if n == 1:
        return [base_hex]
    shades = []
    for i in range(n):
        t = (i / (n - 1)) - 0.5
        factor = 1.0 + t * spread * 2
        shades.append(to_hex((max(0, min(1, r * factor)),
                              max(0, min(1, g * factor)),
                              max(0, min(1, b * factor)))))
    return shades


# ---------- Load data ----------
coex = pd.read_csv(COEX_CSV)
clone_records = []
grand_total = 0
for pi, sid in enumerate(SAMPLE_IDS):
    sub = coex[coex["sample_id"] == sid].sort_values("n_cells", ascending=False).reset_index(drop=True)
    colors = shade_variants(PATIENT_BASE[pi], len(sub))
    for ci, (_, row) in enumerate(sub.iterrows()):
        clone_records.append({
            "patient_idx": pi, "patient_name": SAMPLE_NAMES[pi], "clone_idx": ci,
            "clone_label": f"C{ci+1}", "n_cells": int(row["n_cells"]),
            "n_secA": int(row["n_secA"]), "n_trans": int(row["n_transitioning"]),
            "n_secB": int(row["n_secB"]), "color": colors[ci],
        })
        grand_total += int(row["n_cells"])


def draw_flow(ax, y_src_top, h_src, y_dst_top, h_dst, x0, x1, color, alpha=0.35):
    xm_l = x0 + 0.55 * (x1 - x0)
    xm_r = x0 + 0.45 * (x1 - x0)
    y_src_bot = y_src_top - h_src
    y_dst_bot = y_dst_top - h_dst
    verts = [(x0, y_src_top), (xm_l, y_src_top), (xm_r, y_dst_top), (x1, y_dst_top),
             (x1, y_dst_bot), (xm_r, y_dst_bot), (xm_l, y_src_bot), (x0, y_src_bot),
             (x0, y_src_top)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.LINETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="white",
                           alpha=alpha, linewidth=0.15, zorder=2))


def draw_block(ax, x_center, y_top, height, width, color, edgecolor="white", linewidth=0.4):
    ax.add_patch(FancyBboxPatch((x_center - width / 2, y_top - height), width, height,
                                boxstyle="round,pad=0.006", facecolor=color,
                                edgecolor=edgecolor, linewidth=linewidth, zorder=3))


# ---------- Layout ----------
INTER_PATIENT_GAP = 0.045
INTRA_CLONE_GAP = 0.012
n_patients = len(SAMPLE_IDS)
n_clones_per = [sum(1 for c in clone_records if c["patient_idx"] == pi) for pi in range(n_patients)]
total_intra_gaps = sum(max(0, nc - 1) for nc in n_clones_per) * INTRA_CLONE_GAP
total_inter_gaps = (n_patients - 1) * INTER_PATIENT_GAP
available_h = 1.0 - total_intra_gaps - total_inter_gaps

clone_heights = [cr["n_cells"] / grand_total * available_h for cr in clone_records]

clone_y_tops = []
y_cursor = 1.0
prev_pi = 0
for idx, cr in enumerate(clone_records):
    if cr["patient_idx"] != prev_pi:
        y_cursor -= INTER_PATIENT_GAP
        prev_pi = cr["patient_idx"]
    elif idx > 0:
        y_cursor -= INTRA_CLONE_GAP
    clone_y_tops.append(y_cursor)
    y_cursor -= clone_heights[idx]

EP_GAP = 0.035
ep_totals = {ep: 0 for ep in EPITYPE_ORDER}
for cr in clone_records:
    ep_totals["SecA"] += cr["n_secA"]
    ep_totals["Intermediate"] += cr["n_trans"]
    ep_totals["SecB"] += cr["n_secB"]

available_h_ep = 1.0 - EP_GAP * (len(EPITYPE_ORDER) - 1)
ep_heights = {ep: ep_totals[ep] / grand_total * available_h_ep for ep in EPITYPE_ORDER}

ep_y_tops = {}
y_cursor = 1.0
for ep in EPITYPE_ORDER:
    ep_y_tops[ep] = y_cursor
    y_cursor -= ep_heights[ep] + EP_GAP

# ---------- Plot ----------
fig, ax = plt.subplots(figsize=(4.5, 2.8))
BLOCK_W = 0.05
X_LEFT = 0.12
X_RIGHT = 0.88

clone_cursors = [clone_y_tops[i] for i in range(len(clone_records))]
ep_cursors = {ep: ep_y_tops[ep] for ep in EPITYPE_ORDER}

for idx, cr in enumerate(clone_records):
    for ep, n_key in [("SecA", "n_secA"), ("Intermediate", "n_trans"), ("SecB", "n_secB")]:
        n_flow = cr[n_key]
        if n_flow == 0:
            continue
        frac = n_flow / grand_total
        h_src = frac * available_h
        h_dst = frac * available_h_ep
        draw_flow(ax, clone_cursors[idx], h_src, ep_cursors[ep], h_dst,
                  X_LEFT + BLOCK_W / 2, X_RIGHT - BLOCK_W / 2, color=cr["color"], alpha=0.40)
        clone_cursors[idx] -= h_src
        ep_cursors[ep] -= h_dst

patient_label_positions = {}
for idx, cr in enumerate(clone_records):
    h = clone_heights[idx]
    y_top = clone_y_tops[idx]
    draw_block(ax, X_LEFT, y_top, h, BLOCK_W, cr["color"])
    ax.text(X_LEFT - BLOCK_W / 2 - 0.006, y_top - h / 2, cr["clone_label"],
            fontsize=FN, ha="right", va="center", color="#555555")
    pi = cr["patient_idx"]
    if pi not in patient_label_positions:
        patient_label_positions[pi] = {"y_top": y_top, "y_bot": y_top - h, "n_total": cr["n_cells"]}
    else:
        patient_label_positions[pi]["y_bot"] = y_top - h
        patient_label_positions[pi]["n_total"] += cr["n_cells"]

for pi in range(n_patients):
    info = patient_label_positions[pi]
    y_mid = (info["y_top"] + info["y_bot"]) / 2
    x_brk = X_LEFT - BLOCK_W / 2 - 0.045
    ax.plot([x_brk + 0.006, x_brk, x_brk, x_brk + 0.006],
            [info["y_top"] - 0.005, info["y_top"] - 0.005,
             info["y_bot"] + 0.005, info["y_bot"] + 0.005],
            color=PATIENT_BASE[pi], linewidth=0.5, clip_on=False, solid_capstyle="round")
    ax.text(x_brk - 0.006, y_mid, SAMPLE_NAMES[pi], fontsize=FN, fontweight="bold",
            ha="right", va="center", color=PATIENT_BASE[pi], rotation=90)

for ep in EPITYPE_ORDER:
    draw_block(ax, X_RIGHT, ep_y_tops[ep], ep_heights[ep], BLOCK_W, EPITYPE_COLORS[ep])

ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.axis("off")
plt.subplots_adjust(left=0.01, right=0.99, top=0.90, bottom=0.01)

fig.savefig(OUT_SVG, format="svg", dpi=450, bbox_inches="tight", pad_inches=0.06)
fig.savefig(OUT_PNG, format="png", dpi=450, bbox_inches="tight", pad_inches=0.06)
plt.close(fig)
print(f"  Saved: {OUT_SVG}\n  Saved: {OUT_PNG}")
