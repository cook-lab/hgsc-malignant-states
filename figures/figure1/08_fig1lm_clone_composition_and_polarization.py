#!/usr/bin/env python3
"""
Figure 1l & 1m — Per-clone epithelial state composition and polarization co-existence
=====================================================================================
PURPOSE
    Two panels built from one shared per-clone prep, for three representative
    atlas tumours, one monoclonal and two polyclonal (the same tumours, so l and
    m read together):

      Fig 1l — per CopyKAT CNV subclone, the SecA / Intermediate / SecB state
               COMPOSITION (stacked %), grouped by tumour. Shows clones differ in
               propensity yet each contains all three states.
      Fig 1m — per subclone, the DISTRIBUTION of the secretory-polarization score
               (NMF Factor_2 = SecB) as a violin, grouped by tumour, with dashed
               atlas thresholds (Intermediate p50 / SecB p75). Each clone SPANS
               the SecA->SecB axis: states co-exist WITHIN clones.

    Three renders of the same data: l and m as standalone figures, and l and m
    stacked on one shared clone x-axis (violins on top). The figure page was
    assembled by hand from these renders.

    Also checks, over every completed CopyKAT sample, that the CopyKAT barcodes
    join the final-atlas epithelial table completely (asserted; printed).

INPUTS
    - Per-cell metadata, Factor_2 usage and epithelial state: the final-atlas
        loader figures/_prep/atlas_final_tables.py :: load_epithelial()
        (575,366 epithelial cells; Factor_1..10 from 11d_nmf_usage.csv; `state`
        SecA/Intermediate/SecB/Ciliated; cache output_root/_prep_caches/).
        Atlas polarization thresholds p50/p75 = Factor_2 percentiles among
        non-ciliated epithelial cells (config polarization.partition).
    - CopyKAT outputs (data_root/2026_final_atlas/output/19_cnv/per_sample/<id>/):
        copykat_subclones.csv (barcode, subclone); DONE.txt marks completed
        samples. Epithelial (non-reference) cells = subclone barcodes present in
        the loader's epithelial set (asserted: every non-reference CopyKAT
        barcode is in the final atlas epithelial set; reference cells never are).
        barcodes.csv is read only for is_reference (coverage assertion).
    - Representative samples (CopyKAT per-sample ids):
        SAMPLE_IDS = [16 (monoclonal), 100 (2 clones), 62 (4 clones)].

OUTPUTS
    - figures_dir/figure1/atlas_clone_polarization_and_composition_lm.{svg,png}   (Fig 1l and 1m, stacked as published)
    - output_root/figure1/atlas_copykat_coverage.csv   (per-sample barcode coverage behind the assertion)

MANUSCRIPT PANEL(S): Fig 1l (composition), Fig 1m (violins).

RUNTIME TIER: fast (one cached parquet + ~250 small per-sample CSVs).

NOTE: state labels = loader `state` (SecA / Intermediate / SecB). Ciliated cells
    are not counted in the composition denominator. Factor = config
    polarization.factor (Factor_2). In the code, the *_K / figK / axK names are
    the composition figure (panel l), the *_L / figL / axL names the violin
    figure (panel m), and axl / axm the top (violins) and bottom (composition)
    axes of the stacked render.
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
from matplotlib.colors import to_rgb, to_hex

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import CFG, path  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

# ---------- Paths ----------
PER_SAMPLE = path("data_root", "2026_final_atlas", "output", "19_cnv", "per_sample")
# Stacked render: m (violins, top) over l (composition, bottom) sharing one x-axis.
OUT_LM_SVG = path("figures_dir", "figure1", "atlas_clone_polarization_and_composition_lm.svg")
OUT_LM_PNG = path("figures_dir", "figure1", "atlas_clone_polarization_and_composition_lm.png")
OUT_COV_CSV = path("output_root", "figure1", "atlas_copykat_coverage.csv")

assert os.path.isdir(PER_SAMPLE), f"Missing input: {PER_SAMPLE}"

FACTOR = CFG["polarization"]["factor"]  # Factor_2 (SecB-defining polarization axis)

# ---------- Representative samples ----------
SAMPLE_IDS = ["16", "100", "62"]
SAMPLE_NAMES = ["Patient A\n(monoclonal)", "Patient B", "Patient C"]
PATIENT_BASE = ["#8C8C8C", "#4A7FA5", "#4A967A"]

# State palette matched to Fig 1k (02_fig1bcdek_atlas_study_contribution.py) for consistency.
EPI_COLORS = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55"}
EPI_ORDER = ["SecA", "Intermediate", "SecB"]
MIN_CELLS_PER_CLONE = 20

# ---------- Style ----------
FA, FK, FN = 8, 6, 8
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
    "figure.dpi":       450,
    "savefig.dpi":      450,
})


def shade_variants(base_hex, n, spread=0.28):
    r, g, b = to_rgb(base_hex)
    if n == 1:
        return [base_hex]
    out = []
    for i in range(n):
        t = (i / (n - 1)) - 0.5
        f = 1.0 + t * spread * 2
        out.append(to_hex((max(0, min(1, r * f)), max(0, min(1, g * f)), max(0, min(1, b * f)))))
    return out


# ============================================================================
# SHARED PREP
# ============================================================================
epi_all = load_epithelial(columns=["sample_id_deposited", "state", FACTOR]).rename(
    columns={"sample_id_deposited": "sample_id"})  # CopyKAT runs are keyed by the deposited sample id (Nath 2021 fix)
assert len(epi_all) == 575_366, len(epi_all)
non_cil = epi_all.loc[epi_all["state"] != "Ciliated", FACTOR].to_numpy(float)
P50 = float(np.percentile(non_cil, 50))
P75 = float(np.percentile(non_cil, 75))
print(f"[thresholds] non-ciliated {FACTOR}: p50={P50:.4f}  p75={P75:.4f}")


def load_clones(sid):
    """Non-reference CopyKAT cells of one sample joined to the final-atlas epithelial
    table; returns (aneuploid-cell frame, coverage dict)."""
    sdir = os.path.join(PER_SAMPLE, sid)
    sc = pd.read_csv(os.path.join(sdir, "copykat_subclones.csv"), dtype={"barcode": str})
    bc = pd.read_csv(os.path.join(sdir, "barcodes.csv"), dtype={"barcode": str})
    ref = set(bc.loc[bc["is_reference"], "barcode"])
    nonref = sc[~sc["barcode"].isin(ref)]
    in_epi = nonref["barcode"].isin(epi_all.index)
    cov = {"sample_id": sid, "n_copykat_nonref": len(nonref), "n_in_atlas_epithelial": int(in_epi.sum()),
           "n_reference_in_atlas_epithelial": int(pd.Index(list(ref)).isin(epi_all.index).sum())}
    merged = nonref.set_index("barcode").join(epi_all, how="inner")
    cov["n_sample_id_mismatch"] = int((merged["sample_id"].astype(str) != sid).sum())
    merged = merged[merged["subclone"] != "diploid"].copy()
    merged = merged.rename(columns={FACTOR: "factor2", "state": "epi_std"}).dropna(subset=["factor2"])
    return merged, cov


clones = []   # patient_idx, clone_label, values, frac{}, n
for pi, sid in enumerate(SAMPLE_IDS):
    merged, _ = load_clones(sid)
    keep = [c for c in sorted(merged["subclone"].unique())
            if (merged["subclone"] == c).sum() >= MIN_CELLS_PER_CLONE]
    for ci, cl in enumerate(keep):
        cc = merged[merged["subclone"] == cl]
        comp = cc["epi_std"].astype(str).value_counts()
        tot = comp.reindex(EPI_ORDER).fillna(0).sum()
        clones.append({
            "patient_idx": pi, "clone_label": f"C{ci + 1}",
            "values": cc["factor2"].to_numpy(float), "n": int(len(cc)),
            "frac": {e: (comp.get(e, 0) / tot if tot else 0) for e in EPI_ORDER},
        })
    print(f"  sample {sid}: matched {len(merged)} cells, {len(keep)} clone(s) (>= {MIN_CELLS_PER_CLONE})")

# x positions: clones grouped by tumour, gaps between tumours
GROUP_GAP = 1.0
xs, vcolors, xlabels, ns = [], [], [], []
group_spans = {}
x = 0.0
prev_pi = clones[0]["patient_idx"]
per_patient = {pi: sum(1 for c in clones if c["patient_idx"] == pi) for pi in range(len(SAMPLE_IDS))}
shades = {pi: shade_variants(PATIENT_BASE[pi], per_patient[pi]) for pi in per_patient}
ci_counter = {pi: 0 for pi in per_patient}
for c in clones:
    pi = c["patient_idx"]
    if pi != prev_pi:
        x += GROUP_GAP
        prev_pi = pi
    xs.append(x)
    vcolors.append(shades[pi][ci_counter[pi]])
    ci_counter[pi] += 1
    xlabels.append(c["clone_label"])
    ns.append(c["n"])
    sp = group_spans.setdefault(pi, [x, x])
    sp[1] = x
    x += 1.0
XMIN, XMAX = min(xs) - 0.6, max(xs) + 0.6


def add_group_labels(ax, y_line=-0.30, y_text=-0.40):
    """Tumour bracket + name below an axis (axes-fraction y)."""
    for pi, (x0, x1) in group_spans.items():
        xm = (x0 + x1) / 2
        ax.plot([x0 - 0.3, x1 + 0.3], [y_line, y_line], transform=ax.get_xaxis_transform(),
                color=PATIENT_BASE[pi], lw=0.7, clip_on=False, solid_capstyle="round")
        ax.text(xm, y_text, SAMPLE_NAMES[pi], transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=FN, color=PATIENT_BASE[pi],
                fontweight="bold", linespacing=0.85)


# ---------- CopyKAT barcode coverage (all completed CopyKAT samples) ----------
# Records the barcode coverage of every CopyKAT sample against the final atlas
# epithelial set and asserts that the join is complete.
done = sorted(d for d in os.listdir(PER_SAMPLE)
              if not d.startswith("._") and os.path.isdir(os.path.join(PER_SAMPLE, d))
              and os.path.exists(os.path.join(PER_SAMPLE, d, "DONE.txt"))
              and os.path.exists(os.path.join(PER_SAMPLE, d, "copykat_subclones.csv")))
clone_rows, cov_rows = [], []
for sid in done:
    merged, cov = load_clones(sid)
    cov_rows.append(cov)
cov_df = pd.DataFrame(cov_rows)
cov_df.to_csv(OUT_COV_CSV, index=False)
n_nonref, n_hit = int(cov_df["n_copykat_nonref"].sum()), int(cov_df["n_in_atlas_epithelial"].sum())
print(f"[coverage] {len(cov_df)} CopyKAT samples: {n_hit:,}/{n_nonref:,} non-reference barcodes in the "
      f"final atlas epithelial set ({100 * n_hit / n_nonref:.2f}%); reference barcodes in epithelial set: "
      f"{int(cov_df['n_reference_in_atlas_epithelial'].sum())}; sample_id mismatches: "
      f"{int(cov_df['n_sample_id_mismatch'].sum())}")
assert n_hit == n_nonref, "CopyKAT barcodes missing from the final atlas epithelial set"
assert cov_df["n_sample_id_mismatch"].sum() == 0 and cov_df["n_reference_in_atlas_epithelial"].sum() == 0


# ============================================================================
# Fig 1l + 1m stacked — violins (m, TOP) over state composition (l, BOTTOM),
# sharing ONE CopyKAT clone x-axis.
# ============================================================================
figLM, (axl, axm) = plt.subplots(
    2, 1, sharex=True, figsize=(3.0, 3.10),
    gridspec_kw={"height_ratios": [1, 1], "hspace": 0.08},
)
figLM.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.32)

# --- TOP (panel m): polarization violins, atlas thresholds preserved ---
axl.axhline(P50, ls="--", lw=0.5, color="#888888", zorder=1)
axl.axhline(P75, ls="--", lw=0.5, color="#555555", zorder=1)
parts = axl.violinplot([c["values"] for c in clones], positions=xs, widths=0.82,
                       showmeans=False, showextrema=False, showmedians=False)
for body, col in zip(parts["bodies"], vcolors):
    body.set_facecolor(col); body.set_edgecolor("white"); body.set_linewidth(0.4)
    body.set_alpha(0.92); body.set_zorder(2)
for xi, c in zip(xs, clones):
    med = float(np.median(c["values"]))
    axl.plot([xi - 0.18, xi + 0.18], [med, med], color="#2b2b2b", lw=0.7, zorder=3)
ymax = max(np.percentile(c["values"], 99) for c in clones)
axl.set_ylim(0, max(ymax * 1.05, P75 * 1.15))
axl.set_ylabel(f"SecB polarization\n({FACTOR})", fontsize=FA, labelpad=2, linespacing=0.9)
axl.tick_params(axis="y", length=2, width=0.4, pad=1)
axl.tick_params(axis="x", length=0)          # sharex hides labels; hide top ticks too
for sp in ("top", "right"):
    axl.spines[sp].set_visible(False)
axl.spines["left"].set_linewidth(0.5); axl.spines["bottom"].set_linewidth(0.5)

# --- BOTTOM (panel l): per-clone state composition (stacked %) ---
for xi, c in zip(xs, clones):
    bottom = 0.0
    for e in EPI_ORDER:
        v = c["frac"][e] * 100
        axm.bar(xi, v, bottom=bottom, width=0.82, color=EPI_COLORS[e],
                edgecolor="white", linewidth=0.4)
        bottom += v
axm.set_ylim(0, 100); axm.set_yticks([0, 25, 50, 75, 100])
axm.set_ylabel("Epitype (%)", fontsize=FA, labelpad=2)
axm.set_xlim(XMIN, XMAX)
axm.set_xticks(xs); axm.set_xticklabels(xlabels, fontsize=FK)
axm.tick_params(axis="y", length=2, width=0.4, pad=1)
axm.tick_params(axis="x", length=0, pad=1)
for sp in ("top", "right"):
    axm.spines[sp].set_visible(False)
axm.spines["left"].set_linewidth(0.5); axm.spines["bottom"].set_linewidth(0.5)

# Shared x-axis annotations, ONCE at the very bottom: n= counts under the clone
# labels (angled ~45deg so neighbouring counts don't collide at 3.05" width), then
# the Patient A/B/C brackets + names below that.
for xi, c in zip(xs, clones):
    axm.text(xi, -0.13, f"n={c['n']}", transform=axm.get_xaxis_transform(),
             ha="right", va="top", rotation=45, rotation_mode="anchor",
             fontsize=FN, color="#555555")
add_group_labels(axm, y_line=-0.52, y_text=-0.59)

figLM.savefig(OUT_LM_SVG, format="svg")
figLM.savefig(OUT_LM_PNG, format="png", dpi=450)
plt.close(figLM)
print(f"  Saved: {OUT_LM_SVG}\n  Saved: {OUT_LM_PNG}")
