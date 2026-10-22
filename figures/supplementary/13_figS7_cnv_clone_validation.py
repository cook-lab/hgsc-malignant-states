#!/usr/bin/env python3
"""
Figure S7 — CNV validation: within-clone epitype composition
============================================================

Purpose
    Within-clone epitype composition across all samples with CopyKAT subclone
    inference. Each CNV subclone (>= 20 cells) is a stacked bar split into
    SecA / Intermediate / SecB %. Bars grouped by sample, sorted by sample SecB
    fraction descending; a CNV-verdict ribbon sits above each sample group.

INPUTS
    data_root/2026_final_atlas/output/19_cnv/per_sample/<sample>/
        copykat_subclones.csv   (barcode -> clone_k | diploid; clone assignment and the
                                 silhouette-picked k / monoclonal fallback (best silhouette
                                 < 0.15) are taken as-is; CopyKAT is NOT re-run)
        barcodes.csv            (CopyKAT input cells; is_reference flags the diploid reference)
        copykat_subclone_qc.csv (k, silhouette, monoclonal; reported only)
        DONE.txt                (completed samples)
    figures._prep.atlas_final_tables.load_epithelial()  (epithelial `state` SecA / Intermediate /
        SecB / Ciliated and Factor_2, final atlas)
    figures._prep.atlas_final_tables.load_obs()  (sample_id, study, patient_id, celltype_level1:
        recomputed sample manifest)
    The per-clone compositions, sample verdicts and manifest are computed here:
    epithelial (non-reference) aneuploid cells with a state; samples with >= 20 such cells;
    clones >= 20 cells; verdict = monoclonal if only one clone passes, clonally_driven if
    one clone holds >= 70% of the sample's SecA AND one holds >= 70% of its SecB cells,
    else mixed.

OUTPUTS
    output_root/figures/supplementary/SF7_atlas_cnv_validation.{svg,png}
    output_root/_prep_caches/sf7_cnv/{within_clone_coexistence, per_sample_verdict,
        sample_manifest}.csv  ("transitioning" columns = Intermediate; the manifest follows the
        schema of atlas/05_cnv/as_run/19_cnv_extract.py). within_clone_coexistence.csv is read by
        figures/_prep/text_state_coexistence.py (Discussion, 93.5% of patients);
        per_sample_verdict.csv holds the clonality counts of the legend.

MANUSCRIPT PANEL(S)
    Figure S7.

RUNTIME TIER
    fast (reads 251 small per-sample CSVs + the epithelial cache).
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
from figures._prep.atlas_final_tables import load_epithelial, load_obs  # noqa: E402
from scipy.stats import entropy as scipy_entropy  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

PER_SAMPLE = path("data_root", "2026_final_atlas", "output", "19_cnv", "per_sample")
CACHE_DIR = path("output_root", "_prep_caches", "sf7_cnv")
COEX_CSV = os.path.join(CACHE_DIR, "within_clone_coexistence.csv")
VERDICT_CSV = os.path.join(CACHE_DIR, "per_sample_verdict.csv")
MANIFEST_CSV = os.path.join(CACHE_DIR, "sample_manifest.csv")
OUT_SVG = path("output_root", "figures", "supplementary", "SF7_atlas_cnv_validation.svg")
OUT_PNG = path("output_root", "figures", "supplementary", "SF7_atlas_cnv_validation.png")

MIN_CELLS_PER_CLONE = 20      # clones (and samples) with >= 20 aneuploid epithelial cells
DOMINANCE_THRESH = 0.70       # clonally driven: one clone holds >= 70% of SecA and one >= 70% of SecB
MIN_EPI_CELLS, MIN_REF_CELLS, MAX_REF_CELLS = 200, 50, 500   # 19_cnv_extract.py eligibility

# ============================================================================
# STYLE
# ============================================================================

# Epitype display labels. The CSV columns named *_transitioning hold the Intermediate state.
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
# Per-clone compositions and sample verdicts
# ============================================================================

print("Loading epithelial states (load_epithelial)...", flush=True)
epi = load_epithelial(columns=["state", "Factor_2", "sample_id_deposited", "study", "patient_id"]).rename(
    columns={"sample_id_deposited": "sample_id"})  # CopyKAT runs are keyed by the deposited sample id (Nath 2021 fix, 2026-10-20)
STATE_TO_UPSTREAM = {"SecA": "SecA", "Intermediate": "Transitioning", "SecB": "SecB", "Ciliated": "Ciliated"}
epitype = epi["state"].astype(str).map(STATE_TO_UPSTREAM)
f2 = epi["Factor_2"]

completed = sorted(d for d in os.listdir(PER_SAMPLE)
                   if not d.startswith("._") and os.path.exists(os.path.join(PER_SAMPLE, d, "DONE.txt")))
print(f"  completed CopyKAT samples: {len(completed)}")

clone_rows, qc_rows = [], []
for sample_id in completed:
    sdir = os.path.join(PER_SAMPLE, sample_id)
    sub_p = os.path.join(sdir, "copykat_subclones.csv")
    if not os.path.exists(sub_p):
        continue
    # Read barcodes as strings: some studies (e.g. Nath 2021) have numeric barcodes ("1384"),
    # and the atlas is string-indexed.
    subclones = pd.read_csv(sub_p, dtype={"barcode": str})
    barcodes = pd.read_csv(os.path.join(sdir, "barcodes.csv"), dtype={"barcode": str})
    qc_p = os.path.join(sdir, "copykat_subclone_qc.csv")
    if os.path.exists(qc_p):
        qc_rows.append(pd.read_csv(qc_p).assign(sample_id=sample_id))
    ep = barcodes[~barcodes["is_reference"]][["barcode"]]
    merged = ep.merge(subclones, on="barcode", how="inner").set_index("barcode")
    merged["epitype"] = epitype.reindex(merged.index)
    merged["factor2"] = f2.reindex(merged.index)
    merged = merged.dropna(subset=["epitype", "factor2"])
    aneu = merged[merged["subclone"] != "diploid"]
    if len(aneu) < MIN_CELLS_PER_CLONE:
        continue
    n_a = (aneu["epitype"] == "SecA").sum(); n_b = (aneu["epitype"] == "SecB").sum()
    n_t = (aneu["epitype"] == "Transitioning").sum()
    tot = n_a + n_b + n_t
    sample_secB_frac = n_b / tot if tot > 0 else 0
    for clone_id in aneu["subclone"].unique():
        cc = aneu[aneu["subclone"] == clone_id]
        if len(cc) < MIN_CELLS_PER_CLONE:
            continue
        vc = cc["epitype"].value_counts()
        ca, cb, ct = vc.get("SecA", 0), vc.get("SecB", 0), vc.get("Transitioning", 0)
        ctot = ca + cb + ct
        fa, fb, ft = (ca / ctot, cb / ctot, ct / ctot) if ctot > 0 else (0, 0, 0)
        probs = np.array([fa, ft, fb]); probs = probs[probs > 0]
        clone_rows.append({
            "sample_id": sample_id, "subclone": clone_id, "n_cells": len(cc),
            "n_secA": ca, "n_transitioning": ct, "n_secB": cb,
            "frac_secA": fa, "frac_transitioning": ft, "frac_secB": fb,
            "entropy": scipy_entropy(probs, base=2) if len(probs) > 1 else 0,
            "factor2_mean": cc["factor2"].values.mean(), "factor2_std": cc["factor2"].values.std(),
            "has_both_secA_secB": (ca > 0) and (cb > 0), "sample_secB_frac": sample_secB_frac,
        })
coex = pd.DataFrame(clone_rows)
coex.to_csv(COEX_CSV, index=False)

verdict_rows = []
for sid, grp in coex.groupby("sample_id", sort=False):
    n_cl = len(grp)
    if n_cl == 1:
        v = dict(n_clones=1, dom_clone_secA_frac=1.0, dom_clone_secB_frac=1.0, verdict="monoclonal")
    else:
        ta, tb = grp["n_secA"].sum(), grp["n_secB"].sum()
        da = grp["n_secA"].max() / ta if ta > 0 else 0.0
        db = grp["n_secB"].max() / tb if tb > 0 else 0.0
        v = dict(n_clones=n_cl, dom_clone_secA_frac=da, dom_clone_secB_frac=db,
                 verdict="clonally_driven" if (da >= DOMINANCE_THRESH and db >= DOMINANCE_THRESH) else "mixed")
    v.update(sample_id=sid, sample_secB_frac=grp["sample_secB_frac"].iloc[0],
             within_clone_entropy_mean=float(grp["entropy"].mean()))
    verdict_rows.append(v)
verdict = pd.DataFrame(verdict_rows)[["sample_id", "n_clones", "sample_secB_frac", "dom_clone_secA_frac",
                                      "dom_clone_secB_frac", "within_clone_entropy_mean", "verdict"]]
qc = pd.concat(qc_rows, ignore_index=True) if qc_rows else pd.DataFrame()
if len(qc):
    qc["sample_id"] = qc["sample_id"].astype(str)
    verdict = verdict.merge(qc[["sample_id", "k", "silhouette", "monoclonal"]].rename(
        columns={"k": "copykat_k", "silhouette": "copykat_silhouette", "monoclonal": "copykat_qc_monoclonal"}),
        on="sample_id", how="left")
verdict.to_csv(VERDICT_CSV, index=False)

# sample manifest (19_cnv_extract.py schema) from the final-atlas metadata
obs_m = load_obs(columns=["sample_id_deposited", "celltype_level1"]).rename(columns={"sample_id_deposited": "sample_id"})
n_ref_raw = obs_m[obs_m["celltype_level1"].isin(["Fibroblast", "T/NK cell", "Endothelial"])] \
    .groupby("sample_id", observed=True).size()
st = epi.groupby("sample_id", observed=True)
manifest = pd.DataFrame({
    "study": st["study"].first().astype(str), "patient_id": st["patient_id"].first().astype(str),
    "n_epithelial": st.size(),
})
for lab, colname in [("SecA", "n_secA"), ("Intermediate", "n_transitioning"), ("SecB", "n_secB"), ("Ciliated", "n_ciliated")]:
    manifest[colname] = epi[epi["state"] == lab].groupby("sample_id", observed=True).size().reindex(manifest.index).fillna(0).astype(int)
manifest["n_reference"] = np.minimum(n_ref_raw.reindex(manifest.index).fillna(0).astype(int), MAX_REF_CELLS)
manifest["status"] = np.where((manifest["n_epithelial"] >= MIN_EPI_CELLS) & (manifest["n_reference"] >= MIN_REF_CELLS),
                              "eligible", "skipped")
manifest["selected"] = manifest["status"] == "eligible"
manifest.index.name = "sample_id"
manifest = manifest.reset_index()[["sample_id", "study", "patient_id", "n_epithelial", "n_reference", "n_secA",
                                   "n_transitioning", "n_secB", "n_ciliated", "status", "selected"]]
manifest.to_csv(MANIFEST_CSV, index=False)
vc_ = verdict["verdict"].value_counts().to_dict()
print(f"  recomputed: {len(coex)} clones in {coex['sample_id'].nunique()} samples; verdicts {vc_}; "
      f"manifest {len(manifest)} samples ({int(manifest['selected'].sum())} eligible)")

coex["sample_id"] = pd.to_numeric(coex["sample_id"]).astype(int)   # CopyKAT dirs are e.g. "161.0"
verdict["sample_id"] = pd.to_numeric(verdict["sample_id"]).astype(int)   # CopyKAT dirs are e.g. "161.0"
manifest["sample_id"] = pd.to_numeric(manifest["sample_id"]).astype(int)   # CopyKAT dirs are e.g. "161.0"
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

# Split the samples (kept intact, in SecB order) into NROWS rows of near-equal bar counts (~63
# bars/row): cut at the sample boundary closest to each k/NROWS of the bars.
NROWS = max(4, int(np.ceil(n_bars / 63)))
cum = np.cumsum([len(recs) for _, recs in sample_groups])
cuts, prev_cut = [], 0
for k in range(1, NROWS):
    cut = int(np.argmin(np.abs(cum - k * n_bars / NROWS))) + 1   # cut after this many samples
    cut = min(max(cut, prev_cut + 1), len(sample_groups) - (NROWS - k))
    cuts.append(cut)
    prev_cut = cut
bounds = [0] + cuts + [len(sample_groups)]
rows_of_samples = [sample_groups[a:b] for a, b in zip(bounds[:-1], bounds[1:])]
assert len(rows_of_samples) == NROWS and all(rows_of_samples)
print(f"  Layout: {NROWS} rows across letter page; bars per row "
      f"{[sum(len(r) for _, r in row) for row in rows_of_samples]}")
# One x-scale for every row, so bars are the same width on each row (a shorter row ends early).
X_MAX = max(sum(len(r) for _, r in row) - 1 + INTER_GAP * (len(row) - 1) for row in rows_of_samples)

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

    ax.set_xlim(-0.6, X_MAX + 0.6)
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
