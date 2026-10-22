#!/usr/bin/env python3
"""
Figure S15 -- Secretory polarization score across the Xenium cohort (final SingleR cohort)
==========================================================================================

Purpose
    Per-cell secretory polarization (SecB UCell - SecA UCell, `polarization_UCell`) of EVERY
    secretory epithelial cell (SingleR "Secretory epithelium" = cell_label_singler in
    SecA / Intermediate / SecB), no subsampling, across the three Xenium arms:
        Fallopian tube (FTE) cores   -- one violin per core (15; one core per donor)
        HGSC whole tissues            -- one violin per section (core-8)
        HGSC TMA                      -- one violin per patient, cores pooled (TMA A + TMA B)
    (a) one violin per unit, all arms on one x-axis sorted by the unit median;
    (b) all secretory cells pooled per arm (violin + box).
    Dashed lines = the frozen 06f thresholds that split secretory cells into SecA
    (<= 0.109429) / Intermediate / SecB (>= 0.210501) -- the same thresholds that define
    cell_label_singler. No atlas component (the published panel has none).

Cohort
    TMA: the final TMA core set (spatial/02_qc/06b_final_tma_cores.py via figures/_prep/final_tma_cores.py;
    author decision 2026-10-19), the cores behind every TMA panel: tumour cores passing per-array core QC
    on the SingleR labels, from patients with primary HGSC, treatment-naive at sampling, with >= 500 cells
    on those cores. TMA A (CHUM) 95 + TMA B (COEUR G+H) 283 = 378 patients, 503 cores. Every cohort unit is
    plotted; units with < 5 secretory cells are drawn as points (a violin needs a spread).

Inputs
    $OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz
    $OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv  (codes only; never written)
Outputs  (OUTPUT_ROOT/figures/supplementary/SF13_singler/)
    SuppFig15_xenium_polarization_distribution.{png,svg}
    SF13_units.csv   (de-identified unit code, arm, n secretory cells, median, IQR)
    SF13_arms.csv    (per-arm n units / cells / median / IQR)
    SF13_tests.csv   FTE cores vs HGSC units (2026-10-20, for "higher in HGSC than FTE"): two-sided Mann-Whitney
                     on per-unit median polarization and on per-unit % SecB; units = FTE core (one per donor),
                     whole tissue, TMA patient; units without secretory cells are left out. Copied to
                     final_panels/supplementary/tables/SuppFig15_fte_vs_hgsc_polarization_tests.csv
                     (Results: p = 0.003, per-unit medians, all HGSC units vs FTE cores)
    The folder and table names carry the earlier figure number (SF13 = Figure S15).
MANUSCRIPT PANEL(S): Figure S15a, b; the Results p value above.

Run:  python figures/supplementary/20_figS15_polarization_distribution.py
Runtime: ~1 min.
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT_ROOT = os.environ.get("OUTPUT_ROOT", "output")
SG = os.path.join(OUT_ROOT, "singler_recompute")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from config.study_codes import code_of   # source ID or code -> study code (private key; never written)
from figures._prep.final_tma_cores import final_core_ids   # the final TMA core set
OUT = os.path.join(OUT_ROOT, "figures", "supplementary", "SF13_singler")
os.makedirs(OUT, exist_ok=True)
STEM = "SuppFig15_xenium_polarization_distribution"

SEC = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
T_LOW, T_HIGH = 0.109429, 0.210501          # frozen 06f thresholds (cell_label_singler)
WT8 = ["WT01", "WT02", "WT03", "WT04", "WT05", "WT06", "WT07", "WT08"]
ARMS = ["Fallopian tube", "HGSC whole tissue", "HGSC TMA"]
COL = {"Fallopian tube": "#E07850", "HGSC whole tissue": "#8A5DAF", "HGSC TMA": "#D4A574"}  # Fig 4a
TMA_OBJ = {"tma": "A", "COEUR_TMA_G": "B", "COEUR_TMA_H": "B"}

FA, FK, FN = 6, 5.5, 5
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": FA, "axes.labelsize": FA, "xtick.labelsize": FK, "ytick.labelsize": FK,
    "legend.fontsize": FN, "savefig.dpi": 450, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.4, "xtick.major.width": 0.4, "ytick.major.width": 0.4,
})


def rd(nm, cols):
    return pd.read_csv(os.path.join(SG, f"{nm}_singler_cells.csv.gz"), usecols=cols,
                       dtype={"patient_id": str, "core_id": str})


# ---------------------------------------------------------------- whole tissues
units = []                                                # (code, arm, values)
for s in WT8:
    d = rd(s, ["cell_label_singler", "polarization_UCell"])
    v = d.loc[d.cell_label_singler.isin(SEC), "polarization_UCell"].to_numpy()
    assert not np.isnan(v).any()
    units.append((code_of("Whole tissue", s), "HGSC whole tissue", v))

# ---------------------------------------------------------------- TMA + FTE
cols = ["cell_label_singler", "polarization_UCell", "patient_id", "core_id", "sample_type"]
tma = []
for nm, t in TMA_OBJ.items():
    d = rd(nm, cols)
    d["TMA"] = t
    d["final_core"] = (d.sample_type == "tumour") & d.core_id.isin(final_core_ids(nm))
    d["core_key"] = nm + ":" + d.core_id.astype(str)
    tma.append(d)
tma = pd.concat(tma, ignore_index=True)

tum = tma[tma.final_core]
assert tum.patient_id.notna().all()
n_cores = tum.core_key.nunique()
sec_t = tum[tum.cell_label_singler.isin(SEC)]
assert not sec_t.polarization_UCell.isna().any()
for (t, p), g in tum.groupby(["TMA", "patient_id"]):
    v = g.loc[g.cell_label_singler.isin(SEC), "polarization_UCell"].to_numpy()
    units.append((code_of(f"TMA {t}", p), "HGSC TMA", v))

fte = tma[tma.sample_type == "fallopian"]
for c, g in fte.groupby("core_id"):
    v = g.loc[g.cell_label_singler.isin(SEC), "polarization_UCell"].to_numpy()
    units.append((code_of("FTE core", c), "Fallopian tube", v))

n_by_arm = {a: sum(u[1] == a for u in units) for a in ARMS}
print("units per arm:", n_by_arm, f"| TMA cores: {n_cores}")
assert n_by_arm == {"Fallopian tube": 15, "HGSC whole tissue": 8, "HGSC TMA": 378}, n_by_arm
assert n_cores == 503

# ---------------------------------------------------------------- tables
ut = pd.DataFrame([{"unit": c, "arm": a, "n_secretory_cells": len(v),
                    "median": np.median(v) if len(v) else np.nan,
                    "q25": np.quantile(v, .25) if len(v) else np.nan,
                    "q75": np.quantile(v, .75) if len(v) else np.nan,
                    "pct_SecB": 100 * np.mean(v >= T_HIGH) if len(v) else np.nan} for c, a, v in units])
ut.sort_values(["arm", "unit"]).to_csv(os.path.join(OUT, "SF13_units.csv"), index=False)
pooled = {a: np.concatenate([u[2] for u in units if u[1] == a]) for a in ARMS}
at = pd.DataFrame([{"arm": a, "units": n_by_arm[a], "secretory_cells": len(pooled[a]),
                    "median": np.median(pooled[a]), "q25": np.quantile(pooled[a], .25),
                    "q75": np.quantile(pooled[a], .75),
                    "median_of_unit_medians": ut.loc[ut.arm == a, "median"].median(),
                    "units_lt10_sec_cells": int((ut.loc[ut.arm == a, "n_secretory_cells"] < 10).sum())}
                   for a in ARMS])
at.to_csv(os.path.join(OUT, "SF13_arms.csv"), index=False)

# ---------------------------------------------------------------- FTE vs HGSC (unit level)
uv = ut.dropna(subset=["median"])                     # units with >= 1 secretory cell
fte_u = uv[uv.arm == "Fallopian tube"]
trows = []
for lab, h in [("all HGSC units (whole tissues + TMA patients)", uv[uv.arm != "Fallopian tube"]),
               ("HGSC whole tissues", uv[uv.arm == "HGSC whole tissue"]), ("HGSC TMA patients", uv[uv.arm == "HGSC TMA"])]:
    for metric in ["median", "pct_SecB"]:
        r = mannwhitneyu(h[metric], fte_u[metric], alternative="two-sided")
        trows.append({"comparison": f"{lab} vs FTE cores", "metric": "median polarization per unit" if metric == "median"
                      else "% SecB cells per unit", "n_HGSC_units": len(h), "n_FTE_cores": len(fte_u),
                      "HGSC_median_of_units": round(float(h[metric].median()), 4),
                      "FTE_median_of_units": round(float(fte_u[metric].median()), 4),
                      "AUC_HGSC_gt_FTE": round(float(r.statistic) / (len(h) * len(fte_u)), 3),
                      "p_mann_whitney": float(f"{r.pvalue:.3g}")})
tt = pd.DataFrame(trows)
tt.to_csv(os.path.join(OUT, "SF13_tests.csv"), index=False)
print(tt.to_string(index=False))
print(at.round(4).to_string(index=False))
print("TMA units with < 5 / < 10 secretory cells:",
      int((ut[ut.arm == "HGSC TMA"].n_secretory_cells < 5).sum()),
      int((ut[ut.arm == "HGSC TMA"].n_secretory_cells < 10).sum()))

allv = np.concatenate(list(pooled.values()))
ylo, yhi = np.floor(allv.min() * 20) / 20, np.ceil(allv.max() * 20) / 20
print(f"cells plotted: {len(allv):,}; range {allv.min():.3f} .. {allv.max():.3f}")
for a in ARMS:
    print(f"  {a}: {100 * np.mean(pooled[a] == 0):.1f}% of secretory cells have polarization exactly 0; "
          f"{sum((u[1] == a) and len(u[2]) and np.median(u[2]) == 0 for u in units)} units with median exactly 0")

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(8.5, 2.5))
gs = fig.add_gridspec(1, 2, width_ratios=[6.4, 1.5], wspace=0.06)
ax = fig.add_subplot(gs[0]); axb = fig.add_subplot(gs[1], sharey=ax)


def thresholds(a, label):
    for t in (T_LOW, T_HIGH):
        a.axhline(t, color="0.35", lw=0.45, ls=(0, (3, 2)), zorder=0)
    a.axhline(0, color="black", lw=0.35, alpha=0.3, zorder=0)
    if label:
        x = a.get_xlim()[1]
        for y, s, va in [(T_LOW, "SecA", "top"), ((T_LOW + T_HIGH) / 2, "Int", "center"), (T_HIGH, "SecB", "bottom")]:
            a.text(x, y + (-0.012 if va == "top" else 0.012 if va == "bottom" else 0), s,
                   ha="right", va=va, fontsize=FN - 0.5, color="0.3")


# (a) one violin per unit, sorted by median
srt = sorted(units, key=lambda u: np.median(u[2]) if len(u[2]) else np.inf)
pos = np.arange(1, len(srt) + 1)
big = [(p, u) for p, u in zip(pos, srt) if len(u[2]) >= 5]
vp = ax.violinplot([u[2] for _, u in big], positions=[p for p, _ in big], widths=0.95,
                   showextrema=False, points=200)
for body, (_, u) in zip(vp["bodies"], big):
    body.set_facecolor(COL[u[1]]); body.set_edgecolor("none"); body.set_alpha(0.9)
for p, u in zip(pos, srt):
    if len(u[2]) >= 5:
        m = np.median(u[2]); ax.plot([p - .32, p + .32], [m, m], color="black", lw=0.45, zorder=5)
    elif len(u[2]):
        ax.scatter(np.full(len(u[2]), p), u[2], s=1.2, color=COL[u[1]], edgecolor="none", zorder=4)
ax.set_xlim(0.2, len(srt) + 4)
ax.set_ylim(ylo, yhi)
thresholds(ax, False)
ax.set_xticks([])
ax.set_ylabel("Secretory polarization score\n(SecB − SecA UCell)")
for s in ("top", "right", "bottom"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=2)
ax.annotate("", xy=(len(srt), -0.04), xytext=(1, -0.04), xycoords=("data", "axes fraction"),
            arrowprops=dict(arrowstyle="->", color="grey", lw=0.6))
ax.text(len(srt) / 2, -0.07, "one violin per FTE core / whole tissue / TMA patient, ordered by median polarization (SecA → SecB)", color="grey",
        fontsize=FN - 0.5, ha="center", va="top", style="italic", transform=ax.get_xaxis_transform())
lab = {"Fallopian tube": f"Fallopian tube cores (n={n_by_arm['Fallopian tube']})",
       "HGSC whole tissue": f"HGSC whole tissue (n={n_by_arm['HGSC whole tissue']})",
       "HGSC TMA": f"HGSC TMA patients (n={n_by_arm['HGSC TMA']}; {n_cores} cores)"}
ax.legend(handles=[Patch(facecolor=COL[a], edgecolor="none", label=lab[a]) for a in ARMS],
          loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, frameon=False, handlelength=1.2,
          handleheight=0.8, columnspacing=1.5, borderaxespad=0.2)

# (b) pooled per arm
bp_pos = np.arange(1, 4)
vb = axb.violinplot([pooled[a] for a in ARMS], positions=bp_pos, widths=0.8, showextrema=False, points=300)
for body, a in zip(vb["bodies"], ARMS):
    body.set_facecolor(COL[a]); body.set_edgecolor("none"); body.set_alpha(0.9)
axb.boxplot([pooled[a] for a in ARMS], positions=bp_pos, widths=0.18, showfliers=False, patch_artist=True,
            boxprops=dict(facecolor="white", lw=0.4), medianprops=dict(color="black", lw=0.6),
            whiskerprops=dict(lw=0.4), capprops=dict(lw=0.4))
axb.set_xlim(0.4, 3.6)
thresholds(axb, False)
axb.set_xticks(bp_pos)
axb.set_xticklabels([f"{t}\n{len(pooled[a]):,}" for t, a in zip(["FTE", "WT", "TMA"], ARMS)], fontsize=FN - 0.5)
axb.text(0.5, -0.155, "pooled secretory cells", ha="center", va="top", fontsize=FN - 0.5, color="0.3",
         transform=axb.transAxes)
for y, lab_, va in [(T_LOW, "SecA", "top"), ((T_LOW + T_HIGH) / 2, "Int", "center"), (T_HIGH, "SecB", "bottom")]:
    axb.text(3.65, y + (-0.01 if va == "top" else 0.01 if va == "bottom" else 0), lab_, ha="left", va=va,
             fontsize=FN - 0.5, color="0.3", clip_on=False)
axb.tick_params(axis="y", labelleft=False, length=0)
axb.tick_params(axis="x", length=0)
for s in ("top", "right", "left"):
    axb.spines[s].set_visible(False)

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"{STEM}.{ext}"), dpi=450, bbox_inches="tight", facecolor="white")
plt.close(fig)
for ext in ("png", "svg"):
    f = os.path.join(OUT, f"{STEM}.{ext}")
    print(f"-> {f} ({os.path.getsize(f) / 1e6:.1f} MB)")
