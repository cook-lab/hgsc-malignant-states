#!/usr/bin/env python3
"""
Figure 2b–e — cell cycle, merged composition, patient-matched shifts
====================================================================
PURPOSE
    The four panels of Figure 2b–e, drawn with identical plot heights so they
    align on the page:

      b  cell-cycle phase of SecA vs SecB cells (scored on the final atlas,
         figures/_prep/fig2b_cell_cycle_phase.py)
      c  epithelial-state composition: treatment-naive Prim / Asc / Met
         (cells pooled per site) | primary tumours Pre / Post-chemo
      d  patient-matched SecB % across sites (patients with ascites and >=1 other
         site, >=100 epithelial cells per site; paired Wilcoxon brackets)
      e  patient-matched SecB % of the epithelium before vs after chemotherapy
         (>=100 epithelial cells per time point, samples pooled per time point,
         not site-matched; paired Wilcoxon, exact p shown)

    No panel letters or panel titles (added in layout); group labels
    ("Treatment-naive", "Primary tumours") and "Patient-matched (n = ...)" kept.
    Site names abbreviated on c (Prim, Asc, Met); expanded in the legend.
    DATA SOURCE: all cell metadata (patient, sample, treatment, anatomic/metastatic
    site, state) come from figures/_prep/atlas_final_tables.load_epithelial()
    (final atlas + deposited NMF usage).
      * b: output_root/figure2/fig2b_cell_cycle/cell_cycle_phase_proportions.csv
      * c site bars: treatment-naive cells pooled per final-atlas metastatic_site.
      * c treatment bars: epithelial cells of primary tumours (metastatic_site ==
        "primary", the column of the site bars), pre-treatment vs post-chemotherapy,
        pooled; sample/patient n written to row2_stats.csv.

OUTPUTS (figures_dir/figure2/)
    row2_b_cell_cycle, row2_c_composition, row2_d_paired_sites,
    row2_e_paired_treatment_secb_pct  .{svg,png};
    row2_stats.csv (n, pooled percentages and paired tests of b-e);
    row2_e_variant_secb_pct_stats.csv (panel e: n, n increased, medians, p)

MANUSCRIPT PANEL(S): Fig 2b, 2c, 2d, 2e.
INPUTS
    figures/_prep/atlas_final_tables.py (load_epithelial)
    output_root/figure2/fig2b_cell_cycle/cell_cycle_phase_proportions.csv
        (figures/_prep/fig2b_cell_cycle_phase.py)

RUNTIME TIER: fast (cached loader parquet).
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
from matplotlib.patches import Patch
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import load_epithelial  # noqa: E402

# ---------- Paths ----------
CC_CSV = path("output_root", "figure2", "fig2b_cell_cycle", "cell_cycle_phase_proportions.csv")
OUT = lambda name: path("figures_dir", "figure2", name)  # noqa: E731
assert os.path.exists(CC_CSV), f"Missing: {CC_CSV} (run figures/_prep/fig2b_cell_cycle_phase.py)"
MIN_CELLS_PAIRED = 100

# ---------- Style (Figure 2 rule) ----------
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "svg.fonttype": "none", "pdf.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "savefig.dpi": 450,
})
PALETTE = {"SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55", "Ciliated": "#E07850"}
STACK = ["SecA", "Intermediate", "SecB", "Ciliated"]
CC_COLS = {"SecA": ["#DDD5C8", "#F0C891", "#E6A141"], "SecB": ["#DDD5C8", "#C9B49A", "#9A7D55"]}
LINE_C = "#9A7D55"
# Patient lines (d, e): coloured per segment by direction, atlas cell-type palette (author 2026-10-05):
UP_C, DOWN_C = "#D14E6C", "#87CEFA"   # Smooth muscle red = increase; T/NK blue = decrease


def patient_line(ax, xs, vs):
    """Draw one patient's linked values; each segment red if it rises, blue if it falls."""
    xs, vs = np.asarray(xs, float), np.asarray(vs, float)
    m = ~np.isnan(vs)
    xs, vs = xs[m], vs[m]
    for i in range(len(xs) - 1):
        c = UP_C if vs[i + 1] > vs[i] else DOWN_C
        ax.plot(xs[i:i + 2], vs[i:i + 2], color=c, alpha=0.8, lw=0.6, zorder=2)
    ax.scatter(xs, vs, s=4, color="#8C8C8C", lw=0, zorder=3)

# Shared geometry (inches): identical plot heights and margins so panels align.
AX_H, BOTTOM, TOP, LEFT, RIGHT = 1.55, 0.42, 0.34, 0.42, 0.08
AX_W = {"b": 0.70, "c": 1.95, "d": 1.55, "e": 1.00}


def new_panel(key, extra_right=0.0):
    w = LEFT + AX_W[key] + RIGHT + extra_right
    h = BOTTOM + AX_H + TOP
    fig = plt.figure(figsize=(w, h))
    ax = fig.add_axes([LEFT / w, BOTTOM / h, AX_W[key] / w, AX_H / h])
    return fig, ax


def save(fig, name):
    fig.savefig(OUT(name + ".svg"), format="svg")
    fig.savefig(OUT(name + ".png"), format="png")
    plt.close(fig)


def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def bracket(ax, x1, x2, y, h, text):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.6, color="k", clip_on=False)
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom", fontsize=6)


# ---------- Load ----------
cc = pd.read_csv(CC_CSV, index_col=0).loc[["SecA", "SecB"], ["G1", "S", "G2M"]] * 100

epi = load_epithelial(columns=["patient_id", "sample_id", "treatment_status", "anatomic_site",
                               "metastatic_site", "state", "Factor_2"])
for c in ("patient_id", "sample_id", "treatment_status", "anatomic_site", "metastatic_site", "state"):
    epi[c] = epi[c].astype(str)
epi["site"] = epi["metastatic_site"]

# c treatment bars: primary-tumour epithelium (same site column as the site bars), pre vs post-chemotherapy, pooled.
TRT = ["pre-treatment", "post-chemotherapy"]
ptx = epi[(epi.site == "primary") & epi.treatment_status.isin(TRT)]
trt = pd.crosstab(ptx.treatment_status, ptx.state, normalize="index").reindex(TRT) * 100

SITES = ["primary", "ascites", "metastasis"]
pre_sites = epi[(epi.treatment_status == "pre-treatment") & epi.site.isin(SITES)]
site = pd.crosstab(pre_sites.site, pre_sites.state, normalize="index") * 100
gs = (pre_sites.groupby(["patient_id", "site"], observed=True)
      .agg(n=("state", "size"), secb=("state", lambda x: (x == "SecB").mean()))
      .reset_index())
gs = gs[gs.n >= MIN_CELLS_PAIRED]
multi = gs.groupby("patient_id").site.nunique()
gs = gs[gs.patient_id.isin(multi[multi >= 2].index)]
asc_pts = gs.loc[gs.site == "ascites", "patient_id"].unique()
pw = (gs[gs.patient_id.isin(asc_pts)]
      .pivot(index="patient_id", columns="site", values="secb").reindex(columns=SITES) * 100)

# ---------- Statistics ----------
site_pairs = [("primary", "ascites"), ("ascites", "metastasis"), ("primary", "metastasis")]
rows = []
for a, b in site_pairs:
    d = pw[[a, b]].dropna()
    rows.append({"panel": "d", "comparison": f"{a} vs {b}", "n": len(d),
                 "wilcoxon_p": stats.wilcoxon(d[a], d[b]).pvalue})
for s_ in SITES:
    rows.append({"panel": "c", "comparison": f"{s_} pooled composition (pre-treatment)",
                 "n": int((pre_sites.site == s_).sum()),
                 "n_samples": int(pre_sites.loc[pre_sites.site == s_, "sample_id"].nunique()),
                 "n_patients": int(pre_sites.loc[pre_sites.site == s_, "patient_id"].nunique()),
                 **{f"pct_{k}": round(float(site.loc[s_, k]), 2) for k in STACK}})
    rows.append({"panel": "d", "comparison": f"median SecB % at {s_}",
                 "n": int(pw[s_].notna().sum()), "median": round(float(pw[s_].median()), 2)})
for t_ in TRT:
    sub = ptx[ptx.treatment_status == t_]
    rows.append({"panel": "c", "comparison": f"primary {t_} pooled composition",
                 "n": len(sub), "n_samples": int(sub.sample_id.nunique()),
                 "n_patients": int(sub.patient_id.nunique()),
                 **{f"pct_{k}": round(float(trt.loc[t_, k]), 2) for k in STACK}})
for s_ in ("SecA", "SecB"):
    rows.append({"panel": "b", "comparison": f"{s_} cell-cycle phase (% of cells)",
                 "pct_G1": round(float(cc.loc[s_, "G1"]), 2), "pct_S": round(float(cc.loc[s_, "S"]), 2),
                 "pct_G2M": round(float(cc.loc[s_, "G2M"]), 2)})
st = pd.DataFrame(rows)
st.to_csv(OUT("row2_stats.csv"), index=False)
print(f"ascites-matched patients = {len(pw)}")
p_site = dict(zip(st.comparison, st.wilcoxon_p))

# ---------- b: cell cycle ----------
fig, ax = new_panel("b")
for i, s in enumerate(cc.index):
    bottom = 0
    for j, ph in enumerate(["G1", "S", "G2M"]):
        v = cc.loc[s, ph]
        ax.bar(i, v, bottom=bottom, width=0.72, color=CC_COLS[s][j], lw=0)
        ax.text(i, bottom + v / 2, ph, ha="center", va="center", fontsize=5.5,
                fontweight="bold", color="white" if j == 2 else "black")
        bottom += v
ax.set_xticks([0, 1]); ax.set_xticklabels(cc.index); ax.set_xlim(-0.6, 1.6)
ax.set_ylim(0, 100); ax.set_ylabel("% of cells")
save(fig, "row2_b_cell_cycle")

# ---------- c: merged composition ----------
fig, ax = new_panel("c")
bars = [("Prim", site.loc["primary"]), ("Asc", site.loc["ascites"]),
        ("Met", site.loc["metastasis"]), (None, None),
        ("Pre", trt.loc["pre-treatment"]), ("Post-\nchemo", trt.loc["post-chemotherapy"])]
xt, xl = [], []
for i, (lab, r) in enumerate(bars):
    if r is None:
        continue
    bottom = 0
    for s in STACK:
        v = r[s]
        ax.bar(i, v, bottom=bottom, width=0.75, color=PALETTE[s], lw=0.3, edgecolor="white")
        if v > 9:
            ax.text(i, bottom + v / 2, f"{v:.0f}", ha="center", va="center", fontsize=5.5,
                    color="black" if s == "SecA" else "white")
        bottom += v
    xt.append(i); xl.append(lab)
ax.set_xticks(xt); ax.set_xticklabels(xl); ax.set_ylim(0, 100)
ax.set_ylabel("% of epithelium")
for x0, x1, t in [(0, 2, "Treatment-naïve"), (4, 5, "Primary tumours")]:
    ax.plot([x0 - 0.35, x1 + 0.35], [105, 105], color="k", lw=0.5, clip_on=False)
    ax.text((x0 + x1) / 2, 107, t, ha="center", va="bottom", fontsize=6)
ax.legend(handles=[Patch(color=PALETTE[s], label=s) for s in STACK], ncol=4, frameon=False,
          loc="upper center", bbox_to_anchor=(0.5, -0.17), handlelength=0.9,
          columnspacing=0.7, handletextpad=0.4)
save(fig, "row2_c_composition")

# ---------- d: patient-matched sites ----------
fig, ax = new_panel("d")
x = np.arange(3)
for _, r in pw.iterrows():
    v = r.values; m = ~np.isnan(v)
    patient_line(ax, x, v)
med = pw.median().values
ax.plot(x, med, color="k", lw=1.2, marker="o", ms=4, zorder=5)
for k in range(3):
    ax.text(x[k] + 0.12, med[k], f"{med[k]:.0f}%", fontsize=6, va="center", fontweight="bold")
bracket(ax, 0, 1, 101, 3, stars(p_site["primary vs ascites"]))
bracket(ax, 1, 2, 101, 3, stars(p_site["ascites vs metastasis"]))
bracket(ax, 0, 2, 111, 3, stars(p_site["primary vs metastasis"]))
ax.set_xticks(x); ax.set_xticklabels(["Prim", "Asc", "Met"]); ax.set_xlim(-0.35, 2.55)
ax.set_ylim(0, 118); ax.set_yticks(range(0, 101, 20))
ax.set_ylabel("SecB (% of epithelium)")
ax.text(1, 122, f"Patient-matched (n = {len(pw)})", ha="center", va="bottom", fontsize=6)
save(fig, "row2_d_paired_sites")

print(st.to_string(index=False))

# ---------- e: patient-matched pre/post chemotherapy, SecB % of the epithelium ----------
# Patients with >= 100 epithelial cells at both time points (samples pooled per time point; not
# site-matched); y = % of the patient's epithelial cells (incl. ciliated) called SecB.
tx = epi[epi.treatment_status.isin(["pre-treatment", "post-chemotherapy"])]
gp = (tx.groupby(["patient_id", "treatment_status"], observed=True)["state"]
      .agg(n="size", secb=lambda x: 100 * (x.astype(str) == "SecB").mean()).reset_index())
gp = gp[gp.n >= MIN_CELLS_PAIRED]
prp = (gp.pivot_table(index="patient_id", columns="treatment_status", values="secb")
       .dropna()[["pre-treatment", "post-chemotherapy"]])
p_var = stats.wilcoxon(prp.iloc[:, 0], prp.iloc[:, 1]).pvalue
n_up = int((prp.iloc[:, 1] > prp.iloc[:, 0]).sum())
fig, ax = new_panel("e")
for _, r in prp.iterrows():
    patient_line(ax, [0, 1], r.values)
medv = prp.median().values
ax.plot([0, 1], medv, color="k", lw=1.2, marker="o", ms=4, zorder=5)
ymax = float(np.ceil(prp.values.max() / 10) * 10)
bracket(ax, 0, 1, ymax * 1.02, ymax * 0.03, f"p = {p_var:.2f}")  # exact p (author 2026-10-05)
ax.set_xticks([0, 1]); ax.set_xticklabels(["Pre", "Post-\nchemo"]); ax.set_xlim(-0.4, 1.4)
ax.set_ylim(0, ymax * 1.12)
ax.set_ylabel("SecB (% of epithelium)")
ax.text(0.5, ymax * 1.16, f"Patient-matched (n = {len(prp)})", ha="center", va="bottom", fontsize=6)
save(fig, "row2_e_paired_treatment_secb_pct")
pd.DataFrame([{"panel": "e_variant_secb_pct", "comparison": "pre vs post-chemotherapy", "n": len(prp),
               "wilcoxon_p": p_var, "n_increased": n_up,
               "median_pre": round(float(medv[0]), 2), "median_post": round(float(medv[1]), 2)}]
             ).to_csv(OUT("row2_e_variant_secb_pct_stats.csv"), index=False)
# the panel e test also goes into row2_stats.csv as panel "e"
st = pd.concat([st, pd.DataFrame([{"panel": "e", "comparison": "pre vs post-chemotherapy (SecB % of epithelium)",
                                   "n": len(prp), "wilcoxon_p": p_var, "n_increased": n_up,
                                   "median": round(float(medv[1]), 2)}])], ignore_index=True)
st.to_csv(OUT("row2_stats.csv"), index=False)
print(f"e variant (SecB % of epithelium): n = {len(prp)}, {n_up} increased, "
      f"median {medv[0]:.1f}% -> {medv[1]:.1f}%, Wilcoxon p = {p_var:.4f}")
