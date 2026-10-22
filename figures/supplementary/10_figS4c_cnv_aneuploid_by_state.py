#!/usr/bin/env python3
"""
Figure S4c — CopyKAT aneuploid fraction of the epithelial compartment, with the reference lineages as negative comparators
=========================================================================================================================
PURPOSE
    CopyKAT was run per sample on the epithelial cells plus a fibroblast / T-NK / endothelial reference that was
    supplied as the known-diploid baseline (atlas/05_cnv/as_run/19_cnv_extract.py, 19_cnv_copykat.R). The panel shows
    the per-sample % of CopyKAT-aneuploid cells for the epithelial compartment and the three reference lineages
    (negative comparators). The epithelial states (SecA / Intermediate / SecB / Ciliated) are reported in the table only.

INPUTS (read-only)
    data_root/2026_final_atlas/output/19_cnv/per_sample/*/copykat_prediction.csv  (251 samples; aneuploid/diploid calls;
        CopyKAT is NOT re-run)
    figures._prep.atlas_final_tables.load_obs()        final-atlas celltype_level1, study, sample_id (all cells)
    figures._prep.atlas_final_tables.load_epithelial() epithelial `state` (SecA / Intermediate / SecB / Ciliated)
OUTPUTS
    output_root/figures/supplementary/SuppFig4c_atlas_cnv_aneuploid_epithelial_vs_reference.{svg,png}
    output_root/figures/supplementary/tables/SuppFig4c_cnv_aneuploid_by_epithelial_state.csv  (pooled + per-sample summary)
Panel: per-sample % aneuploid (dots; samples with >= 20 cells of that group) with box (median, IQR), pooled % (bar
label) and n cells / n samples.
MANUSCRIPT PANEL(S)
    Figure S4c.
RUNTIME TIER
    fast (251 per-sample CopyKAT tables + the loader caches).
"""
import glob, os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402
from figures._prep.atlas_final_tables import load_obs, load_epithelial  # noqa: E402

PER_SAMPLE = path("data_root", "2026_final_atlas", "output", "19_cnv", "per_sample")
FIG_DIR = path("output_root", "figures", "supplementary")
TAB_DIR = path("output_root", "figures", "supplementary", "tables")
MIN_CELLS = 20
STATES = ["SecA", "Intermediate", "SecB", "Ciliated"]          # reported in the table only
ORDER = ["Epithelial", "Fibroblast", "Endothelial", "T/NK cell"]   # plotted (author 2026-10-01: annotation validation only)
REF = ["Fibroblast", "Endothelial", "T/NK cell"]   # CopyKAT diploid reference (negative comparators)
DISPLAY = {"Epithelial": "Epithelial (malignant)", "SecA": "SecA", "Intermediate": "Intermediate", "SecB": "SecB", "Ciliated": "Ciliated",
           "Fibroblast": "Fibroblast", "Endothelial": "Endothelial", "T/NK cell": "T/NK cell"}
COL = {"Epithelial": "#E6A141", "SecA": "#E6A141", "Intermediate": "#C08E48", "SecB": "#9A7D55", "Ciliated": "#E07850",
       "Fibroblast": "#BDB3A6", "Endothelial": "#7D4E4E", "T/NK cell": "#87CEFA"}
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.labelsize": 5.5, "ytick.labelsize": 6})

# ---- CopyKAT calls (all 251 samples)
fr = []
for d in sorted(glob.glob(f"{PER_SAMPLE}/*")):
    p = os.path.join(d, "copykat_prediction.csv")
    if not os.path.basename(d).startswith("._") and os.path.exists(p):
        fr.append(pd.read_csv(p))
cnv = pd.concat(fr, ignore_index=True)

# ---- final atlas labels + sample (load_obs), epithelial state (load_epithelial)
obs = load_obs(columns=["celltype_level1", "study", "sample_id_deposited"]).rename(
    columns={"sample_id_deposited": "sample_id"}).astype(str)  # CopyKAT runs are keyed by the deposited sample id (Nath 2021 fix, 2026-10-20)
state = load_epithelial(columns=["state"])["state"].astype(str).rename("schema_nmf")
cnv = cnv[cnv["copykat.pred"].isin(["aneuploid", "diploid"])].set_index("cell.names")
d = cnv.join(obs, how="inner").join(state, how="left")
d = d[(d.celltype_level1 == "Epithelial") | d.celltype_level1.isin(REF)].copy()
d["schema_nmf"] = np.where(d.celltype_level1 == "Epithelial", d.schema_nmf, d.celltype_level1)
d = d.dropna(subset=["schema_nmf"])
d["group"] = np.where(d.celltype_level1 == "Epithelial", "Epithelial", d.celltype_level1)
d["aneuploid"] = d["copykat.pred"] == "aneuploid"
d["sample"] = d.study + "|" + d.sample_id
print(f"cells with a CopyKAT call: {len(d):,} across {d['sample'].nunique()} samples; epithelial {int((d.celltype_level1 == 'Epithelial').sum()):,}")

# ---- pooled and per-sample (plotted groups + epithelial states for the table)
def summarize(key, order):
    pooled = d.groupby(key).aneuploid.agg(n_cells="size", pct_pooled=lambda s: 100 * s.mean()).reindex(order)
    ps = d.groupby(["sample", key]).aneuploid.agg(n="size", pct=lambda s: 100 * s.mean()).reset_index()
    ps = ps[ps.n >= MIN_CELLS]
    summ = ps.groupby(key).pct.agg(n_samples="size", median_per_sample="median",
                                   q1=lambda v: v.quantile(.25), q3=lambda v: v.quantile(.75)).reindex(order)
    return pooled, ps.rename(columns={key: "grp"}), summ
pooled, ps, summ = summarize("group", ORDER)
T = pooled.join(summ)
pS, _, sS = summarize("schema_nmf", STATES); TS = pS.join(sS)
TT = pd.concat([T.assign(level="plotted group"), TS.assign(level="epithelial state (table only)")])
TT.index.name = "group"; TT.insert(0, "display", [DISPLAY[s] for s in TT.index])
TT.round(2).to_csv(f"{TAB_DIR}/SuppFig4c_cnv_aneuploid_by_epithelial_state.csv")
print(TT.round(1).to_string())
# ---- plot
fig, ax = plt.subplots(figsize=(4.6, 1.7))
rng = np.random.default_rng(42)
for i, s in enumerate(ORDER):
    v = ps.loc[ps.grp == s, "pct"].to_numpy()
    ax.boxplot([v], positions=[i], vert=False, widths=0.55, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="white", edgecolor="black", lw=0.5), medianprops=dict(color="black", lw=0.8),
               whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    ax.scatter(v, i + rng.uniform(-0.18, 0.18, len(v)), s=3, color=COL[s], alpha=0.75, lw=0, zorder=3)
    ax.text(102, i, f"{pooled.loc[s, 'pct_pooled']:.1f}%  (n = {int(pooled.loc[s, 'n_cells']):,} cells, "
                    f"{int(summ.loc[s, 'n_samples'])} samples)", va="center", ha="left", fontsize=5.2, clip_on=False)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([DISPLAY[s] for s in ORDER])
for lab, s in zip(ax.get_yticklabels(), ORDER): lab.set_color(COL[s]); lab.set_fontweight("bold")
ax.invert_yaxis(); ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel("Aneuploid cells per sample (%)", fontsize=6)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.text(102, -0.95, "pooled", fontsize=5.2, ha="left", va="center", color="#555555", clip_on=False)
ax.axhline(0.5, color="#999999", lw=0.5, ls=(0, (2, 2)))
ax.text(50, 0.62, "CopyKAT diploid reference (negative comparators)", fontsize=5.2, ha="center", va="top", color="#555555", style="italic")
for ext in ("svg", "png"):
    fig.savefig(f"{FIG_DIR}/SuppFig4c_atlas_cnv_aneuploid_epithelial_vs_reference.{ext}", dpi=450,
                bbox_inches="tight", facecolor="white")
print("saved SuppFig4c")
