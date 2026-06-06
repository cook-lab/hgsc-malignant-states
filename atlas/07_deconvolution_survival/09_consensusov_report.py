#!/usr/bin/env python3
"""
consensusOV TCGA subtype x epitype distribution report
======================================================
HGSC malignant-states atlas backend.

Joins the per-sample consensusOV TCGA subtype calls (08_consensusov_score.R)
with the NMF epitype proportions (07_consensusov_export.py) and asks how
SecA / Intermediate / SecB epithelial proportions distribute across the four
TCGA HGSC subtypes (Differentiated / Immunoreactive / Mesenchymal /
Proliferative) for the primary-untreated cohort. Emits stacked-bar, violin,
per-sample heatmap, and bulk-vs-epi concordance figures + an HTML report.

INPUTS (output_root/07_deconvolution_survival/consensusov/):
  - pseudobulk_metadata.csv
  - consensusov_calls_bulk.csv, consensusov_calls_epi.csv

OUTPUTS (output_root/07_deconvolution_survival/consensusov/):
  - figs/20c_*.svg/pdf
  - tables/20c_per_sample_joined.csv  (KEY cache — Fig 3H), 20c_subtype_composition_summary.csv,
    20c_kruskal_dunn_stats.csv, 20c_bulk_vs_epi_concordance.csv
  - 20_consensusov_report.html

MANUSCRIPT PANELS: Fig 3H (epitype x TCGA subtype dot matrix).

RUNTIME TIER: fast (table join + plotting).

SEEDING: violin strip jitter uses np.random.default_rng(SEED) for determinism.

Usage:
    python 09_consensusov_report.py
"""

import os
import sys
import io
import base64
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.stats as st

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ── Cook Lab style v1.2 ─────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":         8,
    "axes.titlesize":    9,
    "axes.labelsize":    8,
    "xtick.labelsize":   7,
    "ytick.labelsize":   7,
    "legend.fontsize":   6,
    "figure.dpi":        450,
    "savefig.dpi":       450,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "svg.fonttype":      "none",
    "savefig.bbox":      "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# ── Paths ───────────────────────────────────────────────────
DATA_DIR = path("output_root", "07_deconvolution_survival", "consensusov")
FIG_DIR  = os.path.join(DATA_DIR, "figs")
TBL_DIR  = os.path.join(DATA_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TBL_DIR, exist_ok=True)

# ── Palettes ("Transitioning" -> "Intermediate") ────────────
EPITYPE_PALETTE = {
    "SecA":         "#F6D28B",
    "Intermediate": "#7D4E4E",
    "SecB":         "#B8741A",
    "Ciliated":     "#E05A2C",
    "Other_epi":    "#D9C5A2",
}
EPITYPE_ORDER = ["SecA", "Intermediate", "SecB", "Ciliated", "Other_epi"]

SUBTYPE_ORDER = ["DIF", "IMR", "MES", "PRO"]
SUBTYPE_FULL = {"DIF": "Differentiated", "IMR": "Immunoreactive",
                "MES": "Mesenchymal", "PRO": "Proliferative"}
SUBTYPE_PALETTE = {"DIF": "#4F81BD", "IMR": "#9BBB59", "MES": "#C0504D", "PRO": "#8064A2"}


def main():
    print("=" * 70)
    print("  consensusOV TCGA subtype x epitype report")
    print("=" * 70)

    meta = pd.read_csv(os.path.join(DATA_DIR, "pseudobulk_metadata.csv"))
    bulk = pd.read_csv(os.path.join(DATA_DIR, "consensusov_calls_bulk.csv"))
    epi  = pd.read_csv(os.path.join(DATA_DIR, "consensusov_calls_epi.csv"))

    bulk["subtype_bulk"] = bulk["consensusOV"].apply(_canon_subtype)
    epi["subtype_epi"]   = epi["consensusOV"].apply(_canon_subtype)

    bulk_use = bulk[["sample_id", "subtype_bulk", "margin_top1_top2"]].rename(
        columns={"margin_top1_top2": "margin_bulk"})
    epi_use = epi[["sample_id", "subtype_epi", "margin_top1_top2"]].rename(
        columns={"margin_top1_top2": "margin_epi"})

    df = meta.merge(bulk_use, on="sample_id", how="left").merge(epi_use, on="sample_id", how="left")
    df.to_csv(os.path.join(TBL_DIR, "20c_per_sample_joined.csv"), index=False)
    print(f"\nMerged table: {df.shape[0]} samples")

    fig_stacked = _stacked_epitype_by_subtype(df, view="bulk")
    _save(fig_stacked, "20c_stacked_epitype_by_subtype")

    figs_violins, stats_rows = {}, []
    for col, lab in [("pct_secA", "SecA"), ("pct_intermediate", "Intermediate"),
                     ("pct_secB", "SecB")]:
        f, kw, dunn = _violin_by_subtype(df, col, lab, view="bulk")
        figs_violins[col] = f
        _save(f, f"20c_violin_{col}_by_subtype")
        stats_rows.append({"epitype": lab, "metric": col,
                           "kruskal_H": kw.statistic if kw is not None else np.nan,
                           "kruskal_pvalue": kw.pvalue if kw is not None else np.nan})
        if dunn is not None:
            for (a, b), p in dunn.items():
                stats_rows.append({"epitype": lab, "metric": col, "kruskal_H": np.nan,
                                   "kruskal_pvalue": np.nan, "dunn_pair": f"{a} vs {b}",
                                   "dunn_pvalue": p})
    pd.DataFrame(stats_rows).to_csv(os.path.join(TBL_DIR, "20c_kruskal_dunn_stats.csv"), index=False)

    fig_heatmap = _per_sample_heatmap(df, view="bulk")
    _save(fig_heatmap, "20c_per_sample_heatmap")

    fig_conf, conc_df = _bulk_vs_epi_confusion(df)
    _save(fig_conf, "20c_bulk_vs_epi_confusion")
    conc_df.to_csv(os.path.join(TBL_DIR, "20c_bulk_vs_epi_concordance.csv"), index=False)

    summary = _composition_summary(df, view="bulk")
    summary.to_csv(os.path.join(TBL_DIR, "20c_subtype_composition_summary.csv"), index=False)

    _write_html_report(df, summary, conc_df)
    print("\n[done]")
    print(f"  Figures:  {FIG_DIR}/")
    print(f"  Tables:   {TBL_DIR}/")


def _canon_subtype(x):
    if pd.isna(x):
        return np.nan
    return str(x).split("_")[0].upper()


def _save(fig, stem):
    for ext in ("svg", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{stem}.{ext}"))
    plt.close(fig)


def _stacked_epitype_by_subtype(df, view="bulk"):
    col = f"subtype_{view}"
    grp = df.dropna(subset=[col])
    means = grp.groupby(col)[["pct_secA", "pct_intermediate", "pct_secB",
                              "pct_ciliated", "pct_other_epi"]].mean() \
        .reindex([s for s in SUBTYPE_ORDER if s in grp[col].unique()])
    means.columns = ["SecA", "Intermediate", "SecB", "Ciliated", "Other_epi"]
    n_per = grp[col].value_counts().reindex(means.index).astype(int)

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    bottom = np.zeros(len(means))
    for ep in EPITYPE_ORDER:
        ax.bar(range(len(means)), means[ep], bottom=bottom, color=EPITYPE_PALETTE[ep],
               edgecolor="white", linewidth=0.4, label=ep)
        bottom += means[ep].to_numpy()
    ax.set_xticks(range(len(means)))
    ax.set_xticklabels([f"{SUBTYPE_FULL.get(s, s)}\n(n={n_per[s]})" for s in means.index], rotation=0)
    ax.set_ylabel("Mean % of epithelial cells"); ax.set_ylim(0, 100)
    ax.set_title(f"Mean epitype composition by TCGA subtype  ·  {view} pseudobulk")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    return fig


def _violin_by_subtype(df, col, label, view="bulk"):
    s_col = f"subtype_{view}"
    sub = df.dropna(subset=[s_col, col])
    groups = [s for s in SUBTYPE_ORDER if s in sub[s_col].unique()]
    data = [sub.loc[sub[s_col] == s, col].to_numpy() for s in groups]
    kw = st.kruskal(*data) if (len(groups) >= 2 and all(len(d) >= 2 for d in data)) else None
    dunn = _dunn_test(sub, s_col, col, groups) if (kw is not None and kw.pvalue < 0.10
                                                   and len(groups) >= 2) else None

    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    parts = ax.violinplot(data, positions=range(len(groups)),
                          showmeans=False, showmedians=True, widths=0.8)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(SUBTYPE_PALETTE.get(groups[i], "#888"))
        pc.set_alpha(0.6); pc.set_edgecolor("black"); pc.set_linewidth(0.4)
    for k in ("cmedians", "cmaxes", "cmins", "cbars"):
        if k in parts:
            parts[k].set_color("black"); parts[k].set_linewidth(0.6)
    rng = np.random.default_rng(SEED)
    for i, d in enumerate(data):
        ax.scatter(i + rng.uniform(-0.12, 0.12, size=len(d)), d,
                   s=6, color="black", alpha=0.55, linewidths=0)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{SUBTYPE_FULL.get(s, s)}\n(n={len(d)})" for s, d in zip(groups, data)],
                       rotation=0)
    ax.set_ylabel(f"% {label} (per sample)")
    title = f"% {label} per sample by TCGA subtype  ·  {view} pseudobulk"
    if kw is not None:
        title += f"\nKruskal-Wallis H={kw.statistic:.2f}, p={kw.pvalue:.2g}"
    ax.set_title(title)
    fig.tight_layout()
    return fig, kw, dunn


def _dunn_test(df, group_col, value_col, groups):
    df2 = df[[group_col, value_col]].dropna().copy()
    df2["rank"] = df2[value_col].rank(method="average")
    n = len(df2)
    mean_ranks = df2.groupby(group_col)["rank"].mean()
    sizes = df2.groupby(group_col).size()
    pairs = list(combinations(groups, 2))
    raw = []
    for a, b in pairs:
        if a not in mean_ranks or b not in mean_ranks:
            raw.append(np.nan); continue
        se = np.sqrt(n * (n + 1) / 12.0 * (1.0 / sizes[a] + 1.0 / sizes[b]))
        z = (mean_ranks[a] - mean_ranks[b]) / se if se > 0 else 0.0
        raw.append(2 * (1 - st.norm.cdf(abs(z))))
    raw_arr = np.array(raw, dtype=float)
    valid = ~np.isnan(raw_arr)
    pvals = {}
    if valid.sum() > 0:
        adj = np.full_like(raw_arr, np.nan)
        order = np.argsort(raw_arr[valid])
        ranks = np.empty(valid.sum(), dtype=int)
        ranks[order] = np.arange(1, valid.sum() + 1)
        m = valid.sum()
        adj_valid = raw_arr[valid] * m / ranks
        adj_sorted = np.minimum.accumulate(adj_valid[order][::-1])[::-1]
        adj_valid_out = np.empty_like(adj_valid)
        adj_valid_out[order] = adj_sorted
        adj[valid] = np.clip(adj_valid_out, 0, 1)
        for pair, p in zip(pairs, adj):
            pvals[pair] = p
    return pvals


def _per_sample_heatmap(df, view="bulk"):
    s_col = f"subtype_{view}"
    sub = df.dropna(subset=[s_col]).copy().sort_values([s_col, "pct_secA"], ascending=[True, False])
    cols = ["pct_secA", "pct_intermediate", "pct_secB", "pct_ciliated"]
    mat = sub[cols].to_numpy()
    n_samples = mat.shape[0]
    fig, (ax_anno, ax_hm) = plt.subplots(
        nrows=2, ncols=1, sharex=True, figsize=(max(6, n_samples * 0.06), 3.4),
        gridspec_kw=dict(height_ratios=[0.4, 5.0], hspace=0.05))
    sub_colors = sub[s_col].map(SUBTYPE_PALETTE).fillna("#cccccc").to_list()
    ax_anno.imshow(np.array([[matplotlib.colors.to_rgb(c) for c in sub_colors]]), aspect="auto")
    ax_anno.set_yticks([0]); ax_anno.set_yticklabels(["Subtype"]); ax_anno.set_xticks([])
    for sp in ax_anno.spines.values():
        sp.set_visible(False)
    im = ax_hm.imshow(mat.T, aspect="auto", cmap="magma_r", vmin=0, vmax=100)
    ax_hm.set_yticks(range(len(cols)))
    ax_hm.set_yticklabels(["SecA", "Intermediate", "SecB", "Ciliated"])
    ax_hm.set_xticks([]); ax_hm.set_xlabel(f"Samples (n={n_samples}), sorted by subtype")
    cbar = fig.colorbar(im, ax=ax_hm, fraction=0.025, pad=0.02)
    cbar.set_label("% of epithelial cells")
    handles = [plt.Rectangle((0, 0), 1, 1, color=SUBTYPE_PALETTE[s])
               for s in SUBTYPE_ORDER if s in sub[s_col].unique()]
    labels = [SUBTYPE_FULL[s] for s in SUBTYPE_ORDER if s in sub[s_col].unique()]
    ax_anno.legend(handles, labels, ncol=len(labels), bbox_to_anchor=(0.5, 1.7),
                   loc="lower center", frameon=False, handlelength=0.8, columnspacing=1.0)
    return fig


def _bulk_vs_epi_confusion(df):
    sub = df.dropna(subset=["subtype_bulk", "subtype_epi"]).copy()
    if sub.empty:
        fig, ax = plt.subplots(figsize=(3, 2))
        ax.text(0.5, 0.5, "No paired calls", ha="center", va="center"); ax.axis("off")
        return fig, pd.DataFrame()
    subtypes_present = sorted(set(sub["subtype_bulk"].unique()) | set(sub["subtype_epi"].unique()),
                             key=lambda x: SUBTYPE_ORDER.index(x) if x in SUBTYPE_ORDER else 99)
    cm = pd.crosstab(sub["subtype_bulk"], sub["subtype_epi"]).reindex(
        index=subtypes_present, columns=subtypes_present, fill_value=0)
    n_total = int(cm.values.sum())
    n_match = int(np.trace(cm.values))
    fig, ax = plt.subplots(figsize=(3.6, 3.2))
    im = ax.imshow(cm.values, cmap="Blues")
    ax.set_xticks(range(len(cm))); ax.set_xticklabels(cm.columns)
    ax.set_yticks(range(len(cm))); ax.set_yticklabels(cm.index)
    ax.set_xlabel("Epithelial-only call"); ax.set_ylabel("Bulk-tumor call")
    ax.set_title(f"Bulk vs epi-only concordance: {n_match}/{n_total} "
                 f"({100.0*n_match/max(n_total,1):.1f}%)")
    for i in range(len(cm)):
        for j in range(len(cm)):
            v = cm.values[i, j]
            ax.text(j, i, str(v), ha="center", va="center",
                    color="black" if v < cm.values.max() / 2 else "white", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.04)
    fig.tight_layout()
    cm_long = cm.reset_index().melt(id_vars="subtype_bulk", var_name="subtype_epi", value_name="n")
    return fig, cm_long


def _composition_summary(df, view="bulk"):
    s_col = f"subtype_{view}"
    grp = df.dropna(subset=[s_col]).groupby(s_col)
    out = grp.agg(n_samples=("sample_id", "count"), n_patients=("patient_id", "nunique"),
                  mean_pct_secA=("pct_secA", "mean"), mean_pct_intermediate=("pct_intermediate", "mean"),
                  mean_pct_secB=("pct_secB", "mean"), mean_pct_ciliated=("pct_ciliated", "mean"),
                  median_pct_secA=("pct_secA", "median"), median_pct_secB=("pct_secB", "median"))
    out = out.reset_index().rename(columns={s_col: "subtype"})
    out["subtype_full"] = out["subtype"].map(SUBTYPE_FULL)
    return out


def _write_html_report(df, summary, conc):
    out_path = os.path.join(DATA_DIR, "20_consensusov_report.html")

    def _embed(stem):
        with open(os.path.join(FIG_DIR, f"{stem}.svg"), "rb") as fh:
            return base64.b64encode(fh.read()).decode("ascii")

    imgs = {k: _embed(s) for k, s in [
        ("stacked", "20c_stacked_epitype_by_subtype"),
        ("secA", "20c_violin_pct_secA_by_subtype"),
        ("interm", "20c_violin_pct_intermediate_by_subtype"),
        ("secB", "20c_violin_pct_secB_by_subtype"),
        ("heatmap", "20c_per_sample_heatmap"),
        ("confusion", "20c_bulk_vs_epi_confusion"),
    ]}
    n_samples = len(df)
    n_patients = df["patient_id"].nunique()
    studies = sorted(df["study"].dropna().unique().tolist())

    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>consensusOV TCGA subtype mapping</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;font-size:14px;color:#222;max-width:1100px;
      margin:24px auto;padding:0 20px;}}
h1{{font-size:22px;margin-bottom:4px;}}
h2{{font-size:18px;margin-top:28px;border-bottom:1px solid #ddd;padding-bottom:4px;}}
table{{border-collapse:collapse;margin:10px 0;font-size:13px;}}
th,td{{border:1px solid #ccc;padding:4px 8px;text-align:right;}}
th{{background:#f3f3f3;}}
.fig{{margin:18px 0;text-align:center;}} .fig img{{max-width:100%;height:auto;border:1px solid #eee;}}
.caption{{font-size:12px;color:#555;margin-top:4px;}}
.note{{background:#fbf8e6;border-left:3px solid #c4a800;padding:8px 12px;font-size:13px;margin:12px 0;}}
</style></head><body>
<h1>consensusOV TCGA subtype x NMF epitype</h1>
<p><b>Cohort:</b> primary, untreated samples (anatomic_site = adnexa, treatment_status = pre-treatment)<br>
<b>Samples:</b> {n_samples} &nbsp; <b>Patients:</b> {n_patients} &nbsp;
<b>Studies:</b> {", ".join(studies)}<br>
<b>Thresholds:</b> &ge;500 cells/sample, &ge;100 epithelial cells/sample</p>
<div class="note">TCGA subtypes: <b>DIF</b>=Differentiated, <b>IMR</b>=Immunoreactive,
<b>MES</b>=Mesenchymal, <b>PRO</b>=Proliferative. Bulk = all cells per sample;
Epi-only = epithelial-cell pseudobulk (sensitivity check).</div>
<h2>Subtype composition summary (bulk)</h2>
{summary.to_html(index=False, float_format=lambda x: f"{x:.1f}")}
<h2>Mean epitype composition by TCGA subtype</h2>
<div class="fig"><img src="data:image/svg+xml;base64,{imgs['stacked']}">
<div class="caption">Mean per-sample epitype % within each TCGA subtype call (bulk).</div></div>
<h2>Per-epitype distribution across subtypes</h2>
<div class="fig"><img src="data:image/svg+xml;base64,{imgs['secA']}"><div class="caption">% SecA per sample by subtype.</div></div>
<div class="fig"><img src="data:image/svg+xml;base64,{imgs['interm']}"><div class="caption">% Intermediate per sample by subtype.</div></div>
<div class="fig"><img src="data:image/svg+xml;base64,{imgs['secB']}"><div class="caption">% SecB per sample by subtype.</div></div>
<h2>Per-sample heatmap (sorted by subtype)</h2>
<div class="fig"><img src="data:image/svg+xml;base64,{imgs['heatmap']}">
<div class="caption">Each column = one sample; rows = epitype proportions; top strip = TCGA subtype.</div></div>
<h2>Sensitivity: bulk vs epithelial-only concordance</h2>
<div class="fig"><img src="data:image/svg+xml;base64,{imgs['confusion']}">
<div class="caption">Confusion of consensusOV calls between bulk and epi-only pseudobulks.</div></div>
</body></html>"""
    with open(out_path, "w") as fh:
        fh.write(html)
    print(f"  HTML report: {out_path}")


if __name__ == "__main__":
    main()
