#!/usr/bin/env python3
"""
Signature-based TCGA survival with epithelial-fraction adjustment
================================================================
HGSC malignant-states atlas backend.

Rather than deconvolving SecA vs SecB as separate cell types (their profiles are
too similar), this scores each TCGA-OV bulk sample with the SecA and SecB NMF
gene signatures (mean z-score; 177 genes each from step 11e) and uses the
BayesPrism epithelial fraction as a covariate to control for TME contamination.
Tests OS + PFS (full + 5-year) for signature scores and the SecA-SecB
polarization score, alone and adjusted for epithelial purity / stage / age /
platinum. This 22d_signature_scores.csv is the KEY survival cache for Fig 7.

INPUTS:
  - output_root/07_deconvolution_survival/CIBERSORTx_Results.txt (BayesPrism fractions)
  - <data_root>/2026_final_atlas/output/11e_nmf_characterization/11e_gene_classification.csv
  - <data_root>/2026_final_atlas/data/cibersort_data_prev/tcga_ecotyper.txt
  - <data_root>/2026_final_atlas/data/cibersort_data_prev/tcga_hla_clinical.csv

OUTPUTS (output_root/07_deconvolution_survival/):
  - 22d_signature_scores.csv  (KEY cache — Fig 7E/F/G)
  - 22d_km_results.csv, 22d_cox_results.csv
  - 22d_km_*.svg/pdf, 22d_cox_forest_*.svg/pdf, 22d_score_distributions.svg/pdf
  - 22d_signature_survival_report.html

MANUSCRIPT PANELS: Fig 7E (KM OS), 7F (KM PFS), 7G (stepwise Cox forest).

RUNTIME TIER: moderate (lifelines Cox/KM).

SEEDING: deterministic (Cox/KM are deterministic; no RNG draws).

Usage:
    python 05_signature_survival.py
"""

import os
import sys
import warnings
import base64
from io import BytesIO
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

warnings.filterwarnings("ignore")

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test

# ── Paths ─────────────────────────────────────────────────────
OUT_DIR   = path("output_root", "07_deconvolution_survival")
FRACTIONS = os.path.join(OUT_DIR, "CIBERSORTx_Results.txt")
GENE_CLS  = path("data_root", "2026_final_atlas", "output", "11e_nmf_characterization",
                 "11e_gene_classification.csv")
TCGA_EXPR = path("data_root", "2026_final_atlas", "data", "cibersort_data_prev",
                 "tcga_ecotyper.txt")
CLINICAL  = path("data_root", "2026_final_atlas", "data", "cibersort_data_prev",
                 "tcga_hla_clinical.csv")

# ── Style ─────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 450,
})
BLUE, ORANGE, GREY = "#0072B2", "#D55E00", "#999999"


# ============================================================================
# HELPERS
# ============================================================================

def fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    buf.close()
    return b64


def save_fig(fig, name):
    fig.savefig(os.path.join(OUT_DIR, f"{name}.svg"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{name}.pdf"), bbox_inches="tight")
    b64 = fig_to_b64(fig)
    plt.close(fig)
    return b64


def plot_km(df, time_col, event_col, group_col, title, ax, colors=None):
    kmf = KaplanMeierFitter()
    groups = sorted(df[group_col].unique())
    for i, g in enumerate(groups):
        m = df[group_col] == g
        c = colors[i] if colors else None
        kmf.fit(df.loc[m, time_col], df.loc[m, event_col], label=f"{g} (n={m.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=c, lw=1.5)
    ax.set_title(title, fontsize=9, weight="bold")
    ax.set_xlabel("Months"); ax.set_ylabel("Survival probability")
    ax.legend(fontsize=6, loc="lower left"); ax.set_ylim(-0.05, 1.05)
    if len(groups) == 2:
        g1, g2 = [df[df[group_col] == g] for g in groups]
        p = logrank_test(g1[time_col], g2[time_col], g1[event_col], g2[event_col]).p_value
    else:
        p = multivariate_logrank_test(df[time_col], df[group_col], df[event_col]).p_value
    ps = f"p = {p:.4f}" if p >= 0.0001 else "p < 0.0001"
    ax.text(0.98, 0.98, f"Log-rank {ps}", transform=ax.transAxes, ha="right", va="top",
            fontsize=7, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
    return p


def clip_surv(df, tcol, ecol, maxm):
    d = df.copy()
    over = d[tcol] > maxm
    d.loc[over, ecol] = 0
    d.loc[over, tcol] = maxm
    return d


def truncate_id(bc):
    return "-".join(str(bc).replace(".", "-").split("-")[:3])


def pval_fmt(p):
    if p < 0.001:
        return "<0.001"
    if p < 0.01:
        return f"{p:.3f}"
    return f"{p:.4f}"


def pval_badge(p):
    if p < 0.05:
        return f'<span class="sig">p={pval_fmt(p)}</span>'
    if p < 0.1:
        return f'<span class="trend">p={pval_fmt(p)}</span>'
    return f'<span class="ns">p={pval_fmt(p)}</span>'


def img(b64):
    return f'<img src="data:image/png;base64,{b64}" />'


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 65)
    print("  Signature survival with epithelial adjustment")
    print("=" * 65)
    html_parts = []

    print("\n[1] Loading gene signatures...")
    gc = pd.read_csv(GENE_CLS)
    secA_genes = gc[gc["class"] == "SecA-specific"]["gene"].tolist()
    secB_genes = gc[gc["class"] == "SecB-specific"]["gene"].tolist()
    secA_weights = dict(zip(gc[gc["class"] == "SecA-specific"]["gene"],
                            gc[gc["class"] == "SecA-specific"]["secA_loading"]))
    secB_weights = dict(zip(gc[gc["class"] == "SecB-specific"]["gene"],
                            gc[gc["class"] == "SecB-specific"]["secB_loading"]))
    print(f"    SecA genes: {len(secA_genes)} | SecB genes: {len(secB_genes)}")

    print("\n[2] Loading TCGA expression...")
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    print(f"    Expression: {expr.shape[0]} genes x {expr.shape[1]} samples")
    secA_in = [g for g in secA_genes if g in expr.index]
    secB_in = [g for g in secB_genes if g in expr.index]
    print(f"    SecA genes in TCGA: {len(secA_in)}/{len(secA_genes)}")
    print(f"    SecB genes in TCGA: {len(secB_in)}/{len(secB_genes)}")

    print("\n[3] Computing signature scores...")
    expr_z = expr.apply(lambda row: (row - row.mean()) / (row.std() + 1e-10), axis=1)
    secA_score = expr_z.loc[secA_in].mean(axis=0)
    secB_score = expr_z.loc[secB_in].mean(axis=0)
    secA_w = np.array([secA_weights[g] for g in secA_in])
    secB_w = np.array([secB_weights[g] for g in secB_in])
    secA_score_w = (expr_z.loc[secA_in].T * secA_w).T.sum(axis=0) / secA_w.sum()
    secB_score_w = (expr_z.loc[secB_in].T * secB_w).T.sum(axis=0) / secB_w.sum()
    scores = pd.DataFrame({
        "secA_score": secA_score, "secB_score": secB_score,
        "secA_score_weighted": secA_score_w, "secB_score_weighted": secB_score_w,
        "polarization": secA_score - secB_score,
        "polarization_weighted": secA_score_w - secB_score_w,
    })
    print(f"    Correlation(SecA, SecB): r={np.corrcoef(secA_score, secB_score)[0,1]:.4f}")

    print("\n[4] Loading fractions & clinical data...")
    frac = pd.read_csv(FRACTIONS, sep="\t", index_col=0)
    frac = frac[[c for c in frac.columns if c not in ["P-value", "Correlation", "RMSE"]]]
    epi_cols = [c for c in frac.columns if "epithelium" in c.lower()]
    frac["epi_fraction"] = frac[epi_cols].sum(axis=1)
    print(f"    Epithelial fraction: mean={frac['epi_fraction'].mean():.3f}")
    clinical = pd.read_csv(CLINICAL, index_col=0)

    frac_t = {truncate_id(s): s for s in frac.index}
    clin_t = {truncate_id(s): s for s in clinical["PATIENT.NUMBER"]}
    score_t = {truncate_id(s): s for s in scores.index}
    common = set(frac_t) & set(clin_t) & set(score_t)
    print(f"    Three-way match: {len(common)} samples")

    rows = []
    for tid in common:
        r = {"sample_id": tid}
        for col in scores.columns:
            r[col] = scores.loc[score_t[tid], col]
        r["epi_fraction"] = frac.loc[frac_t[tid], "epi_fraction"]
        crow = clinical[clinical["PATIENT.NUMBER"] == clin_t[tid]].iloc[0]
        r["os_months"] = pd.to_numeric(crow.get("OS..Months."), errors="coerce")
        r["os_event"] = 1 if "DECEASED" in str(crow.get("OS", "")) else 0
        r["pfs_months"] = pd.to_numeric(crow.get("PFS..Months."), errors="coerce")
        r["pfs_event"] = 1 if "PROGRESSION" in str(crow.get("PFS", "")) else 0
        r["stage_coded"] = pd.to_numeric(crow.get("Stage..Coded."), errors="coerce")
        r["age"] = pd.to_numeric(crow.get("Age"), errors="coerce")
        r["platinum_coded"] = pd.to_numeric(crow.get("Platinum..Coded."), errors="coerce")
        rows.append(r)
    df = pd.DataFrame(rows).set_index("sample_id")
    df = df[~df.index.duplicated(keep="first")]

    n_os = df[["os_months", "os_event"]].dropna().shape[0]
    e_os = int(df.dropna(subset=["os_months", "os_event"])["os_event"].sum())
    n_pfs = df[["pfs_months", "pfs_event"]].dropna().shape[0]
    e_pfs = int(df.dropna(subset=["pfs_months", "pfs_event"])["pfs_event"].sum())
    print(f"    OS: {n_os} patients ({e_os} events) | PFS: {n_pfs} ({e_pfs})")

    for col, label_lo, label_hi in [
        ("secA_score", "SecA-low", "SecA-high"),
        ("secB_score", "SecB-low", "SecB-high"),
        ("polarization", "SecB-like", "SecA-like"),
        ("polarization_weighted", "SecB-like (w)", "SecA-like (w)"),
    ]:
        med = df[col].median()
        df[f"{col}_group"] = pd.cut(df[col], bins=[-np.inf, med, np.inf],
                                    labels=[label_lo, label_hi])
    t1, t2 = df["polarization"].quantile([1/3, 2/3])
    df["polar_tertile"] = pd.cut(df["polarization"], bins=[-np.inf, t1, t2, np.inf],
                                 labels=["SecB-like", "Intermediate", "SecA-like"])

    df.to_csv(os.path.join(OUT_DIR, "22d_signature_scores.csv"))
    print("    Saved: 22d_signature_scores.csv")

    print("\n[5] Score distributions...")
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.5))
    axes[0].hist(df["secA_score"], bins=30, color=BLUE, alpha=0.7, edgecolor="white")
    axes[0].set_title("SecA Signature Score"); axes[0].set_xlabel("Mean z-score")
    axes[1].hist(df["secB_score"], bins=30, color=ORANGE, alpha=0.7, edgecolor="white")
    axes[1].set_title("SecB Signature Score"); axes[1].set_xlabel("Mean z-score")
    axes[2].hist(df["polarization"], bins=30, color=GREY, alpha=0.7, edgecolor="white")
    axes[2].axvline(0, color="black", ls="--", lw=0.8)
    axes[2].set_title("Polarization (SecA - SecB)"); axes[2].set_xlabel("Score")
    axes[3].scatter(df["secA_score"], df["secB_score"], s=8, alpha=0.5,
                    c=df["epi_fraction"], cmap="viridis")
    axes[3].set_xlabel("SecA score"); axes[3].set_ylabel("SecB score")
    axes[3].set_title("SecA vs SecB (color=epi frac)")
    cb = plt.colorbar(axes[3].collections[0], ax=axes[3], shrink=0.8)
    cb.set_label("Epi fraction", fontsize=7)
    for ax in axes:
        ax.tick_params(labelsize=7)
    fig.tight_layout()
    dist_b64 = save_fig(fig, "22d_score_distributions")

    r_secA_epi = np.corrcoef(df["secA_score"].dropna(), df["epi_fraction"].dropna())[0, 1]
    r_secB_epi = np.corrcoef(df["secB_score"].dropna(), df["epi_fraction"].dropna())[0, 1]
    n_seca_like = (df["polarization"] > 0).sum()
    overview_html = f"""
    <p>Each TCGA-OV bulk sample scored with SecA/SecB NMF gene signatures (177 genes each,
    mean z-score). BayesPrism epithelial fraction used as covariate for TME contamination.</p>
    <h3>Score Distributions</h3>{img(dist_b64)}
    <h3>Summary Statistics</h3><table>
    <tr><th>Metric</th><th>Value</th></tr>
    <tr><td>Samples</td><td>{len(df)}</td></tr>
    <tr><td>SecA-like (polarization &gt; 0)</td><td>{n_seca_like}/{len(df)} ({n_seca_like/len(df)*100:.1f}%)</td></tr>
    <tr><td>Corr(SecA, SecB scores)</td><td>r={np.corrcoef(df['secA_score'],df['secB_score'])[0,1]:.4f}</td></tr>
    <tr><td>Corr(SecA score, epi fraction)</td><td>r={r_secA_epi:.4f}</td></tr>
    <tr><td>Corr(SecB score, epi fraction)</td><td>r={r_secB_epi:.4f}</td></tr>
    <tr><td>OS patients (events)</td><td>{n_os} ({e_os})</td></tr>
    <tr><td>PFS patients (events)</td><td>{n_pfs} ({e_pfs})</td></tr>
    </table>"""
    html_parts.append(("Overview & Score Distributions", overview_html))

    all_km, all_cox = [], []
    for max_mo, fu_label, prefix in [(None, "Full follow-up", ""), (60, "5-year (60 mo)", "_5yr")]:
        print(f"\n{'='*55}\n  Survival: {fu_label}\n{'='*55}")
        dfs = clip_surv(clip_surv(df, "os_months", "os_event", max_mo),
                        "pfs_months", "pfs_event", max_mo) if max_mo else df.copy()

        km_imgs = []
        for ep, tcol, ecol in [("OS", "os_months", "os_event"), ("PFS", "pfs_months", "pfs_event")]:
            ep_label = "Overall Survival" if ep == "OS" else "Progression-Free Survival"
            for gcol, split_name, colors in [
                ("secA_score_group", "SecA Score Median", [BLUE, ORANGE]),
                ("secB_score_group", "SecB Score Median", [BLUE, ORANGE]),
                ("polarization_group", "Polarization Median", [ORANGE, BLUE]),
                ("polar_tertile", "Polarization Tertiles", [ORANGE, GREY, BLUE]),
                ("polarization_weighted_group", "Weighted Polarization Median", [ORANGE, BLUE]),
            ]:
                sub = dfs.dropna(subset=[tcol, ecol, gcol])
                if len(sub) < 20:
                    continue
                fig, ax = plt.subplots(1, 1, figsize=(5, 4))
                p = plot_km(sub, tcol, ecol, gcol,
                            f"{ep_label} by {split_name}\n{fu_label}", ax, colors=colors)
                fig.tight_layout()
                tag = gcol.replace("_group", "").replace("_tertile", "_tert")
                b64 = save_fig(fig, f"22d_km_{ep.lower()}_{tag}{prefix}")
                all_km.append({"endpoint": ep, "grouping": split_name, "pvalue": p,
                               "n": len(sub), "followup": fu_label})
                km_imgs.append((f"{ep} {split_name}", p, len(sub), b64))
                star = "***" if p < 0.05 else ("  ~" if p < 0.1 else "   ")
                print(f"  {star} KM {ep} {split_name}: p={p:.4f} (n={len(sub)})")

        cox_this = []
        for ep, tcol, ecol in [("OS", "os_months", "os_event"), ("PFS", "pfs_months", "pfs_event")]:
            print(f"\n    --- Cox {ep} ({fu_label}) ---")
            sub = dfs.dropna(subset=[tcol, ecol]).copy()
            if len(sub) < 20:
                continue
            for var, vname in [("secA_score", "SecA_score"), ("secB_score", "SecB_score"),
                               ("polarization", "Polarization"),
                               ("polarization_weighted", "Polarization_weighted")]:
                d = sub[[tcol, ecol, var]].dropna()
                if len(d) < 20:
                    continue
                cph = CoxPHFitter()
                try:
                    cph.fit(d, tcol, ecol)
                    s = cph.summary
                    cox_this.append({"endpoint": ep, "model": "univariate", "variable": vname,
                                     "HR": round(np.exp(s["coef"].values[0]), 4),
                                     "HR_CI_low": round(np.exp(s["coef lower 95%"].values[0]), 4),
                                     "HR_CI_high": round(np.exp(s["coef upper 95%"].values[0]), 4),
                                     "p_value": s["p"].values[0], "n": len(d),
                                     "events": int(d[ecol].sum()), "followup": fu_label})
                except Exception as e:
                    print(f"      {vname}: failed -- {e}")
            for var, vname in [("polarization", "Polarization"), ("secA_score", "SecA_score"),
                               ("secB_score", "SecB_score")]:
                d = sub[[tcol, ecol, var, "epi_fraction"]].dropna()
                if len(d) < 20:
                    continue
                cph = CoxPHFitter()
                try:
                    cph.fit(d, tcol, ecol)
                    for v in [var, "epi_fraction"]:
                        s = cph.summary.loc[v]
                        label = f"{vname}+epi" if v == var else "epi_fraction"
                        cox_this.append({"endpoint": ep, "model": "epi-adjusted", "variable": label,
                                         "HR": round(np.exp(s["coef"]), 4),
                                         "HR_CI_low": round(np.exp(s["coef lower 95%"]), 4),
                                         "HR_CI_high": round(np.exp(s["coef upper 95%"]), 4),
                                         "p_value": s["p"], "n": len(d),
                                         "events": int(d[ecol].sum()), "followup": fu_label})
                except Exception as e:
                    print(f"      {vname} (epi-adj): failed -- {e}")
            mv_cols = ["polarization", "epi_fraction"]
            for c in ["stage_coded", "age", "platinum_coded"]:
                if c in sub.columns:
                    mv_cols.append(c)
            d = sub[[tcol, ecol] + mv_cols].dropna()
            if len(d) >= 30:
                cph = CoxPHFitter()
                try:
                    cph.fit(d, tcol, ecol)
                    for v in mv_cols:
                        s = cph.summary.loc[v]
                        cox_this.append({"endpoint": ep, "model": "full-multivariate", "variable": v,
                                         "HR": round(np.exp(s["coef"]), 4),
                                         "HR_CI_low": round(np.exp(s["coef lower 95%"]), 4),
                                         "HR_CI_high": round(np.exp(s["coef upper 95%"]), 4),
                                         "p_value": s["p"], "n": len(d),
                                         "events": int(d[ecol].sum()), "followup": fu_label})
                except Exception as e:
                    print(f"      Full MV failed: {e}")
        all_cox.extend(cox_this)

        forest_imgs = {}
        for ep in ["OS", "PFS"]:
            ep_cox = [r for r in cox_this if r["endpoint"] == ep]
            if not ep_cox:
                continue
            fig, ax = plt.subplots(1, 1, figsize=(8, max(3, len(ep_cox) * 0.4)))
            for i, r in enumerate(reversed(ep_cox)):
                hr, lo, hi = r["HR"], r["HR_CI_low"], r["HR_CI_high"]
                c = ORANGE if r["p_value"] < 0.05 else ("#F0AD4E" if r["p_value"] < 0.1 else GREY)
                ax.plot([lo, hi], [i, i], color=c, lw=2)
                ax.plot(hr, i, "o", color=c, ms=7)
                ax.text(-0.05, i, f"{r['variable']} ({r['model']})",
                        transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6)
                ps = f"p={r['p_value']:.4f}" if r['p_value'] >= 0.0001 else "p<0.0001"
                ax.text(hi + 0.02, i, f"HR={hr:.2f} ({lo:.2f}-{hi:.2f}) {ps}", va="center", fontsize=5.5)
            ax.axvline(1.0, color="black", ls="--", lw=0.8)
            ax.set_xlabel("Hazard Ratio")
            ax.set_title(f"Cox Forest -- {ep} ({fu_label})", fontsize=9, weight="bold")
            ax.set_yticks([]); ax.set_xlim(left=0)
            fig.tight_layout()
            forest_imgs[ep] = save_fig(fig, f"22d_cox_forest_{ep.lower()}{prefix}")

        sec_html = ""
        if max_mo:
            sec_html += f"<p>Survival censored at <b>{max_mo} months</b>.</p>\n"
        sec_html += "<h3>Kaplan-Meier</h3>\n<table>"
        sec_html += "<tr><th>Endpoint</th><th>Grouping</th><th>n</th><th>p-value</th></tr>\n"
        for name, p, n, _ in km_imgs:
            sec_html += (f"<tr><td>{name.split()[0]}</td><td>{' '.join(name.split()[1:])}</td>"
                         f"<td>{n}</td><td>{pval_badge(p)}</td></tr>\n")
        sec_html += "</table>\n<div class=\"fig-row\">\n"
        for _, _, _, b64 in km_imgs:
            sec_html += img(b64) + "\n"
        sec_html += "</div>\n"
        if cox_this:
            sec_html += "<h3>Cox PH</h3>\n<table>"
            sec_html += ("<tr><th>Endpoint</th><th>Model</th><th>Variable</th><th>HR</th>"
                         "<th>95% CI</th><th>p</th><th>n</th></tr>\n")
            for r in cox_this:
                sec_html += (f"<tr><td>{r['endpoint']}</td><td>{r['model']}</td>"
                             f"<td>{r['variable']}</td><td>{r['HR']:.3f}</td>"
                             f"<td>{r['HR_CI_low']:.3f}-{r['HR_CI_high']:.3f}</td>"
                             f"<td>{pval_badge(r['p_value'])}</td><td>{r['n']}</td></tr>\n")
            sec_html += "</table>\n"
        if forest_imgs:
            sec_html += '<h3>Forest Plots</h3>\n<div class="fig-row">\n'
            for b64 in forest_imgs.values():
                sec_html += img(b64) + "\n"
            sec_html += "</div>\n"
        html_parts.append((f"Survival -- {fu_label}", sec_html))

    if all_km:
        pd.DataFrame(all_km).to_csv(os.path.join(OUT_DIR, "22d_km_results.csv"), index=False)
    if all_cox:
        pd.DataFrame(all_cox).to_csv(os.path.join(OUT_DIR, "22d_cox_results.csv"), index=False)

    sig_km = [r for r in all_km if r["pvalue"] < 0.05]
    sig_cox = [r for r in all_cox if r["p_value"] < 0.05]
    trend_km = [r for r in all_km if 0.05 <= r["pvalue"] < 0.1]
    trend_cox = [r for r in all_cox if 0.05 <= r["p_value"] < 0.1]
    findings = '<div class="key-finding">\n<h3>Summary</h3>\n<ul>\n'
    findings += f"<li>KM significant: <b>{len(sig_km)}/{len(all_km)}</b></li>\n"
    findings += f"<li>Cox significant: <b>{len(sig_cox)}/{len(all_cox)}</b></li>\n"
    findings += f"<li>KM trends: <b>{len(trend_km)}</b> | Cox trends: <b>{len(trend_cox)}</b></li>\n</ul>\n"
    if sig_km or sig_cox:
        findings += "<h3>Significant Results</h3>\n<ul>\n"
        for r in sig_km:
            findings += (f"<li><b>KM {r['followup']}</b> {r['endpoint']} "
                         f"{r['grouping']}: p={r['pvalue']:.4f}</li>\n")
        for r in sig_cox:
            findings += (f"<li><b>Cox {r['followup']}</b> {r['endpoint']} {r['variable']} "
                         f"({r['model']}): HR={r['HR']:.3f}, p={r['p_value']:.4f}</li>\n")
        findings += "</ul>\n"
    findings += "</div>\n"
    html_parts.append(("Key Findings", findings))

    toc = "".join(f'<li><a href="#s{i}">{t}</a></li>\n' for i, (t, _) in enumerate(html_parts))
    body = "".join(f'<div class="section" id="s{i}"><h2>{t}</h2>\n{c}\n</div>\n'
                   for i, (t, c) in enumerate(html_parts))
    html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Signature Survival Report</title>
<style>
  body {{ font-family:Helvetica,Arial,sans-serif; font-size:14px; line-height:1.5;
         max-width:1100px; margin:0 auto; padding:20px 30px; background:#fafafa; color:#333; }}
  h1 {{ font-size:22px; border-bottom:3px solid #0072B2; padding-bottom:8px; color:#0072B2; }}
  .subtitle {{ color:#666; font-size:13px; margin-bottom:20px; }}
  h2 {{ font-size:17px; margin:30px 0 12px; border-bottom:1px solid #ddd; padding-bottom:4px; }}
  h3 {{ font-size:14px; margin:18px 0 8px; color:#555; }}
  .section {{ background:#fff; border:1px solid #e0e0e0; border-radius:6px;
              padding:18px 22px; margin-bottom:18px; }}
  table {{ border-collapse:collapse; width:100%; margin:10px 0; font-size:13px; }}
  th,td {{ border:1px solid #ddd; padding:6px 10px; text-align:left; }}
  th {{ background:#f5f5f5; font-weight:600; }}
  .sig {{ color:#D55E00; font-weight:700; }} .trend {{ color:#F0AD4E; font-weight:600; }} .ns {{ color:#999; }}
  img {{ max-width:100%; height:auto; margin:8px 0; border:1px solid #eee; border-radius:4px; }}
  .fig-row {{ display:flex; flex-wrap:wrap; gap:12px; justify-content:center; }}
  .fig-row img {{ max-width:48%; }}
  .toc {{ background:#fff; border:1px solid #e0e0e0; border-radius:6px; padding:12px 20px; margin-bottom:18px; }}
  .toc ul {{ list-style:none; padding-left:0; columns:2; }} .toc a {{ color:#0072B2; text-decoration:none; }}
  .key-finding {{ background:#FFF8E1; border-left:4px solid #F0AD4E; padding:10px 14px; margin:10px 0; border-radius:0 4px 4px 0; }}
</style></head><body>
<h1>Signature-Based Survival Analysis</h1>
<p class="subtitle">HGSC malignant-states atlas &middot; {datetime.now().strftime('%Y-%m-%d %H:%M')}<br>
SecA/SecB gene signature scoring with BayesPrism epithelial fraction adjustment</p>
<div class="toc"><strong>Contents</strong><ul>{toc}</ul></div>
{body}
</body></html>"""
    with open(os.path.join(OUT_DIR, "22d_signature_survival_report.html"), "w") as f:
        f.write(html)
    print("\n    Saved: 22d_signature_survival_report.html")
    print(f"\n{'='*65}\n  Step complete! Output: {OUT_DIR}\n{'='*65}")


if __name__ == "__main__":
    main()
