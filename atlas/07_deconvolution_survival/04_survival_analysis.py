#!/usr/bin/env python3
"""
Survival analysis on deconvolved TCGA fractions
===============================================
HGSC malignant-states atlas backend.

Joins BayesPrism-estimated cell-type fractions (02_bayesprism_deconv.R) with
TCGA-OV clinical data and tests whether deconvolved SecA/SecB epithelial
fractions predict OS and PFS after accounting for immune contamination in bulk
RNA-seq. KM (median/tertile/ratio splits) + univariate/multivariate Cox, at
full follow-up and clipped to 5 years. Emits an HTML report.

INPUTS:
  - output_root/07_deconvolution_survival/CIBERSORTx_Results.txt  (BayesPrism fractions)
  - <data_root>/2026_final_atlas/data/cibersort_data_prev/tcga_hla_clinical.csv

OUTPUTS (output_root/07_deconvolution_survival/):
  - 22c_survival_report.html
  - 22c_km_*.svg/pdf, 22c_cox_forest_*.svg/pdf, 22c_fraction_*.svg/pdf
  - 22c_deconv_clinical_joined.csv, 22c_cox_results.csv, 22c_km_logrank_results.csv

MANUSCRIPT PANELS: supporting Fig 7E/F/G (TCGA survival).

RUNTIME TIER: moderate (lifelines Cox/KM).

SEEDING: deterministic (Cox/KM are deterministic; no RNG draws).

Usage:
    python 04_survival_analysis.py
    python 04_survival_analysis.py --fractions path/to/results.txt
"""

import os
import sys
import argparse
import warnings
import base64
from io import BytesIO
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

warnings.filterwarnings("ignore")

try:
    from lifelines import KaplanMeierFitter, CoxPHFitter
    from lifelines.statistics import logrank_test, multivariate_logrank_test
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False
    print("WARNING: lifelines not installed. pip install lifelines")

# ── Paths ─────────────────────────────────────────────────────
OUT_DIR  = path("output_root", "07_deconvolution_survival")
CLINICAL = path("data_root", "2026_final_atlas", "data", "cibersort_data_prev",
                "tcga_hla_clinical.csv")
DEFAULT_FRACTIONS = os.path.join(OUT_DIR, "CIBERSORTx_Results.txt")

# ── Style (Cook Lab v1.2) ─────────────────────────────────────
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "pdf.fonttype":    42,
    "svg.fonttype":    "none",
    "savefig.dpi":     450,
    "figure.dpi":      150,
})

GROUP_COLORS = {
    "SecA_epithelium":          "#0072B2",
    "Intermediate_epithelium":  "#999999",
    "SecB_epithelium":          "#D55E00",
    "Ciliated_epithelium":      "#E05A2C",
    "T_NK":                     "#6FBFDF",
    "B_Plasma":                 "#5665B6",
    "Macrophage":               "#8FBC8F",
    "DC":                       "#2E8B57",
    "Fibroblast_Stromal":       "#DDD5CA",
    "Endothelial":              "#7D4E4E",
    "Mesothelial":              "#A8A298",
    "Other_immune":             "#6B8E23",
}


# ============================================================================
# HELPERS
# ============================================================================

def save_fig(fig, name):
    fig.savefig(os.path.join(OUT_DIR, f"{name}.svg"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{name}.pdf"), bbox_inches="tight")
    b64 = fig_to_base64(fig)
    plt.close(fig)
    return b64


def fig_to_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()
    return b64


def plot_km(df, time_col, event_col, group_col, title, ax, colors=None):
    kmf = KaplanMeierFitter()
    groups = sorted(df[group_col].unique())
    # colors is a {group_label: color} dict so the palette stays bound to the
    # semantic label (SecA/SecB) rather than the alphabetical sort position.
    for i, group in enumerate(groups):
        mask = df[group_col] == group
        color = colors.get(group) if colors else None
        kmf.fit(df.loc[mask, time_col], df.loc[mask, event_col],
                label=f"{group} (n={mask.sum()})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=color, lw=1.5)
    ax.set_title(title, fontsize=9, weight="bold")
    ax.set_xlabel("Months"); ax.set_ylabel("Survival probability")
    ax.legend(fontsize=6, loc="lower left"); ax.set_ylim(-0.05, 1.05)
    if len(groups) == 2:
        g1, g2 = [df[df[group_col] == g] for g in groups]
        pval = logrank_test(g1[time_col], g2[time_col],
                            g1[event_col], g2[event_col]).p_value
    else:
        pval = multivariate_logrank_test(df[time_col], df[group_col],
                                         df[event_col]).p_value
    pval_str = f"p = {pval:.4f}" if pval >= 0.0001 else "p < 0.0001"
    ax.text(0.98, 0.98, f"Log-rank {pval_str}", transform=ax.transAxes,
            ha="right", va="top", fontsize=7,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="gray", alpha=0.8))
    return pval


def clip_survival(df, time_col, event_col, max_months):
    dc = df.copy()
    over = dc[time_col] > max_months
    dc.loc[over, event_col] = 0
    dc.loc[over, time_col] = max_months
    return dc


def run_km_suite(df, time_col, event_col, group_cols, secB_col,
                 prefix, label_prefix, print_prefix=""):
    results = []
    # Colors keyed by group label (matches the pd.cut() labels) so blue/orange
    # track SecA-low/SecB-high semantics regardless of alphabetical sort order.
    analyses = [
        ("secB_group_median", "SecB Median Split",
         {"SecB-low": "#0072B2", "SecB-high": "#D55E00"}),
        ("secB_group_tertile", "SecB Tertiles",
         {"SecB-low": "#0072B2", "SecB-mid": "#999999", "SecB-high": "#D55E00"}),
        ("ratio_group", "SecB/SecA Ratio",
         {"SecA-dominant": "#0072B2", "SecB-dominant": "#D55E00"}),
        ("ratio_ab_group", "SecA/SecB Ratio",
         {"Low SecA:SecB": "#0072B2", "High SecA:SecB": "#D55E00"}),
        ("secA_prop_group", "SecA Epithelial Proportion",
         {"Low SecA prop": "#0072B2", "High SecA prop": "#D55E00"}),
    ]
    ep_label = "Overall Survival" if "os" in time_col else "Progression-Free Survival"
    ep_short = "OS" if "os" in time_col else "PFS"
    for gcol, split_name, colors in analyses:
        df_sub = df.dropna(subset=[time_col, event_col, gcol])
        if len(df_sub) < 20:
            continue
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        pval = plot_km(df_sub, time_col, event_col, gcol,
                       f"{ep_label} by {split_name}\n{label_prefix}", ax, colors=colors)
        fig.tight_layout()
        tag = gcol.replace("secB_group_", "secB_").replace("ratio_group", "ratio")
        b64 = save_fig(fig, f"22c_km_{ep_short.lower()}_{tag}{prefix}")
        results.append(({
            "endpoint": ep_short, "grouping": tag.replace("secB_", "SecB_"),
            "pvalue": pval, "n": len(df_sub), "followup": label_prefix or "Full",
        }, b64))
        print(f"{print_prefix}    {ep_short} {split_name}: p = {pval:.4f} (n={len(df_sub)})")
    return results


def run_cox_suite(df, time_col, event_col, secA_col, secB_col,
                  label_prefix, print_prefix=""):
    ep_short = "OS" if "os" in time_col else "PFS"
    results = []
    df_cox = df.dropna(subset=[time_col, event_col, secB_col]).copy()
    if len(df_cox) < 20:
        return results

    for var, var_name in [
        (secB_col, "SecB_fraction"), (secA_col, "SecA_fraction"),
        ("log2_secB_secA_ratio", "log2_SecB_SecA_ratio"),
        ("log2_secA_secB_ratio", "log2_SecA_SecB_ratio"),
        ("secA_epi_prop", "SecA_epi_proportion"),
        ("secB_epi_prop", "SecB_epi_proportion"),
    ]:
        if var not in df_cox.columns:
            continue
        df_uni = df_cox[[time_col, event_col, var]].dropna()
        if len(df_uni) < 20:
            continue
        cph = CoxPHFitter()
        try:
            cph.fit(df_uni, duration_col=time_col, event_col=event_col)
            s = cph.summary
            hr = np.exp(s["coef"].values[0])
            ci_lo = np.exp(s["coef lower 95%"].values[0])
            ci_hi = np.exp(s["coef upper 95%"].values[0])
            pval = s["p"].values[0]
            results.append({
                "endpoint": ep_short, "model": "univariate", "variable": var_name,
                "HR": round(hr, 4), "HR_CI_low": round(ci_lo, 4),
                "HR_CI_high": round(ci_hi, 4), "p_value": pval,
                "n": len(df_uni), "events": int(df_uni[event_col].sum()),
                "followup": label_prefix or "Full",
            })
            print(f"{print_prefix}      {var_name}: HR={hr:.3f} "
                  f"({ci_lo:.3f}-{ci_hi:.3f}), p={pval:.4f}")
        except Exception as e:
            print(f"{print_prefix}      {var_name}: Cox failed -- {e}")

    mv_vars = [secB_col]
    # Multivariate adjustment matches the manuscript model (epithelial fraction +
    # stage + age). Platinum sensitivity is intentionally NOT a covariate here
    # (author-approved; see manuscript-corrections). Platinum is still loaded and
    # used for subgroup analyses elsewhere.
    covar_map = {"stage_coded": "Stage", "age": "Age"}
    for var in covar_map:
        if var in df_cox.columns:
            mv_vars.append(var)
    df_mv = df_cox[[time_col, event_col] + mv_vars].dropna()
    if len(df_mv) >= 30:
        cph_mv = CoxPHFitter()
        try:
            cph_mv.fit(df_mv, duration_col=time_col, event_col=event_col)
            for var in mv_vars:
                s = cph_mv.summary.loc[var]
                results.append({
                    "endpoint": ep_short, "model": "multivariate",
                    "variable": covar_map.get(var, var),
                    "HR": round(np.exp(s["coef"]), 4),
                    "HR_CI_low": round(np.exp(s["coef lower 95%"]), 4),
                    "HR_CI_high": round(np.exp(s["coef upper 95%"]), 4),
                    "p_value": s["p"], "n": len(df_mv),
                    "events": int(df_mv[event_col].sum()), "followup": label_prefix or "Full",
                })
            print(f"{print_prefix}      Multivariate (n={len(df_mv)}): "
                  f"SecB HR={np.exp(cph_mv.summary.loc[secB_col, 'coef']):.3f}, "
                  f"p={cph_mv.summary.loc[secB_col, 'p']:.4f}")
        except Exception as e:
            print(f"{print_prefix}      Multivariate Cox failed: {e}")
    return results


def make_forest_plot(cox_results, endpoint, prefix=""):
    ep = [r for r in cox_results if r["endpoint"] == endpoint]
    if not ep:
        return None
    fig, ax = plt.subplots(1, 1, figsize=(7, max(3, len(ep) * 0.5)))
    for i, r in enumerate(reversed(ep)):
        hr, lo, hi = r["HR"], r["HR_CI_low"], r["HR_CI_high"]
        color = "#D55E00" if r["p_value"] < 0.05 else "#999999"
        ax.plot([lo, hi], [i, i], color=color, lw=2)
        ax.plot(hr, i, "o", color=color, markersize=8)
        ax.text(-0.05, i, f"{r['variable']} ({r['model']})",
                transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7)
        pstr = f"p={r['p_value']:.4f}" if r['p_value'] >= 0.0001 else "p<0.0001"
        ax.text(hi + 0.02, i, f"HR={hr:.2f} ({lo:.2f}-{hi:.2f}) {pstr}",
                va="center", fontsize=6)
    ax.axvline(1.0, color="black", ls="--", lw=0.8)
    ax.set_xlabel("Hazard Ratio")
    suffix = f" ({prefix})" if prefix else ""
    ax.set_title(f"Cox PH Forest Plot -- {endpoint}{suffix}", fontsize=9, weight="bold")
    ax.set_yticks([]); ax.set_xlim(left=0)
    fig.tight_layout()
    tag = f"_{prefix.lower().replace(' ', '_').replace('-', '')}" if prefix else ""
    return save_fig(fig, f"22c_cox_forest_{endpoint.lower()}{tag}")


# ============================================================================
# HTML REPORT
# ============================================================================

def build_html_report(sections, out_path):
    toc, body = "", ""
    for i, (title, content) in enumerate(sections):
        anchor = f"section-{i}"
        toc += f'<li><a href="#{anchor}">{title}</a></li>\n'
        body += f'<div class="section" id="{anchor}">\n<h2>{title}</h2>\n{content}\n</div>\n'
    html = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<title>TCGA Deconvolution Survival Report</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: Helvetica, Arial, sans-serif; font-size: 14px; line-height: 1.5;
          color: #333; max-width: 1100px; margin: 0 auto; padding: 20px 30px; background: #fafafa; }}
  h1 {{ font-size: 22px; border-bottom: 3px solid #0072B2; padding-bottom: 8px;
        margin-bottom: 6px; color: #0072B2; }}
  .subtitle {{ color: #666; font-size: 13px; margin-bottom: 20px; }}
  h2 {{ font-size: 17px; margin: 30px 0 12px 0; border-bottom: 1px solid #ddd; padding-bottom: 4px; }}
  h3 {{ font-size: 14px; margin: 18px 0 8px 0; color: #555; }}
  .section {{ background: #fff; border: 1px solid #e0e0e0; border-radius: 6px;
              padding: 18px 22px; margin-bottom: 18px; }}
  table {{ border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 13px; }}
  th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: left; }}
  th {{ background: #f5f5f5; font-weight: 600; }}
  img {{ max-width: 100%; height: auto; margin: 8px 0; border: 1px solid #eee; border-radius: 4px; }}
  .fig-row {{ display: flex; flex-wrap: wrap; gap: 12px; justify-content: center; }}
  .fig-row img {{ max-width: 48%; }}
  .toc {{ background: #fff; border: 1px solid #e0e0e0; border-radius: 6px;
          padding: 12px 20px; margin-bottom: 18px; }}
  .toc ul {{ list-style: none; padding-left: 0; columns: 2; }}
  .badge {{ display: inline-block; padding: 2px 8px; border-radius: 3px; font-size: 11px; font-weight: 600; }}
  .badge-sig {{ background: #D55E00; color: white; }}
  .badge-ns {{ background: #e0e0e0; color: #666; }}
  .badge-trend {{ background: #F0AD4E; color: white; }}
  .key-finding {{ background: #FFF8E1; border-left: 4px solid #F0AD4E;
                  padding: 10px 14px; margin: 10px 0; border-radius: 0 4px 4px 0; }}
</style></head><body>
<h1>TCGA Deconvolution Survival Analysis</h1>
<p class="subtitle">HGSC malignant-states atlas &middot;
Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
<div class="toc"><strong>Contents</strong><ul>{toc}</ul></div>
{body}
</body></html>"""
    with open(out_path, "w") as f:
        f.write(html)


def img_tag(b64):
    return f'<img src="data:image/png;base64,{b64}" />'


def pval_badge(p):
    if p < 0.05:
        return f'<span class="badge badge-sig">p={p:.4f}</span>'
    if p < 0.1:
        return f'<span class="badge badge-trend">p={p:.4f}</span>'
    return f'<span class="badge badge-ns">p={p:.4f}</span>'


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(description="Survival analysis on deconvolved fractions")
    parser.add_argument("--fractions", default=DEFAULT_FRACTIONS)
    args = parser.parse_args()

    print("=" * 65)
    print("  Survival analysis on deconvolved TCGA fractions")
    print("=" * 65)
    html_sections = []

    print("\n[1] Loading data...")
    if not os.path.exists(args.fractions):
        print(f"  ERROR: Results not found at {args.fractions}")
        sys.exit(1)
    fractions = pd.read_csv(args.fractions, sep="\t", index_col=0)
    print(f"    Fractions: {fractions.shape[0]} x {fractions.shape[1]}")
    qc_cols = ["P-value", "Correlation", "RMSE"]
    frac_cols = [c for c in fractions.columns if c not in qc_cols]
    fractions = fractions[frac_cols]
    print(f"    Cell types: {', '.join(frac_cols)}")
    clinical = pd.read_csv(CLINICAL, index_col=0)
    print(f"    Clinical: {clinical.shape[0]} x {clinical.shape[1]}")

    print("\n[2] Matching sample IDs...")

    def truncate_tcga_id(barcode):
        s = str(barcode).replace(".", "-")
        parts = s.split("-")
        return "-".join(parts[:3]) if len(parts) >= 3 else s

    frac_trunc = {truncate_tcga_id(s): s for s in fractions.index}
    clin_trunc = {truncate_tcga_id(s): s for s in clinical["PATIENT.NUMBER"]}
    matched_trunc = set(frac_trunc) & set(clin_trunc)
    print(f"    Matched: {len(matched_trunc)}")
    if not matched_trunc:
        print("    ERROR: No samples matched!")
        sys.exit(1)

    fractions_m = fractions.loc[[frac_trunc[t] for t in matched_trunc]].copy()
    fractions_m.index = [truncate_tcga_id(s) for s in fractions_m.index]
    clinical_m = clinical[clinical["PATIENT.NUMBER"].isin(
        [clin_trunc[t] for t in matched_trunc])].copy()
    clinical_m.index = clinical_m["PATIENT.NUMBER"].apply(truncate_tcga_id)
    clinical_m = clinical_m.drop(columns=["PATIENT.NUMBER"])
    for d in [fractions_m, clinical_m]:
        if d.index.duplicated().any():
            d = d[~d.index.duplicated(keep="first")]
    df = fractions_m.join(clinical_m, how="inner")
    print(f"    Joined: {len(df)} samples")

    print("\n[3] Parsing survival endpoints...")
    if "OS..Months." in df.columns and "OS" in df.columns:
        df["os_months"] = pd.to_numeric(df["OS..Months."], errors="coerce")
        df["os_event"] = df["OS"].apply(lambda x: 1 if "DECEASED" in str(x) else 0)
        n_os = df[["os_months", "os_event"]].dropna().shape[0]
        e_os = int(df.dropna(subset=["os_months", "os_event"])["os_event"].sum())
        print(f"    OS: {n_os} patients ({e_os} events)")
    else:
        df["os_months"] = np.nan; df["os_event"] = np.nan; n_os = e_os = 0
    if "PFS..Months." in df.columns and "PFS" in df.columns:
        df["pfs_months"] = pd.to_numeric(df["PFS..Months."], errors="coerce")
        df["pfs_event"] = df["PFS"].apply(lambda x: 1 if "PROGRESSION" in str(x) else 0)
        n_pfs = df[["pfs_months", "pfs_event"]].dropna().shape[0]
        e_pfs = int(df.dropna(subset=["pfs_months", "pfs_event"])["pfs_event"].sum())
        print(f"    PFS: {n_pfs} patients ({e_pfs} events)")
    else:
        df["pfs_months"] = np.nan; df["pfs_event"] = np.nan; n_pfs = e_pfs = 0

    print("\n[4] Computing derived variables...")
    epi_cols = [c for c in frac_cols if "epithelium" in c.lower()]
    secA_col = next((c for c in frac_cols if "seca" in c.lower()), None)
    secB_col = next((c for c in frac_cols if "secb" in c.lower()), None)
    trans_col = [c for c in frac_cols if "intermediate" in c.lower()]
    if not (secA_col and secB_col):
        print("    ERROR: SecA/SecB columns not found")
        sys.exit(1)
    print(f"    SecA: {secA_col},  SecB: {secB_col}")

    pseudo = 1e-6
    df["secB_secA_ratio"] = (df[secB_col] + pseudo) / (df[secA_col] + pseudo)
    df["log2_secB_secA_ratio"] = np.log2(df["secB_secA_ratio"])
    df["secA_secB_ratio"] = (df[secA_col] + pseudo) / (df[secB_col] + pseudo)
    df["log2_secA_secB_ratio"] = np.log2(df["secA_secB_ratio"])
    df["total_epi"] = df[epi_cols].sum(axis=1)
    df["secA_epi_prop"] = df[secA_col] / (df[epi_cols].sum(axis=1) + pseudo)
    df["secB_epi_prop"] = df[secB_col] / (df[epi_cols].sum(axis=1) + pseudo)

    secB_median = df[secB_col].median()
    df["secB_group_median"] = pd.cut(df[secB_col], bins=[-np.inf, secB_median, np.inf],
                                     labels=["SecB-low", "SecB-high"])
    t1, t2 = df[secB_col].quantile([1/3, 2/3])
    df["secB_group_tertile"] = pd.cut(df[secB_col], bins=[-np.inf, t1, t2, np.inf],
                                      labels=["SecB-low", "SecB-mid", "SecB-high"])
    ratio_median = df["log2_secB_secA_ratio"].median()
    df["ratio_group"] = pd.cut(df["log2_secB_secA_ratio"], bins=[-np.inf, ratio_median, np.inf],
                               labels=["SecA-dominant", "SecB-dominant"])
    ratio_ab_median = df["log2_secA_secB_ratio"].median()
    df["ratio_ab_group"] = pd.cut(df["log2_secA_secB_ratio"], bins=[-np.inf, ratio_ab_median, np.inf],
                                  labels=["Low SecA:SecB", "High SecA:SecB"])
    secA_prop_median = df["secA_epi_prop"].median()
    df["secA_prop_group"] = pd.cut(df["secA_epi_prop"], bins=[-np.inf, secA_prop_median, np.inf],
                                   labels=["Low SecA prop", "High SecA prop"])

    if "Stage..Coded." in df.columns:
        df["stage_coded"] = pd.to_numeric(df["Stage..Coded."], errors="coerce")
    if "Age" in df.columns:
        df["age"] = pd.to_numeric(df["Age"], errors="coerce")
    if "Platinum..Coded." in df.columns:
        df["platinum_coded"] = pd.to_numeric(df["Platinum..Coded."], errors="coerce")

    df.to_csv(os.path.join(OUT_DIR, "22c_deconv_clinical_joined.csv"))
    print(f"    Saved: 22c_deconv_clinical_joined.csv ({len(df)} samples)")

    if not HAS_LIFELINES:
        print("\n  lifelines not installed -- skipping survival")
        return

    frac_summary_rows = [{
        "Cell Type": ct, "Mean": f"{df[ct].mean():.4f}", "Median": f"{df[ct].median():.4f}",
        "SD": f"{df[ct].std():.4f}", "Min": f"{df[ct].min():.4f}", "Max": f"{df[ct].max():.4f}",
    } for ct in frac_cols]
    frac_summary_html = pd.DataFrame(frac_summary_rows).to_html(index=False, classes="", border=0)
    overview_html = f"""
    <p>BayesPrism deconvolution of <b>{len(df)} TCGA-OV samples</b> using an
    atlas-derived scRNA-seq reference with <b>{len(frac_cols)} cell types</b>.</p>
    <h3>Deconvolved Cell Type Fractions</h3>{frac_summary_html}
    <h3>Key Parameters</h3><table>
    <tr><td>SecB median cutoff</td><td>{secB_median:.4f}</td></tr>
    <tr><td>SecB tertile cutoffs</td><td>{t1:.4f}, {t2:.4f}</td></tr>
    <tr><td>OS patients (events)</td><td>{n_os} ({e_os})</td></tr>
    <tr><td>PFS patients (events)</td><td>{n_pfs} ({e_pfs})</td></tr>
    </table>"""
    html_sections.append(("Overview & Fractions", overview_html))

    print("\n[5] Generating fraction heatmap & correlations...")
    frac_data = df[frac_cols].sort_values(secB_col, ascending=False)
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    sns.heatmap(frac_data.T, cmap="YlOrRd", xticklabels=False, yticklabels=True, ax=ax,
                cbar_kws={"shrink": 0.6, "label": "Estimated fraction"})
    ax.set_title("BayesPrism Estimated Cell Type Fractions (TCGA-OV)", fontsize=9, weight="bold")
    ax.set_ylabel("")
    fig.tight_layout()
    heatmap_b64 = save_fig(fig, "22c_fraction_heatmap")

    tme_cols = [c for c in frac_cols if c not in [secA_col, secB_col] + (trans_col or [])]
    corr_b64 = None
    if tme_cols:
        corr_rows = []
        for tme in tme_cols:
            for epi, en in [(secA_col, "SecA"), (secB_col, "SecB")]:
                corr_rows.append({"epithelial": en, "TME": tme,
                                  "pearson_r": round(df[[epi, tme]].dropna().corr().iloc[0, 1], 4)})
        corr_pivot = pd.DataFrame(corr_rows).pivot(index="TME", columns="epithelial", values="pearson_r")
        fig, ax = plt.subplots(1, 1, figsize=(5, max(3, len(tme_cols) * 0.4)))
        sns.heatmap(corr_pivot, annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                    vmin=-0.5, vmax=0.5, ax=ax,
                    cbar_kws={"shrink": 0.7, "label": "Pearson r"})
        ax.set_title("Epithelial vs TME Fraction Correlations", fontsize=9, weight="bold")
        fig.tight_layout()
        corr_b64 = save_fig(fig, "22c_fraction_correlations")
    viz_html = f'<h3>Fraction Heatmap (sorted by SecB)</h3>\n{img_tag(heatmap_b64)}\n'
    if corr_b64:
        viz_html += f'<h3>SecA/SecB vs TME Correlations</h3>\n{img_tag(corr_b64)}\n'
    html_sections.append(("Fraction Visualizations", viz_html))

    all_km, all_cox = [], []
    for max_mo, label, prefix in [(None, "Full follow-up", ""), (60, "5-year (60 months)", "_5yr")]:
        print(f"\n{'='*50}\n  Survival: {label}\n{'='*50}")
        if max_mo is not None:
            df_s = clip_survival(clip_survival(df, "os_months", "os_event", max_mo),
                                 "pfs_months", "pfs_event", max_mo)
            n_os_c = df_s.dropna(subset=["os_months", "os_event"]).shape[0]
            e_os_c = int(df_s.dropna(subset=["os_months", "os_event"])["os_event"].sum())
            n_pfs_c = df_s.dropna(subset=["pfs_months", "pfs_event"]).shape[0]
            e_pfs_c = int(df_s.dropna(subset=["pfs_months", "pfs_event"])["pfs_event"].sum())
        else:
            df_s = df

        km_section_imgs = []
        for time_col, event_col in [("os_months", "os_event"), ("pfs_months", "pfs_event")]:
            for rd, b64 in run_km_suite(df_s, time_col, event_col, None, secB_col,
                                        prefix, label, "  "):
                all_km.append(rd)
                km_section_imgs.append((rd, b64))

        cox_this = []
        for time_col, event_col in [("os_months", "os_event"), ("pfs_months", "pfs_event")]:
            ep_short = "OS" if "os" in time_col else "PFS"
            print(f"\n    --- Cox {ep_short} ({label}) ---")
            cr = run_cox_suite(df_s, time_col, event_col, secA_col, secB_col, label, "  ")
            cox_this.extend(cr)
            all_cox.extend(cr)

        forest_imgs = {}
        for ep in ["OS", "PFS"]:
            b64 = make_forest_plot(cox_this, ep, label)
            if b64:
                forest_imgs[ep] = b64

        sec_html = ""
        if max_mo is not None:
            sec_html += (f"<p>Survival censored at <b>{max_mo} months</b>. "
                         f"OS: {n_os_c} patients ({e_os_c} events), "
                         f"PFS: {n_pfs_c} patients ({e_pfs_c} events).</p>\n")
        sec_html += "<h3>Kaplan-Meier Log-Rank Tests</h3>\n<table>"
        sec_html += "<tr><th>Endpoint</th><th>Grouping</th><th>n</th><th>p-value</th></tr>\n"
        for rd, _ in km_section_imgs:
            sec_html += (f"<tr><td>{rd['endpoint']}</td><td>{rd['grouping']}</td>"
                         f"<td>{rd['n']}</td><td>{pval_badge(rd['pvalue'])}</td></tr>\n")
        sec_html += "</table>\n<h3>Kaplan-Meier Curves</h3>\n<div class=\"fig-row\">\n"
        for _, b64 in km_section_imgs:
            sec_html += img_tag(b64) + "\n"
        sec_html += "</div>\n"
        if cox_this:
            sec_html += "<h3>Cox Proportional Hazards</h3>\n<table>"
            sec_html += ("<tr><th>Endpoint</th><th>Model</th><th>Variable</th><th>HR</th>"
                         "<th>95% CI</th><th>p-value</th><th>n</th><th>Events</th></tr>\n")
            for r in cox_this:
                sec_html += (f"<tr><td>{r['endpoint']}</td><td>{r['model']}</td>"
                             f"<td>{r['variable']}</td><td>{r['HR']:.3f}</td>"
                             f"<td>{r['HR_CI_low']:.3f}-{r['HR_CI_high']:.3f}</td>"
                             f"<td>{pval_badge(r['p_value'])}</td><td>{r['n']}</td>"
                             f"<td>{r['events']}</td></tr>\n")
            sec_html += "</table>\n"
        if forest_imgs:
            sec_html += '<h3>Forest Plots</h3>\n<div class="fig-row">\n'
            for ep, b64 in forest_imgs.items():
                sec_html += img_tag(b64) + "\n"
            sec_html += "</div>\n"
        html_sections.append((f"Survival -- {label}", sec_html))

    if all_km:
        pd.DataFrame(all_km).to_csv(os.path.join(OUT_DIR, "22c_km_logrank_results.csv"), index=False)
    if all_cox:
        pd.DataFrame(all_cox).to_csv(os.path.join(OUT_DIR, "22c_cox_results.csv"), index=False)

    build_html_report(html_sections, os.path.join(OUT_DIR, "22c_survival_report.html"))
    print("\n    Saved: 22c_survival_report.html")
    print(f"\n{'='*65}\n  Step complete! Output: {OUT_DIR}\n{'='*65}")


if __name__ == "__main__":
    main()
