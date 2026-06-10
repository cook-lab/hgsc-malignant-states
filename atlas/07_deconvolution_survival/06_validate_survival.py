#!/usr/bin/env python3
"""
Rigorous validation of the SecA/SecB polarization survival signal
=================================================================
HGSC malignant-states atlas backend.

Validates the 05_signature_survival.py finding (polarization HR~0.5, p<0.01) with:
  1. Permutation test (10,000 label shuffles -> empirical p)
  2. Bootstrap CI on the Cox HR (1,000 resamples)
  3. Alternative scoring (ssGSEA-like rank-based; top-50 / top-25 restricted)
  4. Leave-one-out influence (does one patient drive the result?)
  5. Subgroup analyses (platinum status, stage, age)
  6. Comparison vs the merged-BayesPrism fractions (if 03 has run)

INPUTS (output_root/07_deconvolution_survival/):
  - 22d_signature_scores.csv (from 05_signature_survival.py)
  - bayesprism_merged_results.txt (optional, from 03_bayesprism_merged.R)
  - CIBERSORTx_Results.txt
  - <data_root>/2026_final_atlas/output/11e_nmf_characterization/11e_gene_classification.csv
  - <data_root>/2026_final_atlas/data/cibersort_data_prev/{tcga_ecotyper.txt,tcga_hla_clinical.csv}

OUTPUTS (output_root/07_deconvolution_survival/):
  - 22f_validation_report.html
  - 22f_permutation_null.csv, 22f_bootstrap_hr.csv
  - 22f_{permutation_null,bootstrap_hr,loo_influence,merged_comparison}.svg/pdf

MANUSCRIPT PANELS: robustness support for Fig 7E/F/G.

RUNTIME TIER: heavy (10k permutations + 1k bootstraps + LOO Cox refits).

SEEDING: permutation / bootstrap / subgroup RNGs all use config SEED.

Usage:
    python 06_validate_survival.py
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
from config.config import path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

from lifelines import CoxPHFitter

# ── Paths ─────────────────────────────────────────────────────
OUT_DIR     = path("output_root", "07_deconvolution_survival")
SCORES      = os.path.join(OUT_DIR, "22d_signature_scores.csv")
# Prefer the regenerated classification (atlas/03/01b_gene_classification.py); fall
# back to the deposited copy if a clean re-run has not produced it (audit A4/H4).
_GC_RECOMPUTED = path("output_root", "03_epithelial_nmf", "11e_gene_classification.csv")
_GC_DEPOSITED  = path("data_root", "2026_final_atlas", "output", "11e_nmf_characterization",
                      "11e_gene_classification.csv")
GENE_CLS    = _GC_RECOMPUTED if os.path.exists(_GC_RECOMPUTED) else _GC_DEPOSITED
TCGA_EXPR   = path("data_root", "2026_final_atlas", "data", "cibersort_data_prev",
                   "tcga_ecotyper.txt")
CLINICAL    = path("data_root", "2026_final_atlas", "data", "cibersort_data_prev",
                   "tcga_hla_clinical.csv")
MERGED_FRAC = os.path.join(OUT_DIR, "bayesprism_merged_results.txt")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial"],
    "font.size": 8, "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 450,
})
BLUE, ORANGE, GREY = "#0072B2", "#D55E00", "#999999"


def fig_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    buf.close()
    return b64


def save_fig(fig, name):
    fig.savefig(os.path.join(OUT_DIR, f"{name}.svg"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{name}.pdf"), bbox_inches="tight")
    b64 = fig_b64(fig)
    plt.close(fig)
    return b64


def img(b64):
    return f'<img src="data:image/png;base64,{b64}" />'


def pval_badge(p):
    if p < 0.05:
        return f'<span style="color:#D55E00;font-weight:700">p={p:.4f}</span>'
    if p < 0.1:
        return f'<span style="color:#F0AD4E;font-weight:600">p={p:.4f}</span>'
    return f'<span style="color:#999">p={p:.4f}</span>'


def truncate_id(bc):
    return "-".join(str(bc).replace(".", "-").split("-")[:3])


def main():
    print("=" * 65)
    print("  Validation of survival signal")
    print("=" * 65)
    html = []

    df = pd.read_csv(SCORES, index_col=0)
    print(f"\n[0] Loaded: {len(df)} samples")
    df_os = df.dropna(subset=["os_months", "os_event"]).copy()
    df_pfs = df.dropna(subset=["pfs_months", "pfs_event"]).copy()
    print(f"    OS valid: {len(df_os)} ({int(df_os['os_event'].sum())} events)")
    print(f"    PFS valid: {len(df_pfs)} ({int(df_pfs['pfs_event'].sum())} events)")

    # ── 1. Reproduce original ─────────────────────────────────
    print("\n[1] Reproducing original result...")
    orig_results = {}
    for label, dfs, tcol, ecol in [("OS", df_os, "os_months", "os_event"),
                                   ("PFS", df_pfs, "pfs_months", "pfs_event")]:
        cph = CoxPHFitter()
        cph.fit(dfs[[tcol, ecol, "polarization"]], tcol, ecol)
        s = cph.summary.loc["polarization"]
        orig_results[label] = {"hr": np.exp(s["coef"]), "p": s["p"], "coef": s["coef"],
                               "n": len(dfs), "events": int(dfs[ecol].sum())}
        print(f"    {label}: HR={np.exp(s['coef']):.4f}, p={s['p']:.6f}")

    # ── 2. Permutation test ───────────────────────────────────
    N_PERM = 10000
    print(f"\n[2] Permutation test ({N_PERM} shuffles)...")
    perm_results = {}
    for label, dfs, tcol, ecol in [("OS", df_os, "os_months", "os_event"),
                                   ("PFS", df_pfs, "pfs_months", "pfs_event")]:
        observed_coef = orig_results[label]["coef"]
        null_coefs = []
        rng = np.random.RandomState(SEED)
        for i in range(N_PERM):
            perm_df = dfs[[tcol, ecol]].copy()
            perm_df["polarization"] = rng.permutation(dfs["polarization"].values)
            cph = CoxPHFitter(penalizer=0.01)
            try:
                cph.fit(perm_df, tcol, ecol)
                null_coefs.append(cph.summary.loc["polarization", "coef"])
            except Exception:
                null_coefs.append(0.0)
            if (i + 1) % 2000 == 0:
                print(f"    {label}: {i+1}/{N_PERM}")
        null_coefs = np.array(null_coefs)
        emp_p = (np.sum(np.abs(null_coefs) >= np.abs(observed_coef)) + 1) / (N_PERM + 1)
        perm_results[label] = {"emp_p": emp_p, "null_coefs": null_coefs, "observed": observed_coef}
        print(f"    {label}: empirical p = {emp_p:.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, label in zip(axes, ["OS", "PFS"]):
        null = perm_results[label]["null_coefs"]
        obs = perm_results[label]["observed"]
        ax.hist(null, bins=50, color=GREY, alpha=0.7, edgecolor="white", density=True)
        ax.axvline(obs, color=ORANGE, lw=2, ls="--", label=f"Observed (coef={obs:.3f})")
        ax.set_title(f"{label} Permutation Test (n={N_PERM:,})", fontsize=9, weight="bold")
        ax.set_xlabel("Cox coefficient (polarization)"); ax.set_ylabel("Density")
        ax.text(0.02, 0.95, f"Empirical p = {perm_results[label]['emp_p']:.4f}",
                transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(fc="white", ec="gray", alpha=0.8))
        ax.legend(fontsize=7)
    fig.tight_layout()
    perm_b64 = save_fig(fig, "22f_permutation_null")
    pd.DataFrame({"OS_null": perm_results["OS"]["null_coefs"],
                  "PFS_null": perm_results["PFS"]["null_coefs"]}).to_csv(
        os.path.join(OUT_DIR, "22f_permutation_null.csv"), index=False)
    perm_html = f"""
    <p>Permutation test: shuffle polarization {N_PERM:,} times, refit Cox.</p>
    <table>
    <tr><th>Endpoint</th><th>Observed coef</th><th>HR</th><th>Parametric p</th><th>Empirical p</th></tr>
    <tr><td>OS</td><td>{orig_results['OS']['coef']:.4f}</td><td>{orig_results['OS']['hr']:.4f}</td>
        <td>{pval_badge(orig_results['OS']['p'])}</td><td>{pval_badge(perm_results['OS']['emp_p'])}</td></tr>
    <tr><td>PFS</td><td>{orig_results['PFS']['coef']:.4f}</td><td>{orig_results['PFS']['hr']:.4f}</td>
        <td>{pval_badge(orig_results['PFS']['p'])}</td><td>{pval_badge(perm_results['PFS']['emp_p'])}</td></tr>
    </table>{img(perm_b64)}"""
    html.append(("Permutation Test", perm_html))

    # ── 3. Bootstrap CI ───────────────────────────────────────
    N_BOOT = 1000
    print(f"\n[3] Bootstrap ({N_BOOT} resamples)...")
    boot_results = {}
    for label, dfs, tcol, ecol in [("OS", df_os, "os_months", "os_event"),
                                   ("PFS", df_pfs, "pfs_months", "pfs_event")]:
        boot_hrs = []
        rng = np.random.RandomState(SEED)
        for i in range(N_BOOT):
            idx = rng.choice(len(dfs), size=len(dfs), replace=True)
            boot_df = dfs.iloc[idx][[tcol, ecol, "polarization"]].reset_index(drop=True)
            cph = CoxPHFitter(penalizer=0.01)
            try:
                cph.fit(boot_df, tcol, ecol)
                boot_hrs.append(np.exp(cph.summary.loc["polarization", "coef"]))
            except Exception:
                pass
        boot_hrs = np.array(boot_hrs)
        ci_lo, ci_hi = np.percentile(boot_hrs, [2.5, 97.5])
        boot_results[label] = {"hrs": boot_hrs, "ci_lo": ci_lo, "ci_hi": ci_hi,
                               "median": np.median(boot_hrs)}
        print(f"    {label}: median HR={np.median(boot_hrs):.4f}, 95% CI=[{ci_lo:.4f}, {ci_hi:.4f}]")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, label in zip(axes, ["OS", "PFS"]):
        hrs = boot_results[label]["hrs"]
        ci_lo, ci_hi = boot_results[label]["ci_lo"], boot_results[label]["ci_hi"]
        ax.hist(hrs, bins=40, color=BLUE, alpha=0.7, edgecolor="white")
        ax.axvline(1.0, color="black", ls="--", lw=1)
        ax.axvline(ci_lo, color=ORANGE, ls=":", lw=1.5, label=f"2.5%: {ci_lo:.3f}")
        ax.axvline(ci_hi, color=ORANGE, ls=":", lw=1.5, label=f"97.5%: {ci_hi:.3f}")
        ax.axvline(orig_results[label]["hr"], color=ORANGE, ls="-", lw=2,
                   label=f"Point est: {orig_results[label]['hr']:.3f}")
        ax.set_title(f"{label} Bootstrap HR (n={N_BOOT})", fontsize=9, weight="bold")
        ax.set_xlabel("Hazard Ratio"); ax.legend(fontsize=7)
    fig.tight_layout()
    boot_b64 = save_fig(fig, "22f_bootstrap_hr")
    pd.DataFrame({"OS_HR": boot_results["OS"]["hrs"],
                  "PFS_HR": boot_results["PFS"]["hrs"]}).to_csv(
        os.path.join(OUT_DIR, "22f_bootstrap_hr.csv"), index=False)
    boot_html = f"""
    <p>Non-parametric bootstrap: resample patients with replacement {N_BOOT} times.</p>
    <table>
    <tr><th>Endpoint</th><th>Point HR</th><th>Bootstrap median HR</th><th>Bootstrap 95% CI</th><th>CI excludes 1?</th></tr>
    <tr><td>OS</td><td>{orig_results['OS']['hr']:.4f}</td><td>{boot_results['OS']['median']:.4f}</td>
        <td>[{boot_results['OS']['ci_lo']:.4f}, {boot_results['OS']['ci_hi']:.4f}]</td>
        <td>{'<b style="color:#D55E00">YES</b>' if boot_results['OS']['ci_hi'] < 1 else 'No'}</td></tr>
    <tr><td>PFS</td><td>{orig_results['PFS']['hr']:.4f}</td><td>{boot_results['PFS']['median']:.4f}</td>
        <td>[{boot_results['PFS']['ci_lo']:.4f}, {boot_results['PFS']['ci_hi']:.4f}]</td>
        <td>{'<b style="color:#D55E00">YES</b>' if boot_results['PFS']['ci_hi'] < 1 else 'No'}</td></tr>
    </table>{img(boot_b64)}"""
    html.append(("Bootstrap Confidence Intervals", boot_html))

    # ── 4. Alternative scoring methods ────────────────────────
    print("\n[4] Alternative scoring methods...")
    gc = pd.read_csv(GENE_CLS)
    secA_genes = gc[gc["class"] == "SecA-specific"]["gene"].tolist()
    secB_genes = gc[gc["class"] == "SecB-specific"]["gene"].tolist()
    expr = pd.read_csv(TCGA_EXPR, sep="\t", index_col=0)
    secA_in = [g for g in secA_genes if g in expr.index]
    secB_in = [g for g in secB_genes if g in expr.index]

    expr_rank = expr.rank(axis=0, ascending=True, pct=True)
    ssgsea_polar = expr_rank.loc[secA_in].mean(axis=0) - expr_rank.loc[secB_in].mean(axis=0)
    expr_z = expr.apply(lambda row: (row - row.mean()) / (row.std() + 1e-10), axis=1)

    def topn_polar(n):
        a = [g for g in gc[gc["class"] == "SecA-specific"].nlargest(n, "secA_loading")["gene"]
             if g in expr.index]
        b = [g for g in gc[gc["class"] == "SecB-specific"].nlargest(n, "secB_loading")["gene"]
             if g in expr.index]
        return expr_z.loc[a].mean(axis=0) - expr_z.loc[b].mean(axis=0)

    alt_methods = {
        "Mean z-score (177 genes)": df["polarization"],
        "ssGSEA-like (rank-based)": ssgsea_polar,
        "Top 50 genes (z-score)": topn_polar(50),
        "Top 25 genes (z-score)": topn_polar(25),
    }
    alt_results = []
    for method_name, sc in alt_methods.items():
        for label, dfs, tcol, ecol in [("OS", df_os, "os_months", "os_event"),
                                       ("PFS", df_pfs, "pfs_months", "pfs_event")]:
            common_idx = dfs.index.intersection(sc.index)
            if len(common_idx) < 20:
                score_map = {truncate_id(s): sc[s] for s in sc.index}
                test_df = dfs[[tcol, ecol]].copy()
                test_df["score"] = [score_map.get(truncate_id(s), np.nan) for s in test_df.index]
                test_df = test_df.dropna()
            else:
                test_df = dfs[[tcol, ecol]].copy()
                test_df["score"] = sc.reindex(test_df.index)
                test_df = test_df.dropna()
            if len(test_df) < 20:
                continue
            cph = CoxPHFitter()
            try:
                cph.fit(test_df, tcol, ecol)
                s = cph.summary.loc["score"]
                alt_results.append({"method": method_name, "endpoint": label,
                                    "HR": np.exp(s["coef"]), "p_value": s["p"], "n": len(test_df)})
                star = "***" if s["p"] < 0.05 else ("  ~" if s["p"] < 0.1 else "   ")
                print(f"  {star} {method_name} {label}: HR={np.exp(s['coef']):.3f}, p={s['p']:.4f}")
            except Exception as e:
                print(f"      {method_name} {label}: failed -- {e}")

    alt_html = "<table><tr><th>Scoring Method</th><th>Endpoint</th><th>HR</th><th>p-value</th><th>n</th></tr>\n"
    for r in alt_results:
        alt_html += (f"<tr><td>{r['method']}</td><td>{r['endpoint']}</td><td>{r['HR']:.4f}</td>"
                     f"<td>{pval_badge(r['p_value'])}</td><td>{r['n']}</td></tr>\n")
    alt_html += "</table>\n<p><i>All methods test polarization (SecA - SecB) as continuous Cox predictor.</i></p>"
    html.append(("Alternative Scoring Methods", alt_html))

    # ── 5. Leave-one-out influence ────────────────────────────
    print("\n[5] Leave-one-out influence check...")
    loo_results = {}
    for label, dfs, tcol, ecol in [("OS", df_os, "os_months", "os_event"),
                                   ("PFS", df_pfs, "pfs_months", "pfs_event")]:
        loo_hrs, loo_ps = [], []
        for i in range(len(dfs)):
            loo_df = dfs.drop(dfs.index[i])[[tcol, ecol, "polarization"]]
            cph = CoxPHFitter()
            try:
                cph.fit(loo_df, tcol, ecol)
                loo_hrs.append(np.exp(cph.summary.loc["polarization", "coef"]))
                loo_ps.append(cph.summary.loc["polarization", "p"])
            except Exception:
                loo_hrs.append(np.nan); loo_ps.append(np.nan)
        loo_hrs = np.array(loo_hrs); loo_ps = np.array(loo_ps)
        loo_results[label] = {"hrs": loo_hrs, "ps": loo_ps,
                              "n_flip": int(np.sum(loo_ps >= 0.05)), "n_total": len(dfs)}
        print(f"    {label}: HR range [{np.nanmin(loo_hrs):.4f}, {np.nanmax(loo_hrs):.4f}], "
              f"{loo_results[label]['n_flip']}/{len(dfs)} drops lose significance")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, label in zip(axes, ["OS", "PFS"]):
        hrs = loo_results[label]["hrs"]
        ax.plot(sorted(hrs), np.linspace(0, 1, len(hrs)), color=BLUE, lw=1.5)
        ax.axvline(orig_results[label]["hr"], color=ORANGE, ls="--", lw=1.5,
                   label=f"Full model: {orig_results[label]['hr']:.3f}")
        ax.axvline(1.0, color="black", ls=":", lw=0.8)
        ax.set_xlabel("HR (leave-one-out)"); ax.set_ylabel("Cumulative proportion")
        ax.set_title(f"{label} LOO Influence "
                     f"({loo_results[label]['n_flip']}/{loo_results[label]['n_total']} lose p<0.05)",
                     fontsize=9, weight="bold")
        ax.legend(fontsize=7)
    fig.tight_layout()
    loo_b64 = save_fig(fig, "22f_loo_influence")
    loo_html = f"""
    <p>Remove each patient one at a time, refit Cox.</p>
    <table>
    <tr><th>Endpoint</th><th>Full HR</th><th>LOO HR range</th><th>Drops losing p&lt;0.05</th><th>Stable?</th></tr>
    <tr><td>OS</td><td>{orig_results['OS']['hr']:.4f}</td>
        <td>[{np.nanmin(loo_results['OS']['hrs']):.4f}, {np.nanmax(loo_results['OS']['hrs']):.4f}]</td>
        <td>{loo_results['OS']['n_flip']}/{loo_results['OS']['n_total']}</td>
        <td>{'<b style="color:#D55E00">YES</b>' if loo_results['OS']['n_flip']/loo_results['OS']['n_total'] < 0.1 else 'Fragile'}</td></tr>
    <tr><td>PFS</td><td>{orig_results['PFS']['hr']:.4f}</td>
        <td>[{np.nanmin(loo_results['PFS']['hrs']):.4f}, {np.nanmax(loo_results['PFS']['hrs']):.4f}]</td>
        <td>{loo_results['PFS']['n_flip']}/{loo_results['PFS']['n_total']}</td>
        <td>{'<b style="color:#D55E00">YES</b>' if loo_results['PFS']['n_flip']/loo_results['PFS']['n_total'] < 0.1 else 'Fragile'}</td></tr>
    </table>{img(loo_b64)}"""
    html.append(("Leave-One-Out Influence", loo_html))

    # ── 6. Merged BayesPrism comparison (if available) ────────
    if os.path.exists(MERGED_FRAC):
        print("\n[6] Comparing with merged BayesPrism fractions...")
        mfrac = pd.read_csv(MERGED_FRAC, sep="\t", index_col=0)
        mfrac_t = {truncate_id(s): s for s in mfrac.index}
        orig_frac = pd.read_csv(os.path.join(OUT_DIR, "CIBERSORTx_Results.txt"),
                                sep="\t", index_col=0)
        epi_cols_orig = [c for c in orig_frac.columns if "epithelium" in c.lower()]
        orig_frac["epi_total"] = orig_frac[epi_cols_orig].sum(axis=1)
        orig_t = {truncate_id(s): s for s in orig_frac.index}
        common_t = set(orig_t) & set(mfrac_t)
        if common_t:
            epi_orig = [orig_frac.loc[orig_t[t], "epi_total"] for t in common_t]
            sec_col = [c for c in mfrac.columns if "secretory" in c.lower()]
            cil_col = [c for c in mfrac.columns if "ciliated" in c.lower()]
            if sec_col and cil_col:
                epi_merged = [mfrac.loc[mfrac_t[t], sec_col[0]] + mfrac.loc[mfrac_t[t], cil_col[0]]
                              for t in common_t]
            else:
                epi_cols_m = [c for c in mfrac.columns if "epithelium" in c.lower()]
                epi_merged = [mfrac.loc[mfrac_t[t], epi_cols_m].sum() for t in common_t]
            r_val = np.corrcoef(epi_orig, epi_merged)[0, 1]
            print(f"    Epithelial fraction correlation (orig vs merged): r={r_val:.4f}")
            fig, ax = plt.subplots(1, 1, figsize=(5, 5))
            ax.scatter(epi_orig, epi_merged, s=10, alpha=0.5, color=BLUE)
            ax.plot([0, 1], [0, 1], "k--", lw=0.8)
            ax.set_xlabel("Original (4 epi types summed)")
            ax.set_ylabel("Merged (Secretory + Ciliated)")
            ax.set_title(f"Epithelial Fraction: Original vs Merged (r={r_val:.3f})",
                         fontsize=9, weight="bold")
            fig.tight_layout()
            merged_b64 = save_fig(fig, "22f_merged_comparison")
            merged_survival = ""
            for label, dfs, tcol, ecol in [("OS", df_os, "os_months", "os_event"),
                                           ("PFS", df_pfs, "pfs_months", "pfs_event")]:
                test_df = dfs[[tcol, ecol, "polarization"]].copy()
                test_df["epi_merged"] = ([
                    mfrac.loc[mfrac_t[truncate_id(s)], sec_col[0]] +
                    mfrac.loc[mfrac_t[truncate_id(s)], cil_col[0]]
                    if truncate_id(s) in mfrac_t else np.nan
                    for s in test_df.index] if sec_col and cil_col else np.nan)
                test_df = test_df.dropna()
                if len(test_df) >= 20:
                    cph = CoxPHFitter()
                    cph.fit(test_df, tcol, ecol)
                    s = cph.summary.loc["polarization"]
                    merged_survival += (f"<tr><td>{label}</td><td>{np.exp(s['coef']):.4f}</td>"
                                        f"<td>{pval_badge(s['p'])}</td><td>{len(test_df)}</td></tr>\n")
            merged_html = f"""
            <p>Comparison of epithelial fraction estimates: original (4 subtypes summed) vs
            merged (Secretory + Ciliated) BayesPrism runs.</p>{img(merged_b64)}
            <h3>Survival with Merged Epi Fraction as Covariate</h3>
            <table><tr><th>Endpoint</th><th>Polarization HR</th><th>p-value</th><th>n</th></tr>
            {merged_survival}</table>"""
            html.append(("Merged BayesPrism Comparison", merged_html))
        else:
            print("    No common samples found")
    else:
        print("\n[6] Merged BayesPrism not available (run 03 first)")
        html.append(("Merged BayesPrism Comparison",
                     "<p><i>Merged BayesPrism run not yet complete. Re-run after 03 finishes.</i></p>"))

    # ── 7. Subgroup analyses ──────────────────────────────────
    print("\n[7] Subgroup analyses...")
    subgroup_results = []

    def subgroup_cox(sub, name):
        for label, tcol, ecol in [("OS", "os_months", "os_event"),
                                  ("PFS", "pfs_months", "pfs_event")]:
            sub_s = sub.dropna(subset=[tcol, ecol, "polarization"])
            if len(sub_s) < 15:
                continue
            cph = CoxPHFitter(penalizer=0.01)
            try:
                cph.fit(sub_s[[tcol, ecol, "polarization"]], tcol, ecol)
                s = cph.summary.loc["polarization"]
                subgroup_results.append({"subgroup": name, "endpoint": label,
                                         "HR": np.exp(s["coef"]), "p_value": s["p"], "n": len(sub_s)})
                star = "***" if s["p"] < 0.05 else "   "
                print(f"    {star} {name} {label}: HR={np.exp(s['coef']):.3f}, p={s['p']:.4f} (n={len(sub_s)})")
            except Exception:
                pass

    for plat_val, plat_name in [(1, "Platinum-sensitive"), (2, "Platinum-resistant")]:
        subgroup_cox(df[df["platinum_coded"] == plat_val], plat_name)
    for stage, sname in [(3, "Stage III"), (4, "Stage IV")]:
        subgroup_cox(df[df["stage_coded"] == stage], sname)
    age_med = df["age"].median()
    subgroup_cox(df[df["age"] < age_med], f"Age < {age_med:.0f}")
    subgroup_cox(df[df["age"] >= age_med], f"Age >= {age_med:.0f}")

    sub_html = "<table><tr><th>Subgroup</th><th>Endpoint</th><th>HR</th><th>p-value</th><th>n</th></tr>\n"
    for r in subgroup_results:
        sub_html += (f"<tr><td>{r['subgroup']}</td><td>{r['endpoint']}</td><td>{r['HR']:.4f}</td>"
                     f"<td>{pval_badge(r['p_value'])}</td><td>{r['n']}</td></tr>\n")
    sub_html += "</table>\n<p><i>All subgroup tests use polarization (SecA-SecB) as univariate Cox predictor.</i></p>"
    html.append(("Subgroup Analyses", sub_html))

    # ── 8. Summary ────────────────────────────────────────────
    n_alt_sig = sum(1 for r in alt_results if r["p_value"] < 0.05)
    n_sub_sig = sum(1 for r in subgroup_results if r["p_value"] < 0.05)
    summary_html = ('<div style="background:#FFF8E1;border-left:4px solid #F0AD4E;padding:12px;'
                    'margin:10px 0;border-radius:4px">\n<h3>Validation Summary</h3>\n<ul>\n')
    summary_html += (f"<li><b>Permutation test</b>: OS emp-p={perm_results['OS']['emp_p']:.4f}, "
                     f"PFS emp-p={perm_results['PFS']['emp_p']:.4f}</li>\n")
    summary_html += (f"<li><b>Bootstrap 95% CI</b>: OS [{boot_results['OS']['ci_lo']:.3f}, "
                     f"{boot_results['OS']['ci_hi']:.3f}], PFS [{boot_results['PFS']['ci_lo']:.3f}, "
                     f"{boot_results['PFS']['ci_hi']:.3f}]</li>\n")
    summary_html += f"<li><b>Alternative scoring</b>: {n_alt_sig}/{len(alt_results)} significant</li>\n"
    summary_html += (f"<li><b>LOO stability</b>: OS {loo_results['OS']['n_flip']}/{loo_results['OS']['n_total']} "
                     f"lose sig, PFS {loo_results['PFS']['n_flip']}/{loo_results['PFS']['n_total']} lose sig</li>\n")
    summary_html += f"<li><b>Subgroups significant</b>: {n_sub_sig}/{len(subgroup_results)}</li>\n"
    summary_html += "</ul></div>\n"
    html.append(("Validation Summary", summary_html))

    # ── Build HTML ────────────────────────────────────────────
    toc = "".join(f'<li><a href="#s{i}">{t}</a></li>\n' for i, (t, _) in enumerate(html))
    body = "".join(f'<div style="background:#fff;border:1px solid #e0e0e0;border-radius:6px;'
                   f'padding:18px 22px;margin-bottom:18px" id="s{i}"><h2>{t}</h2>\n{c}\n</div>\n'
                   for i, (t, c) in enumerate(html))
    report = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Survival Validation Report</title>
<style>
  body {{ font-family:Helvetica,Arial,sans-serif; font-size:14px; line-height:1.5;
         max-width:1100px; margin:0 auto; padding:20px 30px; background:#fafafa; color:#333; }}
  h1 {{ font-size:22px; border-bottom:3px solid #0072B2; padding-bottom:8px; color:#0072B2; }}
  h2 {{ font-size:17px; margin:30px 0 12px; border-bottom:1px solid #ddd; padding-bottom:4px; }}
  h3 {{ font-size:14px; margin:18px 0 8px; color:#555; }}
  table {{ border-collapse:collapse; width:100%; margin:10px 0; font-size:13px; }}
  th,td {{ border:1px solid #ddd; padding:6px 10px; text-align:left; }}
  th {{ background:#f5f5f5; font-weight:600; }}
  img {{ max-width:100%; border:1px solid #eee; border-radius:4px; margin:8px 0; }}
  .toc {{ background:#fff; border:1px solid #e0e0e0; border-radius:6px; padding:12px 20px; margin-bottom:18px; }}
  .toc ul {{ list-style:none; padding:0; columns:2; }} .toc a {{ color:#0072B2; text-decoration:none; }}
</style></head><body>
<h1>Survival Signal Validation</h1>
<p style="color:#666;font-size:13px">HGSC malignant-states atlas &middot; {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
<div class="toc"><strong>Contents</strong><ul>{toc}</ul></div>
{body}
</body></html>"""
    with open(os.path.join(OUT_DIR, "22f_validation_report.html"), "w") as f:
        f.write(report)
    print("\n    Saved: 22f_validation_report.html")
    print(f"\n{'='*65}\n  Step complete!\n{'='*65}")


if __name__ == "__main__":
    main()
