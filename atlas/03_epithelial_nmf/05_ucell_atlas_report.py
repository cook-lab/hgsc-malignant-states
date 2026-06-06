#!/usr/bin/env python3
"""
Atlas UCell scoring report & cutoff analysis
============================================
HGSC malignant-states atlas backend.

Summarises the cross-platform UCell scores (from 04_ucell_atlas_scoring.R)
across the NMF-defined SecA / Intermediate / SecB epithelial populations:
distributions, per-label violins, SecA-vs-SecB scatter, candidate polarization
cutoffs, per-patient consistency, and a ROC/Youden cutoff analysis used to
guide organoid classification. Emits a self-contained HTML report.

INPUTS (output_root/03_epithelial_nmf/ucell_atlas/):
  - atlas_ucell_scores.csv
  - atlas_secretory_metadata.csv

OUTPUTS (output_root/03_epithelial_nmf/ucell_atlas/):
  - 18c_ucell_atlas_report.html
  - atlas_ucell_summary_stats.csv
  - 18c_A..F_*.svg/pdf

MANUSCRIPT PANELS: supporting analysis for the cross-platform polarization
  scoring (Fig 3B / SF11 cutoff derivation).

RUNTIME TIER: fast (CSV summarisation + plotting).

SEEDING: plot subsamples use random_state=SEED for determinism.

Usage:
    python 05_ucell_atlas_report.py
"""

import os
import sys
import warnings
import base64
from io import BytesIO

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
OUT_DIR    = path("output_root", "03_epithelial_nmf", "ucell_atlas")
SCORES_CSV = os.path.join(OUT_DIR, "atlas_ucell_scores.csv")
META_CSV   = os.path.join(OUT_DIR, "atlas_secretory_metadata.csv")

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

# NMF label colors ("Transitioning" -> "Intermediate")
NMF_COLORS = {
    "SecA epithelium":          "#0072B2",
    "Intermediate epithelium":  "#999999",
    "SecB epithelium":          "#D55E00",
}
NMF_ORDER = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]

# ── 1. Load data ─────────────────────────────────────────────
print("=" * 60)
print("  Atlas UCell report & cutoff analysis")
print("=" * 60)

print("\n[1] Loading scores and metadata...")
scores = pd.read_csv(SCORES_CSV)
meta   = pd.read_csv(META_CSV, index_col=0)
df = meta.join(scores.set_index("barcode"), how="inner")
print(f"    Merged data: {len(df):,} cells")
for label in NMF_ORDER:
    print(f"    {label}: {(df['celltype_nmf'] == label).sum():,}")

# ── 2. Summary statistics ────────────────────────────────────
print("\n[2] Computing summary statistics...")
stats_rows = []
for label in NMF_ORDER:
    sub = df[df["celltype_nmf"] == label]
    for score_col, score_name in [("SecA_UCell", "SecA_UCell"),
                                   ("SecB_UCell", "SecB_UCell"),
                                   ("sec_polarization", "Polarization")]:
        vals = sub[score_col]
        stats_rows.append({
            "celltype_nmf": label, "score": score_name, "n": len(vals),
            "mean": vals.mean(), "median": vals.median(), "sd": vals.std(),
            "p5": np.percentile(vals, 5), "p25": np.percentile(vals, 25),
            "p75": np.percentile(vals, 75), "p95": np.percentile(vals, 95),
            "min": vals.min(), "max": vals.max(),
            "iqr": np.percentile(vals, 75) - np.percentile(vals, 25),
        })
stats_df = pd.DataFrame(stats_rows)
stats_path = os.path.join(OUT_DIR, "atlas_ucell_summary_stats.csv")
stats_df.to_csv(stats_path, index=False, float_format="%.6f")
print("    Saved: atlas_ucell_summary_stats.csv")


# ── Figure helpers ───────────────────────────────────────────
def save_fig(fig, name):
    fig.savefig(os.path.join(OUT_DIR, f"{name}.svg"), bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{name}.pdf"), bbox_inches="tight")
    plt.close(fig)


def fig_to_base64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", dpi=150)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()
    return b64


print("\n[3] Generating figures...")

# ── A: Score distributions ───────────────────────────────────
fig_a, axes_a = plt.subplots(1, 3, figsize=(12, 3.5))
for ax, (col, title) in zip(axes_a, [("SecA_UCell", "SecA UCell Score"),
                                      ("SecB_UCell", "SecB UCell Score"),
                                      ("sec_polarization", "Polarization (SecB - SecA)")]):
    for label in NMF_ORDER:
        sub = df[df["celltype_nmf"] == label][col]
        ax.hist(sub, bins=80, density=True, alpha=0.5,
                color=NMF_COLORS[label], label=label.replace(" epithelium", ""))
        from scipy.stats import gaussian_kde
        if sub.std() > 0:
            kde = gaussian_kde(sub, bw_method=0.05)
            x = np.linspace(sub.min(), sub.max(), 300)
            ax.plot(x, kde(x), color=NMF_COLORS[label], lw=1.5)
    ax.set_xlabel(title); ax.set_ylabel("Density"); ax.legend(fontsize=6)
fig_a.suptitle("UCell Score Distributions by NMF Label", fontsize=10, weight="bold")
fig_a.tight_layout()
b64_a = fig_to_base64(fig_a)
save_fig(fig_a, "18c_A_score_distributions")
print("    A: Score distributions")

# ── B: Violin + box ──────────────────────────────────────────
fig_b, axes_b = plt.subplots(1, 3, figsize=(12, 4))
for ax, (col, title) in zip(axes_b, [("SecA_UCell", "SecA UCell"),
                                      ("SecB_UCell", "SecB UCell"),
                                      ("sec_polarization", "Polarization")]):
    palette = [NMF_COLORS[l] for l in NMF_ORDER]
    parts = ax.violinplot([df[df["celltype_nmf"] == l][col].values for l in NMF_ORDER],
                          showmeans=False, showmedians=False, showextrema=False)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(palette[i]); pc.set_alpha(0.6)
    bp = ax.boxplot([df[df["celltype_nmf"] == l][col].values for l in NMF_ORDER],
                    widths=0.15, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.5))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(palette[i]); patch.set_alpha(0.8)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["SecA", "Interm.", "SecB"], fontsize=7)
    ax.set_ylabel(title)
    for i, label in enumerate(NMF_ORDER):
        med = df[df["celltype_nmf"] == label][col].median()
        ax.text(i + 1, med, f" {med:.3f}", va="center", ha="left", fontsize=5.5,
                color="black", weight="bold")
fig_b.suptitle("UCell Score Violins by NMF Label", fontsize=10, weight="bold")
fig_b.tight_layout()
b64_b = fig_to_base64(fig_b)
save_fig(fig_b, "18c_B_score_violins")
print("    B: Violin plots")

# ── C: SecA vs SecB scatter ──────────────────────────────────
fig_c, ax_c = plt.subplots(1, 1, figsize=(6, 5.5))
for label in reversed(NMF_ORDER):
    sub = df[df["celltype_nmf"] == label].sample(frac=1, random_state=SEED)
    n_plot = min(len(sub), 30000)
    sub_plot = sub.sample(n=n_plot, random_state=SEED) if len(sub) > n_plot else sub
    ax_c.scatter(sub_plot["SecA_UCell"], sub_plot["SecB_UCell"],
                 c=NMF_COLORS[label], s=0.3, alpha=0.3, rasterized=True,
                 label=f"{label.replace(' epithelium', '')} (n={len(sub):,})")
ax_c.plot([0, 1], [0, 1], "k--", lw=0.8, alpha=0.5)
ax_c.set_xlabel("SecA UCell Score"); ax_c.set_ylabel("SecB UCell Score")
ax_c.set_title("SecA vs SecB UCell Scores (Atlas Secretory Cells)")
ax_c.legend(fontsize=6, markerscale=8)
fig_c.tight_layout()
b64_c = fig_to_base64(fig_c)
save_fig(fig_c, "18c_C_secA_vs_secB_scatter")
print("    C: SecA vs SecB scatter")

# ── D: Polarization distribution with candidate cutoffs ──────
fig_d, ax_d = plt.subplots(1, 1, figsize=(8, 4))
for label in NMF_ORDER:
    sub = df[df["celltype_nmf"] == label]["sec_polarization"]
    ax_d.hist(sub, bins=100, density=True, alpha=0.5,
              color=NMF_COLORS[label], label=label.replace(" epithelium", ""))

interm = df[df["celltype_nmf"] == "Intermediate epithelium"]["sec_polarization"]
secA_vals = df[df["celltype_nmf"] == "SecA epithelium"]["sec_polarization"]
secB_vals = df[df["celltype_nmf"] == "SecB epithelium"]["sec_polarization"]
cut_low = (np.percentile(secA_vals, 75) + np.percentile(interm, 25)) / 2
cut_high = (np.percentile(interm, 75) + np.percentile(secB_vals, 25)) / 2

ax_d.axvline(cut_low, color="#0072B2", ls="--", lw=1.2,
             label=f"Candidate SecA/Interm cutoff: {cut_low:.4f}")
ax_d.axvline(cut_high, color="#D55E00", ls="--", lw=1.2,
             label=f"Candidate Interm/SecB cutoff: {cut_high:.4f}")
ax_d.axvline(0, color="black", ls=":", lw=0.8, alpha=0.5)
ax_d.set_xlabel("Polarization Score (SecB - SecA)"); ax_d.set_ylabel("Density")
ax_d.set_title("Polarization Score Distributions with Candidate Cutoffs")
ax_d.legend(fontsize=6)
fig_d.tight_layout()
b64_d = fig_to_base64(fig_d)
save_fig(fig_d, "18c_D_polarization_by_subtype")
print("    D: Polarization distributions")

# ── E: Per-patient consistency ───────────────────────────────
if "patient_id" in df.columns:
    patient_stats = df.groupby(["patient_id", "celltype_nmf"]).agg(
        median_SecA=("SecA_UCell", "median"),
        median_SecB=("SecB_UCell", "median"),
        median_pol=("sec_polarization", "median"),
        n_cells=("SecA_UCell", "count"),
    ).reset_index()
    min_cells = 50
    patient_stats_filt = patient_stats[patient_stats["n_cells"] >= min_cells]
    fig_e, axes_e = plt.subplots(1, 3, figsize=(14, 4))
    for ax, (col, title) in zip(axes_e, [("median_SecA", "Median SecA UCell"),
                                          ("median_SecB", "Median SecB UCell"),
                                          ("median_pol", "Median Polarization")]):
        for label in NMF_ORDER:
            sub = patient_stats_filt[patient_stats_filt["celltype_nmf"] == label]
            ax.scatter(range(len(sub)), sub[col].sort_values().values,
                       c=NMF_COLORS[label], s=12, alpha=0.7,
                       label=label.replace(" epithelium", ""))
        ax.set_xlabel("Patient (ranked)"); ax.set_ylabel(title); ax.legend(fontsize=5)
    fig_e.suptitle(f"Per-Patient Score Consistency (n>={min_cells} cells/patient/label)",
                   fontsize=10, weight="bold")
    fig_e.tight_layout()
    b64_e = fig_to_base64(fig_e)
    save_fig(fig_e, "18c_E_per_patient")
    print("    E: Per-patient consistency")
else:
    b64_e = None
    print("    E: Skipped (no patient_id)")

# ── F: ROC-like cutoff analysis ──────────────────────────────
secA_secB = df[df["celltype_nmf"].isin(["SecA epithelium", "SecB epithelium"])].copy()
secA_secB["is_SecB"] = (secA_secB["celltype_nmf"] == "SecB epithelium").astype(int)
fpr, tpr, thresholds = roc_curve(secA_secB["is_SecB"], secA_secB["sec_polarization"])
roc_auc = auc(fpr, tpr)
optimal_idx = np.argmax(tpr - fpr)
optimal_threshold = thresholds[optimal_idx]

fig_f, (ax_roc, ax_thresh) = plt.subplots(1, 2, figsize=(11, 4.5))
ax_roc.plot(fpr, tpr, color="#D55E00", lw=2, label=f"AUC = {roc_auc:.4f}")
ax_roc.plot([0, 1], [0, 1], "k--", lw=0.8)
ax_roc.scatter(fpr[optimal_idx], tpr[optimal_idx], c="red", s=60, zorder=5,
               label=f"Optimal: {optimal_threshold:.4f}")
ax_roc.set_xlabel("False Positive Rate"); ax_roc.set_ylabel("True Positive Rate")
ax_roc.set_title("ROC: SecA vs SecB Classification by Polarization")
ax_roc.legend(fontsize=7)

thresh_range = np.linspace(secA_secB["sec_polarization"].quantile(0.01),
                           secA_secB["sec_polarization"].quantile(0.99), 200)
sensitivities, specificities = [], []
for t in thresh_range:
    tp = ((secA_secB["sec_polarization"] >= t) & (secA_secB["is_SecB"] == 1)).sum()
    fn = ((secA_secB["sec_polarization"] < t) & (secA_secB["is_SecB"] == 1)).sum()
    tn = ((secA_secB["sec_polarization"] < t) & (secA_secB["is_SecB"] == 0)).sum()
    fp = ((secA_secB["sec_polarization"] >= t) & (secA_secB["is_SecB"] == 0)).sum()
    sensitivities.append(tp / (tp + fn) if (tp + fn) > 0 else 0)
    specificities.append(tn / (tn + fp) if (tn + fp) > 0 else 0)
ax_thresh.plot(thresh_range, sensitivities, color="#D55E00", lw=1.5, label="Sensitivity (SecB)")
ax_thresh.plot(thresh_range, specificities, color="#0072B2", lw=1.5, label="Specificity (SecA)")
ax_thresh.axvline(optimal_threshold, color="red", ls="--", lw=1,
                  label=f"Optimal: {optimal_threshold:.4f}")
ax_thresh.set_xlabel("Polarization Cutoff"); ax_thresh.set_ylabel("Rate")
ax_thresh.set_title("Sensitivity/Specificity vs Polarization Threshold")
ax_thresh.legend(fontsize=7)
fig_f.tight_layout()
b64_f = fig_to_base64(fig_f)
save_fig(fig_f, "18c_F_cutoff_roc")
print("    F: ROC & threshold analysis")

# ── Confusion matrix at candidate cutoffs ────────────────────
df["predicted"] = "Intermediate"
df.loc[df["sec_polarization"] < cut_low, "predicted"] = "SecA"
df.loc[df["sec_polarization"] >= cut_high, "predicted"] = "SecB"

pred_order = ["SecA", "Intermediate", "SecB"]
true_labels = [l.replace(" epithelium", "") for l in NMF_ORDER]
confusion = pd.crosstab(df["celltype_nmf"].str.replace(" epithelium", ""),
                        df["predicted"], margins=True)
for c in pred_order:
    if c not in confusion.columns:
        confusion[c] = 0
confusion = confusion.reindex(index=true_labels + ["All"], columns=pred_order + ["All"],
                              fill_value=0)
confusion_pct = confusion.copy().astype(float)
for label in true_labels:
    row_total = confusion.loc[label, "All"]
    if row_total > 0:
        for c in pred_order:
            confusion_pct.loc[label, c] = confusion.loc[label, c] / row_total * 100

# ── Build HTML report ────────────────────────────────────────
print("\n[10] Building HTML report...")
stats_html = stats_df.round(4).to_html(index=False, classes="stats-table")
confusion_html = confusion.to_html(classes="stats-table")
confusion_pct_html = confusion_pct.round(1).to_html(classes="stats-table")

html_content = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="UTF-8">
<title>Atlas UCell Scoring Report</title>
<style>
  body {{ font-family: Helvetica, Arial, sans-serif; max-width: 1200px;
          margin: 0 auto; padding: 20px; background: #fafafa; color: #333; }}
  h1 {{ color: #222; border-bottom: 3px solid #0072B2; padding-bottom: 8px; }}
  h2 {{ color: #0072B2; margin-top: 40px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }}
  h3 {{ color: #555; }}
  .figure {{ text-align: center; margin: 20px 0; background: white; padding: 15px;
             border-radius: 6px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
  .figure img {{ max-width: 100%; height: auto; }}
  .stats-table {{ border-collapse: collapse; margin: 15px auto; font-size: 12px; }}
  .stats-table th, .stats-table td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: right; }}
  .stats-table th {{ background: #0072B2; color: white; font-weight: bold; }}
  .stats-table tr:nth-child(even) {{ background: #f5f5f5; }}
  .highlight {{ background: #fff3cd; padding: 12px; border-left: 4px solid #D55E00;
                margin: 15px 0; border-radius: 4px; }}
  .key-value {{ display: inline-block; background: #e8f4fd; padding: 4px 10px;
                border-radius: 4px; margin: 3px; font-size: 12px; }}
  .key-value strong {{ color: #0072B2; }}
  .two-col {{ display: flex; gap: 20px; flex-wrap: wrap; }}
  .two-col > div {{ flex: 1; min-width: 400px; }}
</style></head><body>

<h1>Atlas UCell Scoring Report</h1>
<p>UCell scores for atlas secretory epithelial cells (NMF-defined
SecA / Intermediate / SecB), scored with the gene-space-aligned noBCAM
signatures matching the organoid/xenium pipeline.</p>

<div class="highlight">
  <strong>Key numbers:</strong><br>
  <span class="key-value"><strong>Total cells:</strong> {len(df):,}</span>
  <span class="key-value"><strong>SecA:</strong> {(df["celltype_nmf"]=="SecA epithelium").sum():,}</span>
  <span class="key-value"><strong>Intermediate:</strong> {(df["celltype_nmf"]=="Intermediate epithelium").sum():,}</span>
  <span class="key-value"><strong>SecB:</strong> {(df["celltype_nmf"]=="SecB epithelium").sum():,}</span>
  <span class="key-value"><strong>ROC AUC (SecA vs SecB):</strong> {roc_auc:.4f}</span>
  <span class="key-value"><strong>Optimal polarization cutoff:</strong> {optimal_threshold:.4f}</span>
</div>

<h2>1. Score Distributions by NMF Label</h2>
<div class="figure"><img src="data:image/png;base64,{b64_a}" alt="Score distributions"></div>
<h2>2. Violin + Box Plots</h2>
<div class="figure"><img src="data:image/png;base64,{b64_b}" alt="Violin plots"></div>
<h2>3. SecA vs SecB Scatter</h2>
<div class="figure"><img src="data:image/png;base64,{b64_c}" alt="Scatter plot"></div>
<h2>4. Polarization Distributions with Candidate Cutoffs</h2>
<div class="figure"><img src="data:image/png;base64,{b64_d}" alt="Polarization distributions"></div>
<div class="highlight">
  <strong>Candidate cutoffs (percentile midpoint method):</strong><br>
  <span class="key-value"><strong>SecA/Intermediate:</strong> {cut_low:.4f}</span>
  <span class="key-value"><strong>Intermediate/SecB:</strong> {cut_high:.4f}</span>
</div>
<h2>5. Per-Patient Consistency</h2>
{"<div class='figure'><img src='data:image/png;base64," + b64_e + "' alt='Per-patient'></div>" if b64_e else "<p><em>Patient information not available.</em></p>"}
<h2>6. ROC & Threshold Analysis (SecA vs SecB)</h2>
<div class="figure"><img src="data:image/png;base64,{b64_f}" alt="ROC and threshold"></div>
<div class="highlight">
  <strong>ROC-optimal cutoff (Youden's J):</strong>
  <span class="key-value"><strong>Polarization threshold:</strong> {optimal_threshold:.4f}</span>
  <span class="key-value"><strong>AUC:</strong> {roc_auc:.4f}</span>
</div>
<h2>7. Classification at Candidate Cutoffs</h2>
<div class="two-col">
  <div><h3>Confusion Matrix (counts)</h3>{confusion_html}</div>
  <div><h3>Confusion Matrix (% of true label)</h3>{confusion_pct_html}</div>
</div>
<h2>8. Summary Statistics</h2>
{stats_html}
</body></html>
"""

report_path = os.path.join(OUT_DIR, "18c_ucell_atlas_report.html")
with open(report_path, "w") as f:
    f.write(html_content)
print(f"    Saved: {report_path}")

print(f"\n{'='*60}")
print("  Step complete!")
print(f"  Report: {report_path}")
print(f"{'='*60}")
