#!/usr/bin/env python3
"""
SecA vs SecB differential communication (NMF labels)
====================================================
HGSC malignant-states atlas backend.

Analyses the global LIANA+ results from 01_cellcomm_nmf.py to identify
differential ligand-receptor interactions between SecA and SecB epithelial
cells (incoming + outgoing), differential TME communication partners, and the
SecA-enriched / SecB-enriched / shared L-R pair lists. Emits CSV tables, bar
figures, and an HTML report. Curated autocrine merge -> Supp Data 7.

INPUTS:
  - output_root/06_cellcomm/tables/17b_liana_global.csv  (from 01_cellcomm_nmf.py)

OUTPUTS (output_root/06_cellcomm/17c_secA_secB_communication_nmf/):
  - tables/17c_*.csv
  - figs/17c_*.svg/pdf
  - 17c_secA_secB_communication_nmf.html

MANUSCRIPT PANELS: Fig 5F (autocrine shift), Supp Data 7 (autocrine L-R pairs).

RUNTIME TIER: fast (post-hoc analysis of the LIANA table).

SEEDING: deterministic (no stochastic step).

Usage:
    python 02_secA_secB_communication_nmf.py
"""

import os
import sys
import base64
import io
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

warnings.filterwarnings("ignore")

# ── Cook Lab style v1.2 ─────────────────────────────────────
plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":       8,
    "axes.titlesize":  9,
    "axes.labelsize":  8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6,
    "figure.dpi":      150,
    "savefig.dpi":     450,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "svg.fonttype":    "none",
    "savefig.bbox":    "tight",
    "axes.spines.top":   False,
    "axes.spines.right": False,
})

# ── Paths ────────────────────────────────────────────────────
LIANA_CSV = path("output_root", "06_cellcomm", "tables", "17b_liana_global.csv")
OUT_DIR   = path("output_root", "06_cellcomm", "17c_secA_secB_communication_nmf")
FIG_DIR   = os.path.join(OUT_DIR, "figs")
TABLE_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# ── NMF-based epithelial schema ("Transitioning" -> "Intermediate") ──
V2_MAP = {
    "SecA epithelium":          "SecA",
    "Intermediate epithelium":  "Intermediate",
    "SecB epithelium":          "SecB",
    "Ciliated epithelial cell": "Ciliated",
}
V2_PALETTE = {
    "SecA":         "#E6A141", "SecB":         "#B8741A",
    "Intermediate": "#C49A5E", "Ciliated":     "#E05A2C",
}

# ── Level-1 compartment mapping (NMF epithelial + level2 non-epithelial) ──
LEVEL1_MAP = {
    "SecA epithelium":          "Epithelial",
    "Intermediate epithelium":  "Epithelial",
    "SecB epithelium":          "Epithelial",
    "Ciliated epithelial cell": "Epithelial",
    "Mesothelial cell":                          "Mesothelial",
    "Hypoxic mesothelial cell":                  "Mesothelial",
    "Activated fibroblast":                      "Fibroblast",
    "Myo-fibroblastic cancer-associated fibroblast": "Fibroblast",
    "PI16+ universal fibroblast":                "Fibroblast",
    "Cycling fibroblast":                        "Fibroblast",
    "Hypoxic inflammatory cancer-associated fibroblast": "Fibroblast",
    "Ovarian stromal cell":                      "Fibroblast",
    "Ovarian steroidogenic cell":                "Fibroblast",
    "Schwann cell":                              "Fibroblast",
    "Contractile smooth muscle cell":            "Smooth muscle",
    "Stress-response smooth muscle cell":        "Smooth muscle",
    "Inflammatory fibroblast-like smooth muscle cell": "Smooth muscle",
    "Pericyte":                                  "Pericyte",
    "Angiogenic endothelial cell":               "Endothelial",
    "Venous endothelial cell":                   "Endothelial",
    "Arterial endothelial cell":                 "Endothelial",
    "Lymphatic endothelial cell":                "Endothelial",
    "Cycling endothelial cell":                  "Endothelial",
    "CD8 effector/exhausted T cell":             "T/NK cell",
    "CD8 effector T cell":                       "T/NK cell",
    "CD8 tissue-resident memory T cell":         "T/NK cell",
    "CD4 naive T cell":                          "T/NK cell",
    "CD4 Regulatory T cell":                     "T/NK cell",
    "CD4 regulatory T cell":                     "T/NK cell",
    "Cycling T/NK cell":                         "T/NK cell",
    "Quiescent T cell":                          "T/NK cell",
    "MAIT cell":                                 "T/NK cell",
    "NK cell":                                   "T/NK cell",
    "CD56bright NK cell":                        "T/NK cell",
    "CD56dim NK cell":                           "T/NK cell",
    "Metallothionein-high stress-response T cell": "T/NK cell",
    "Metallothionein-stress T cell":             "T/NK cell",
    "gdT cell":                                  "T/NK cell",
    "IFN-activated B cell":                      "B cell",
    "Activated B cell":                          "B cell",
    "Cycling B cell":                            "B cell",
    "Plasma cell":                               "Plasma cell",
    "Cycling plasma cell":                       "Plasma cell",
    "C1Q tissue-resident macrophage":            "Macrophage",
    "Cycling C1Q+ tissue-resident macrophage":   "Macrophage",
    "Monocyte-derived macrophage":               "Macrophage",
    "Inflammatory macrophage":                   "Macrophage",
    "Hypoxic macrophage":                        "Macrophage",
    "Cycling macrophage":                        "Macrophage",
    "Monocyte":                                  "Macrophage",
    "Classical monocyte":                        "Macrophage",
    "Type 1 DC":                                 "DC",
    "Type 2 DC":                                 "DC",
    "Conventional dendritic cell type 1":        "DC",
    "Conventional dendritic cell type 2":        "DC",
    "Mature DC":                                 "DC",
    "Mature dendritic cell":                     "DC",
    "Plasmacytoid DC":                           "DC",
    "Plasmacytoid dendritic cell":               "DC",
    "Neutrophil":                                "Neutrophil",
    "Mast cell":                                 "Mast cell",
    "Cycling mast cell":                         "Mast cell",
    "Hematopoietic stem cell":                   "Other",
}
COMPARTMENT_ORDER = [
    "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
    "T/NK cell", "B cell", "Plasma cell", "Macrophage", "DC", "Neutrophil", "Mast cell",
]
SIG_THRESH = 0.05
EPI_LABELS = {"SecA", "SecB", "Intermediate", "Ciliated"}

# ── 1. Load & annotate ───────────────────────────────────────
print("=" * 60)
print("  SecA vs SecB Communication (NMF labels)")
print("=" * 60)

print("\n[1] Loading global LIANA results (NMF labels)...")
df = pd.read_csv(LIANA_CSV)
print(f"    {len(df):,} total interactions")

df["source_v2"]     = df["source"].map(V2_MAP)
df["target_v2"]     = df["target"].map(V2_MAP)
df["source_level1"] = df["source"].map(LEVEL1_MAP).fillna("Unknown")
df["target_level1"] = df["target"].map(LEVEL1_MAP).fillna("Unknown")

sig = df[df["magnitude_rank"] <= SIG_THRESH].copy()
print(f"    {len(sig):,} significant (magnitude_rank <= {SIG_THRESH})")
sig["lr_pair"] = sig["ligand_complex"] + " -> " + sig["receptor_complex"]

# ── 2. SecA vs SecB differential L-R pairs ───────────────────
print("\n[2] Computing differential L-R pairs (SecA vs SecB)...")

def compute_pole_lr_stats(sig_df, pole, direction):
    sub = sig_df[sig_df["target_v2"] == pole] if direction == "incoming" \
        else sig_df[sig_df["source_v2"] == pole]
    return sub.groupby("lr_pair").agg(
        n_interactions=("lr_pair", "size"),
        mean_magnitude_rank=("magnitude_rank", "mean"),
        mean_lr_means=("lr_means", "mean"),
        mean_lrscore=("lrscore", "mean"),
    ).reset_index()

diff_rows = []
for direction in ["incoming", "outgoing"]:
    seca_stats = compute_pole_lr_stats(sig, "SecA", direction)
    secb_stats = compute_pole_lr_stats(sig, "SecB", direction)
    merged = seca_stats.merge(secb_stats, on="lr_pair", how="outer",
                              suffixes=("_SecA", "_SecB")).fillna(0)
    merged["n_SecA"] = merged["n_interactions_SecA"]
    merged["n_SecB"] = merged["n_interactions_SecB"]
    merged["log2fc"] = np.log2((merged["n_SecB"] + 1) / (merged["n_SecA"] + 1))
    merged["delta_lrscore"] = merged["mean_lrscore_SecB"] - merged["mean_lrscore_SecA"]
    merged["direction"] = direction
    diff_rows.append(merged)

diff_lr = pd.concat(diff_rows, ignore_index=True)
diff_lr["abs_log2fc"] = diff_lr["log2fc"].abs()

for d in ["incoming", "outgoing"]:
    sub = diff_lr[diff_lr["direction"] == d]
    sub = sub[(sub["n_SecA"] >= 3) | (sub["n_SecB"] >= 3)]
    sub.nlargest(30, "abs_log2fc").to_csv(
        os.path.join(TABLE_DIR, f"17c_differential_lr_{d}_top30.csv"), index=False)
    print(f"    {d}: {len(sub)} L-R pairs, saved top 30")
diff_lr.to_csv(os.path.join(TABLE_DIR, "17c_differential_lr_all.csv"), index=False)

# ── 3. Differential TME communication partners ───────────────
print("\n[3] Computing differential TME partners...")
partner_rows = []
for direction in ["incoming", "outgoing"]:
    if direction == "incoming":
        sub_a = sig[(sig["target_v2"] == "SecA") & (~sig["source_v2"].isin(EPI_LABELS))]
        sub_b = sig[(sig["target_v2"] == "SecB") & (~sig["source_v2"].isin(EPI_LABELS))]
        partner_col = "source_level1"
    else:
        sub_a = sig[(sig["source_v2"] == "SecA") & (~sig["target_v2"].isin(EPI_LABELS))]
        sub_b = sig[(sig["source_v2"] == "SecB") & (~sig["target_v2"].isin(EPI_LABELS))]
        partner_col = "target_level1"
    count_a = sub_a[partner_col].value_counts().rename("n_SecA")
    count_b = sub_b[partner_col].value_counts().rename("n_SecB")
    merged = pd.concat([count_a, count_b], axis=1).fillna(0).astype(int)
    merged["prop_SecA"] = merged["n_SecA"] / max(merged["n_SecA"].sum(), 1) * 100
    merged["prop_SecB"] = merged["n_SecB"] / max(merged["n_SecB"].sum(), 1) * 100
    merged["delta_prop"] = merged["prop_SecB"] - merged["prop_SecA"]
    merged["log2fc"] = np.log2((merged["n_SecB"] + 1) / (merged["n_SecA"] + 1))
    merged["direction"] = direction
    merged.index.name = "compartment"
    partner_rows.append(merged.reset_index())
partner_df = pd.concat(partner_rows, ignore_index=True)
partner_df = partner_df[partner_df["compartment"].isin(COMPARTMENT_ORDER)]
partner_df.to_csv(os.path.join(TABLE_DIR, "17c_differential_tme_partners.csv"), index=False)
print(f"    Saved: {len(partner_df)} compartment x direction entries")

# ── 4. Top L-R pairs per pole ────────────────────────────────
print("\n[4] Identifying SecA-unique, SecB-unique, and shared top L-R pairs...")
for direction in ["incoming", "outgoing"]:
    sub = diff_lr[diff_lr["direction"] == direction].copy()
    sub = sub[(sub["n_SecA"] >= 3) | (sub["n_SecB"] >= 3)]
    seca_enriched = sub[sub["log2fc"] < -1].nsmallest(20, "log2fc")
    secb_enriched = sub[sub["log2fc"] > 1].nlargest(20, "log2fc")
    shared = sub[(sub["log2fc"].abs() < 0.5) &
                 (sub["n_SecA"] >= 5) & (sub["n_SecB"] >= 5)].copy()
    shared["_total"] = shared["n_SecA"] + shared["n_SecB"]
    shared = shared.nlargest(20, "_total").drop(columns=["_total"])
    seca_enriched.to_csv(os.path.join(TABLE_DIR, f"17c_seca_enriched_lr_{direction}.csv"), index=False)
    secb_enriched.to_csv(os.path.join(TABLE_DIR, f"17c_secb_enriched_lr_{direction}.csv"), index=False)
    shared.to_csv(os.path.join(TABLE_DIR, f"17c_shared_lr_{direction}.csv"), index=False)
    print(f"    {direction}: {len(seca_enriched)} SecA-enriched, "
          f"{len(secb_enriched)} SecB-enriched, {len(shared)} shared")

# ── 5. Interaction counts summary ────────────────────────────
print("\n[5] Computing interaction count summary...")
count_summary = {}
for pole in ["SecA", "SecB", "Intermediate", "Ciliated"]:
    n_inc = len(sig[sig["target_v2"] == pole])
    n_out = len(sig[sig["source_v2"] == pole])
    count_summary[pole] = {"incoming": n_inc, "outgoing": n_out, "total": n_inc + n_out}
count_df = pd.DataFrame(count_summary).T
count_df.index.name = "nmf_group"
count_df.to_csv(os.path.join(TABLE_DIR, "17c_interaction_counts_by_nmf_group.csv"))
print(f"    Interaction counts:\n{count_df}")

# ── 6. Figures ───────────────────────────────────────────────
print("\n[6] Generating figures...")
fig_data = {}

def fig_to_base64(fig, dpi=150):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("utf-8")
    plt.close(fig)
    return b64

# 6a. Differential L-R bar plots
for direction in ["incoming", "outgoing"]:
    sub = diff_lr[diff_lr["direction"] == direction].copy()
    sub = sub[(sub["n_SecA"] >= 3) | (sub["n_SecB"] >= 3)]
    plot_df = pd.concat([sub.nsmallest(15, "log2fc"), sub.nlargest(15, "log2fc")]) \
        .drop_duplicates("lr_pair").sort_values("log2fc")
    fig, ax = plt.subplots(figsize=(6, 8))
    colors = [V2_PALETTE["SecA"] if v < 0 else V2_PALETTE["SecB"] for v in plot_df["log2fc"]]
    ax.barh(range(len(plot_df)), plot_df["log2fc"].values, color=colors,
            edgecolor="#333333", linewidth=0.3)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["lr_pair"].values, fontsize=6)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("log2FC (SecB / SecA)", fontsize=8)
    ax.set_title(f"Differential L-R pairs — {direction} (NMF labels)",
                 fontsize=10, fontweight="bold")
    for i, (_, row) in enumerate(plot_df.iterrows()):
        label = f"A:{int(row['n_SecA'])} B:{int(row['n_SecB'])}"
        x_pos = row["log2fc"]
        ha = "left" if x_pos >= 0 else "right"
        ax.text(x_pos + (0.05 if x_pos >= 0 else -0.05), i, label,
                va="center", ha=ha, fontsize=4.5, color="#666666")
    ax.legend(handles=[Patch(facecolor=V2_PALETTE["SecA"], label="SecA-enriched"),
                       Patch(facecolor=V2_PALETTE["SecB"], label="SecB-enriched")],
              fontsize=7, frameon=False, loc="lower right")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, f"17c_diff_lr_{direction}.svg"), format="svg",
                bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, f"17c_diff_lr_{direction}.pdf"), format="pdf",
                bbox_inches="tight")
    fig_data[f"diff_lr_{direction}"] = fig_to_base64(fig)

# 6b. TME partner bars
for direction in ["incoming", "outgoing"]:
    sub = partner_df[partner_df["direction"] == direction].copy()
    sub = sub[sub["compartment"].isin(COMPARTMENT_ORDER)]
    sub["compartment"] = pd.Categorical(sub["compartment"],
                                        categories=COMPARTMENT_ORDER, ordered=True)
    sub = sub.sort_values("compartment")
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(sub))
    w = 0.35
    ax.bar(x - w/2, sub["prop_SecA"].values, w, color=V2_PALETTE["SecA"],
           edgecolor="#333333", linewidth=0.3, label="SecA")
    ax.bar(x + w/2, sub["prop_SecB"].values, w, color=V2_PALETTE["SecB"],
           edgecolor="#333333", linewidth=0.3, label="SecB")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["compartment"].values, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("% of pole's interactions", fontsize=8)
    ax.set_title(f"TME partner distribution — {direction} (NMF labels)",
                 fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, f"17c_tme_partners_{direction}.svg"), format="svg",
                bbox_inches="tight")
    fig_data[f"tme_partners_{direction}"] = fig_to_base64(fig)

# 6c. Interaction count summary bar
fig, ax = plt.subplots(figsize=(5, 3))
groups = ["SecA", "SecB", "Intermediate", "Ciliated"]
x = np.arange(len(groups))
inc_vals = [count_summary.get(g, {}).get("incoming", 0) for g in groups]
out_vals = [count_summary.get(g, {}).get("outgoing", 0) for g in groups]
w = 0.35
ax.bar(x - w/2, inc_vals, w, color=[V2_PALETTE.get(g, "#ccc") for g in groups],
       edgecolor="#333333", linewidth=0.3, label="Incoming", alpha=0.7)
ax.bar(x + w/2, out_vals, w, color=[V2_PALETTE.get(g, "#ccc") for g in groups],
       edgecolor="#333333", linewidth=0.3, label="Outgoing", alpha=1.0)
ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
ax.set_ylabel("Significant interactions", fontsize=8)
ax.set_title("Interaction counts by NMF group", fontsize=10, fontweight="bold")
ax.legend(fontsize=7, frameon=False)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "17c_interaction_counts.svg"), format="svg",
            bbox_inches="tight")
fig_data["interaction_counts"] = fig_to_base64(fig)

# ── 7. HTML report ───────────────────────────────────────────
print("\n[7] Generating HTML report...")

def img_tag(key, width="100%"):
    if key in fig_data:
        return f'<img src="data:image/png;base64,{fig_data[key]}" style="width:{width};">'
    return "<p><em>Figure not generated</em></p>"

html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>SecA vs SecB Cell Communication (NMF labels)</title>
<style>
    body {{ font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto;
           padding: 20px; background: #fafafa; color: #333; }}
    h1 {{ color: #333; border-bottom: 2px solid #E6A141; padding-bottom: 8px; }}
    h2 {{ color: #555; margin-top: 30px; border-bottom: 1px solid #ddd; padding-bottom: 5px; }}
    .fig {{ text-align: center; margin: 15px 0; }}
    .fig img {{ border: 1px solid #ddd; border-radius: 4px; }}
    .schema {{ background: #fff; padding: 12px; border-left: 4px solid #E6A141; margin: 10px 0; }}
    .two-col {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
</style></head><body>

<h1>SecA vs SecB Cell Communication (NMF labels)</h1>
<p>LIANA+ results re-analysed using NMF-derived epithelial labels.
Significance: magnitude_rank &le; {SIG_THRESH}.</p>

<div class="schema">
<strong>NMF Schema:</strong>
SecA = NMF Factor_2 &lt; p50 |
Intermediate = p50 &le; Factor_2 &lt; p75 |
SecB = Factor_2 &ge; p75 |
Ciliated = celltype_level2
</div>

<h2>1. Interaction Counts by NMF Group</h2>
<div class="fig">{img_tag("interaction_counts", "50%")}</div>

<h2>2. Differential L-R Pairs (SecA vs SecB)</h2>
<h3>Incoming signals (TME -> Epithelial)</h3>
<div class="fig">{img_tag("diff_lr_incoming", "70%")}</div>
<h3>Outgoing signals (Epithelial -> TME)</h3>
<div class="fig">{img_tag("diff_lr_outgoing", "70%")}</div>

<h2>3. TME Partner Distribution</h2>
<div class="two-col">
<div><h3>Incoming</h3><div class="fig">{img_tag("tme_partners_incoming", "100%")}</div></div>
<div><h3>Outgoing</h3><div class="fig">{img_tag("tme_partners_outgoing", "100%")}</div></div>
</div>
</body></html>
"""

html_path = os.path.join(OUT_DIR, "17c_secA_secB_communication_nmf.html")
with open(html_path, "w") as f:
    f.write(html)
print(f"    Saved: {html_path}")
print("\nDone!")
