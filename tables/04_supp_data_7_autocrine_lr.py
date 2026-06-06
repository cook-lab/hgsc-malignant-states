#!/usr/bin/env python3
"""
Supp Data 7 — SecA/SecB autocrine ligand-receptor pairs
=======================================================
HGSC malignant-states atlas · supplemental table generator.

PURPOSE
    Re-analyse the global LIANA+ result (canonical step 17b, run on NMF-derived
    epithelial labels) to derive SecA-vs-SecB differential ligand-receptor
    communication, and emit the curated autocrine L-R table for Supp Data 7.
    Analytical logic preserved verbatim from canonical step 17c; this migration
    re-routes paths via config, standardizes the "Transitioning"->"Intermediate"
    label, and adds the curated SD7 autocrine export + tidy supplemental write.

INPUTS  (under output_root)
    - 17_cellcomm_nmf/tables/17b_liana_global.csv   (global LIANA L-R table; from step 17b)

OUTPUTS  (under output_root/17c_secA_secB_communication_nmf/tables/ and supplemental/)
    - 17c_differential_lr_all.csv, 17c_differential_lr_{incoming,outgoing}_top30.csv
    - 17c_{seca,secb}_enriched_lr_{incoming,outgoing}.csv, 17c_shared_lr_{incoming,outgoing}.csv
    - 17c_differential_tme_partners.csv, 17c_interaction_counts_by_nmf_group.csv
    - supplemental/Supplemental_Table_7_autocrine_LR_pairs.csv   <- Supp Data 7 (curated)
    Diagnostic figures (svg/pdf) are also written, unchanged.

MANUSCRIPT PANEL(S)
    Supp Data 7 (autocrine L-R pairs SecA/SecB). The 17b_liana_global cache also
    drives Fig 5F (autocrine dumbbell).

RUNTIME TIER
    fast (operates on the precomputed 17b LIANA table).
"""

import io
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- central config (tables/ is 1 level below repo root) ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
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

# ── Paths (config-resolved) ──────────────────────────────────
LIANA_CSV = path("output_root", "17_cellcomm_nmf", "tables", "17b_liana_global.csv")
OUT_DIR   = path("output_root", "17c_secA_secB_communication_nmf")
FIG_DIR   = os.path.join(OUT_DIR, "figs")
TABLE_DIR = os.path.join(OUT_DIR, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# ── NMF-based epithelial schema ("Intermediate" replaces legacy "Transitioning") ──
V2_MAP = {
    "SecA epithelium":          "SecA",
    "Intermediate epithelium":  "Intermediate",
    "Transitioning epithelium": "Intermediate",   # legacy alias -> standardized label
    "SecB epithelium":          "SecB",
    "Ciliated epithelial cell": "Ciliated",
}

V2_PALETTE = {
    "SecA":         "#E6A141",
    "SecB":         "#B8741A",
    "Intermediate": "#C49A5E",
    "Ciliated":     "#E05A2C",
}

LEVEL1_MAP = {
    # Epithelial (NMF labels)
    "SecA epithelium":          "Epithelial",
    "Intermediate epithelium":  "Epithelial",
    "Transitioning epithelium": "Epithelial",
    "SecB epithelium":          "Epithelial",
    "Ciliated epithelial cell": "Epithelial",
    # Mesothelial
    "Mesothelial cell":                          "Mesothelial",
    "Hypoxic mesothelial cell":                  "Mesothelial",
    # Fibroblast
    "Activated fibroblast":                      "Fibroblast",
    "Myo-fibroblastic cancer-associated fibroblast": "Fibroblast",
    "PI16+ universal fibroblast":                "Fibroblast",
    "Cycling fibroblast":                        "Fibroblast",
    "Hypoxic inflammatory cancer-associated fibroblast": "Fibroblast",
    "Ovarian stromal cell":                      "Fibroblast",
    "Ovarian steroidogenic cell":                "Fibroblast",
    "Schwann cell":                              "Fibroblast",
    # Smooth muscle
    "Contractile smooth muscle cell":            "Smooth muscle",
    "Stress-response smooth muscle cell":        "Smooth muscle",
    "Inflammatory fibroblast-like smooth muscle cell": "Smooth muscle",
    # Pericyte
    "Pericyte":                                  "Pericyte",
    # Endothelial
    "Angiogenic endothelial cell":               "Endothelial",
    "Venous endothelial cell":                   "Endothelial",
    "Arterial endothelial cell":                 "Endothelial",
    "Lymphatic endothelial cell":                "Endothelial",
    "Cycling endothelial cell":                  "Endothelial",
    # T/NK
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
    # B cell
    "IFN-activated B cell":                      "B cell",
    "Activated B cell":                          "B cell",
    "Cycling B cell":                            "B cell",
    # Plasma
    "Plasma cell":                               "Plasma cell",
    "Cycling plasma cell":                       "Plasma cell",
    # Macrophage
    "C1Q tissue-resident macrophage":            "Macrophage",
    "Cycling C1Q+ tissue-resident macrophage":   "Macrophage",
    "Monocyte-derived macrophage":               "Macrophage",
    "Inflammatory macrophage":                   "Macrophage",
    "Hypoxic macrophage":                        "Macrophage",
    "Cycling macrophage":                        "Macrophage",
    "Monocyte":                                  "Macrophage",
    "Classical monocyte":                        "Macrophage",
    # DC
    "Type 1 DC":                                 "DC",
    "Type 2 DC":                                 "DC",
    "Conventional dendritic cell type 1":        "DC",
    "Conventional dendritic cell type 2":        "DC",
    "Mature DC":                                 "DC",
    "Mature dendritic cell":                     "DC",
    "Plasmacytoid DC":                           "DC",
    "Plasmacytoid dendritic cell":               "DC",
    # Neutrophil / Mast
    "Neutrophil":                                "Neutrophil",
    "Mast cell":                                 "Mast cell",
    "Cycling mast cell":                         "Mast cell",
    # Other
    "Hematopoietic stem cell":                   "Other",
}

COMPARTMENT_ORDER = [
    "Mesothelial", "Fibroblast", "Smooth muscle", "Pericyte", "Endothelial",
    "T/NK cell", "B cell", "Plasma cell", "Macrophage", "DC",
    "Neutrophil", "Mast cell",
]

SIG_THRESH = 0.05
EPI_LABELS = {"SecA", "SecB", "Intermediate", "Ciliated"}

# ============================================================================
# 1. LOAD & ANNOTATE
# ============================================================================
print("=" * 60)
print("  Supp Data 7 — SecA vs SecB Communication (NMF labels)")
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

# ============================================================================
# 2. SecA vs SecB DIFFERENTIAL L-R PAIRS
# ============================================================================
print("\n[2] Computing differential L-R pairs (SecA vs SecB)...")

def compute_pole_lr_stats(sig_df, pole, direction):
    if direction == "incoming":
        sub = sig_df[sig_df["target_v2"] == pole]
    else:
        sub = sig_df[sig_df["source_v2"] == pole]
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

# ============================================================================
# 3. DIFFERENTIAL TME COMMUNICATION PARTNERS
# ============================================================================
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

# ============================================================================
# 4. TOP L-R PAIRS PER POLE
# ============================================================================
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

# ============================================================================
# 5. INTERACTION COUNTS SUMMARY
# ============================================================================
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

# ============================================================================
# 6. SUPP DATA 7 — curated autocrine L-R table
# ============================================================================
# Autocrine = ligand source pole and receptor target pole are the SAME epithelial
# pole. Curate SecA-autocrine and SecB-autocrine significant pairs into one table.
print("\n[6] Curating Supp Data 7 autocrine L-R table...")

autocrine_rows = []
for pole in ["SecA", "SecB"]:
    auto = sig[(sig["source_v2"] == pole) & (sig["target_v2"] == pole)].copy()
    grp = auto.groupby(["lr_pair", "ligand_complex", "receptor_complex"]).agg(
        n_interactions=("lr_pair", "size"),
        mean_magnitude_rank=("magnitude_rank", "mean"),
        mean_lr_means=("lr_means", "mean"),
        mean_lrscore=("lrscore", "mean"),
    ).reset_index()
    grp.insert(0, "epithelial_pole", pole)
    autocrine_rows.append(grp)

supp7 = pd.concat(autocrine_rows, ignore_index=True)
supp7 = supp7.sort_values(["epithelial_pole", "mean_lrscore"],
                          ascending=[True, False]).reset_index(drop=True)
supp7_out = path("output_root", "supplemental",
                 "Supplemental_Table_7_autocrine_LR_pairs.csv")
supp7.to_csv(supp7_out, index=False)
print(f"    Saved: {supp7_out}  ({len(supp7):,} autocrine pairs)")

# ============================================================================
# 7. FIGURES  (diagnostic; logic unchanged)
# ============================================================================
print("\n[7] Generating figures...")

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
    ax.set_title(f"Differential L-R pairs - {direction} (NMF labels)", fontsize=10, fontweight="bold")
    for i, (_, row) in enumerate(plot_df.iterrows()):
        x_pos = row["log2fc"]
        ha = "left" if x_pos >= 0 else "right"
        offset = 0.05 if x_pos >= 0 else -0.05
        ax.text(x_pos + offset, i, f"A:{int(row['n_SecA'])} B:{int(row['n_SecB'])}",
                va="center", ha=ha, fontsize=4.5, color="#666666")
    ax.legend(handles=[Patch(facecolor=V2_PALETTE["SecA"], label="SecA-enriched"),
                       Patch(facecolor=V2_PALETTE["SecB"], label="SecB-enriched")],
              fontsize=7, frameon=False, loc="lower right")
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, f"17c_diff_lr_{direction}.svg"), format="svg", bbox_inches="tight")
    fig.savefig(os.path.join(FIG_DIR, f"17c_diff_lr_{direction}.pdf"), format="pdf", bbox_inches="tight")
    plt.close(fig)

for direction in ["incoming", "outgoing"]:
    sub = partner_df[partner_df["direction"] == direction].copy()
    sub["compartment"] = pd.Categorical(sub["compartment"], categories=COMPARTMENT_ORDER, ordered=True)
    sub = sub.sort_values("compartment")
    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(sub)); w = 0.35
    ax.bar(x - w/2, sub["prop_SecA"].values, w, color=V2_PALETTE["SecA"],
           edgecolor="#333333", linewidth=0.3, label="SecA")
    ax.bar(x + w/2, sub["prop_SecB"].values, w, color=V2_PALETTE["SecB"],
           edgecolor="#333333", linewidth=0.3, label="SecB")
    ax.set_xticks(x); ax.set_xticklabels(sub["compartment"].values, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("% of pole's interactions", fontsize=8)
    ax.set_title(f"TME partner distribution - {direction} (NMF labels)", fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, frameon=False)
    plt.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, f"17c_tme_partners_{direction}.svg"), format="svg", bbox_inches="tight")
    plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 3))
groups = ["SecA", "SecB", "Intermediate", "Ciliated"]
x = np.arange(len(groups)); w = 0.35
inc_vals = [count_summary.get(g, {}).get("incoming", 0) for g in groups]
out_vals = [count_summary.get(g, {}).get("outgoing", 0) for g in groups]
ax.bar(x - w/2, inc_vals, w, color=[V2_PALETTE.get(g, "#ccc") for g in groups],
       edgecolor="#333333", linewidth=0.3, label="Incoming", alpha=0.7)
ax.bar(x + w/2, out_vals, w, color=[V2_PALETTE.get(g, "#ccc") for g in groups],
       edgecolor="#333333", linewidth=0.3, label="Outgoing", alpha=1.0)
ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
ax.set_ylabel("Significant interactions", fontsize=8)
ax.set_title("Interaction counts by NMF group", fontsize=10, fontweight="bold")
ax.legend(fontsize=7, frameon=False)
plt.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "17c_interaction_counts.svg"), format="svg", bbox_inches="tight")
plt.close(fig)

print("\nDone!")
