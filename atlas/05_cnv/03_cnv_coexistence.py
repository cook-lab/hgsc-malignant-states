#!/usr/bin/env python3
"""
CNV — within-clone SecA/SecB coexistence analysis
=================================================
HGSC malignant-states atlas backend.

The key manuscript analysis: shows that SecA and SecB cells coexist within the
same CopyKAT CNV subclones, falsifying a purely clonal explanation for SecB.
Per clone (>=20 epithelial cells) computes SecA/SecB/Intermediate proportions,
the Factor_2 distribution, and epitype entropy; assigns a per-sample verdict
(monoclonal / clonally_driven / mixed).

INPUTS (output_root/05_cnv/):
  - per_sample/*/copykat_subclones.csv, per_sample/*/barcodes.csv
  - output_root/03_epithelial_nmf/11d_nmf_usage.csv (Factor_2)
  - <data_root>/2026_final_atlas/celltype_h5ad/hgsc_atlas_final_epithelial.h5ad (UMAP)

OUTPUTS (output_root/05_cnv/):
  - tables/within_clone_coexistence.csv   (KEY cache — Fig 1J, SF7)
  - tables/per_sample_verdict.csv
  - figs/19_cnv_*.svg/pdf

MANUSCRIPT PANELS: Fig 1J (alluvial cache), SF7 (within-clone bars).

RUNTIME TIER: moderate (per-sample CopyKAT outputs + UMAP read).

SEEDING: deterministic (plot subsamples per-subtype but cap-based; no RNG draws
  affecting analytical results).

Usage:
    python 03_cnv_coexistence.py
"""

import os
import sys
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import entropy as scipy_entropy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path  # noqa: E402

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

# ── Paths ────────────────────────────────────────────────────
CNV_DIR    = path("output_root", "05_cnv")
PER_SAMPLE = os.path.join(CNV_DIR, "per_sample")
TBL_DIR    = os.path.join(CNV_DIR, "tables")
FIG_DIR    = os.path.join(CNV_DIR, "figs")
NMF_CSV    = path("output_root", "03_epithelial_nmf", "11d_nmf_usage.csv")
H5AD_EPI   = path("data_root", "2026_final_atlas", "celltype_h5ad",
                  "hgsc_atlas_final_epithelial.h5ad")
os.makedirs(TBL_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

EPITYPE_PALETTE = {
    "SecA": "#E6A141", "Intermediate": "#CF8C2E",
    "SecB": "#B8741A", "Ciliated": "#E05A2C",
}

MIN_CELLS_PER_CLONE = 20
# Verdict thresholds (see per-sample verdict logic below)
DOMINANCE_THRESH = 0.70


def compute_sample_verdict(clone_rows_for_sample):
    df = pd.DataFrame(clone_rows_for_sample)
    n_clones = len(df)
    if n_clones == 0:
        return None
    if n_clones == 1:
        return {
            "n_clones": 1,
            "dom_clone_secA_frac": 1.0,
            "dom_clone_secB_frac": 1.0,
            "within_clone_entropy_mean": float(df["entropy"].mean()),
            "verdict": "monoclonal",
        }
    total_secA = df["n_secA"].sum()
    total_secB = df["n_secB"].sum()
    dom_secA = df["n_secA"].max() / total_secA if total_secA > 0 else 0.0
    dom_secB = df["n_secB"].max() / total_secB if total_secB > 0 else 0.0
    verdict = "clonally_driven" if (dom_secA >= DOMINANCE_THRESH and
                                    dom_secB >= DOMINANCE_THRESH) else "mixed"
    return {
        "n_clones": n_clones,
        "dom_clone_secA_frac": dom_secA,
        "dom_clone_secB_frac": dom_secB,
        "within_clone_entropy_mean": float(df["entropy"].mean()),
        "verdict": verdict,
    }


def main():
    print("=" * 65)
    print("  CNV — Within-Clone SecA/SecB Coexistence")
    print("=" * 65)

    print("\n[1] Loading NMF Factor 2 scores...")
    f2 = pd.read_csv(NMF_CSV, index_col=0)["Factor_2"]
    print(f"  {len(f2):,} cells with Factor 2 scores")

    print("\n[2] Loading UMAP coordinates...")
    adata = ad.read_h5ad(H5AD_EPI, backed="r")
    umap_key = "X_umap_local" if "X_umap_local" in adata.obsm else "X_umap"
    umap_df = pd.DataFrame(adata.obsm[umap_key], index=adata.obs.index,
                           columns=["UMAP1", "UMAP2"])
    del adata
    print(f"  {len(umap_df):,} cells with UMAP coords")

    print("\n[3] Analyzing within-clone coexistence...")
    completed = [s for s in sorted(os.listdir(PER_SAMPLE))
                 if os.path.isdir(os.path.join(PER_SAMPLE, s))
                 and os.path.exists(os.path.join(PER_SAMPLE, s, "DONE.txt"))]
    print(f"  Completed samples: {len(completed)}")

    clone_rows, sample_rows, example_samples = [], [], []
    for sample_id in completed:
        sample_dir = os.path.join(PER_SAMPLE, sample_id)
        subclone_path = os.path.join(sample_dir, "copykat_subclones.csv")
        if not os.path.exists(subclone_path):
            continue
        subclones = pd.read_csv(subclone_path)
        barcodes = pd.read_csv(os.path.join(sample_dir, "barcodes.csv"))

        epi = barcodes[~barcodes["is_reference"]].copy()
        merged = epi.merge(subclones, on="barcode", how="inner").set_index("barcode")
        merged["factor2"] = f2.reindex(merged.index)
        merged = merged.dropna(subset=["factor2"])
        merged_aneuploid = merged[merged["subclone"] != "diploid"]
        if len(merged_aneuploid) < MIN_CELLS_PER_CLONE:
            continue

        n_secA = (merged_aneuploid["epitype"] == "SecA").sum()
        n_secB = (merged_aneuploid["epitype"] == "SecB").sum()
        n_interm = (merged_aneuploid["epitype"] == "Intermediate").sum()
        total = n_secA + n_secB + n_interm
        sample_secB_frac = n_secB / total if total > 0 else 0

        n_informative_clones = 0
        for clone_id in merged_aneuploid["subclone"].unique():
            clone_cells = merged_aneuploid[merged_aneuploid["subclone"] == clone_id]
            if len(clone_cells) < MIN_CELLS_PER_CLONE:
                continue
            n_informative_clones += 1
            cs = clone_cells["epitype"].value_counts()
            clone_secA = cs.get("SecA", 0)
            clone_secB = cs.get("SecB", 0)
            clone_interm = cs.get("Intermediate", 0)
            clone_total = clone_secA + clone_secB + clone_interm

            frac_secA = clone_secA / clone_total if clone_total > 0 else 0
            frac_secB = clone_secB / clone_total if clone_total > 0 else 0
            frac_interm = clone_interm / clone_total if clone_total > 0 else 0

            probs = np.array([frac_secA, frac_interm, frac_secB])
            probs = probs[probs > 0]
            ent = scipy_entropy(probs, base=2) if len(probs) > 1 else 0

            f2_vals = clone_cells["factor2"].values
            clone_rows.append({
                "sample_id": sample_id, "subclone": clone_id,
                "n_cells": len(clone_cells), "n_secA": clone_secA,
                "n_intermediate": clone_interm, "n_secB": clone_secB,
                "frac_secA": frac_secA, "frac_intermediate": frac_interm,
                "frac_secB": frac_secB, "entropy": ent,
                "factor2_mean": f2_vals.mean(), "factor2_std": f2_vals.std(),
                "has_both_secA_secB": (clone_secA > 0) and (clone_secB > 0),
                "sample_secB_frac": sample_secB_frac,
            })

        sample_rows.append({
            "sample_id": sample_id, "n_clones": n_informative_clones,
            "n_cells": len(merged_aneuploid), "sample_secB_frac": sample_secB_frac,
        })
        if n_informative_clones >= 2 and len(example_samples) < 6:
            example_samples.append(sample_id)

    clone_df = pd.DataFrame(clone_rows)
    sample_df = pd.DataFrame(sample_rows)
    print(f"  Samples with >=2 informative clones: {len(sample_df)}")
    print(f"  Total informative clones: {len(clone_df)}")

    if len(clone_df) > 0:
        n_with_both = clone_df["has_both_secA_secB"].sum()
        print(f"  Clones with BOTH SecA and SecB: {n_with_both}/{len(clone_df)} "
              f"({100*n_with_both/len(clone_df):.1f}%)")
        print(f"  Median within-clone entropy: {clone_df['entropy'].median():.3f}")

    clone_df.to_csv(os.path.join(TBL_DIR, "within_clone_coexistence.csv"), index=False)
    print("\n  Saved within_clone_coexistence.csv")

    verdict_rows = []
    if len(clone_df) > 0:
        for sid, grp in clone_df.groupby("sample_id"):
            v = compute_sample_verdict(grp.to_dict("records"))
            if v is None:
                continue
            v["sample_id"] = sid
            v["sample_secB_frac"] = grp["sample_secB_frac"].iloc[0]
            verdict_rows.append(v)
        verdict_df = pd.DataFrame(verdict_rows)[
            ["sample_id", "n_clones", "sample_secB_frac",
             "dom_clone_secA_frac", "dom_clone_secB_frac",
             "within_clone_entropy_mean", "verdict"]
        ]
        verdict_df.to_csv(os.path.join(TBL_DIR, "per_sample_verdict.csv"), index=False)
        print("  Saved per_sample_verdict.csv")
        for k, v in verdict_df["verdict"].value_counts().items():
            print(f"    {k}: {v}")

    print("\n[4] Generating figures...")
    if len(clone_df) == 0:
        print("  No informative clones — skipping figures.")
        return

    # ── Panel A: Example UMAP (subclone vs epitype) ──────────
    if example_samples:
        n_examples = min(4, len(example_samples))
        fig, axes = plt.subplots(n_examples, 2, figsize=(8, 3.5 * n_examples))
        if n_examples == 1:
            axes = axes.reshape(1, 2)
        clone_colors = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
        for row_i, sid in enumerate(example_samples[:n_examples]):
            sample_dir = os.path.join(PER_SAMPLE, sid)
            subclones = pd.read_csv(os.path.join(sample_dir, "copykat_subclones.csv"))
            barcodes = pd.read_csv(os.path.join(sample_dir, "barcodes.csv"))
            epi = barcodes[~barcodes["is_reference"]]
            merged = epi.merge(subclones, on="barcode", how="inner")
            merged = merged[merged["subclone"] != "diploid"].set_index("barcode")
            common = merged.index.intersection(umap_df.index)
            if len(common) < 50:
                continue
            merged = merged.loc[common]
            u1 = umap_df.loc[common, "UMAP1"].values
            u2 = umap_df.loc[common, "UMAP2"].values

            ax = axes[row_i, 0]
            for ci, clone_id in enumerate(sorted(merged["subclone"].unique())):
                mask = merged["subclone"].values == clone_id
                ax.scatter(u1[mask], u2[mask], s=1,
                           c=clone_colors[ci % len(clone_colors)], alpha=0.6,
                           rasterized=True, label=clone_id)
            ax.set_title(f"{sid} — CNV Subclone", fontsize=7, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            ax.legend(fontsize=5, markerscale=5, frameon=False)

            ax = axes[row_i, 1]
            for ep in ["SecA", "Intermediate", "SecB"]:
                mask = merged["epitype"].values == ep
                if mask.any():
                    ax.scatter(u1[mask], u2[mask], s=1,
                               c=EPITYPE_PALETTE.get(ep, "#999"),
                               alpha=0.6, rasterized=True, label=ep)
            ax.set_title(f"{sid} — Epitype", fontsize=7, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            ax.legend(fontsize=5, markerscale=5, frameon=False)

        fig.suptitle("CNV Subclone vs SecA/SecB Epitype — Example Samples",
                     fontsize=10, fontweight="bold", y=1.01)
        plt.tight_layout()
        for fmt in ["svg", "pdf"]:
            fig.savefig(os.path.join(FIG_DIR, f"19_cnv_example_umaps.{fmt}"),
                        bbox_inches="tight")
        plt.close()
        print("  Saved example_umaps")

    # ── Panel B: Within-clone stacked bars ────────────────────
    if len(sample_df) > 0:
        top_samples = sample_df.nlargest(10, "n_clones")["sample_id"].tolist()
        plot_clones = clone_df[clone_df["sample_id"].isin(top_samples)].copy()
        plot_clones["label"] = plot_clones["sample_id"] + "\n" + plot_clones["subclone"]
        fig, ax = plt.subplots(figsize=(max(8, len(plot_clones) * 0.4), 4))
        bar_x = np.arange(len(plot_clones))
        bottom = np.zeros(len(plot_clones))
        for ep, color in [("frac_secA", "#E6A141"),
                          ("frac_intermediate", "#CF8C2E"),
                          ("frac_secB", "#B8741A")]:
            vals = plot_clones[ep].values * 100
            ax.bar(bar_x, vals, bottom=bottom, color=color, width=0.7,
                   edgecolor="#333333", linewidth=0.3)
            bottom += vals
        ax.set_xticks(bar_x)
        ax.set_xticklabels(plot_clones["label"].values, fontsize=5, rotation=45, ha="right")
        ax.set_ylabel("Proportion (%)"); ax.set_ylim(0, 100)
        ax.set_title("Within-Clone Epitype Proportions", fontweight="bold")
        ax.legend(handles=[Patch(facecolor="#E6A141", label="SecA"),
                           Patch(facecolor="#CF8C2E", label="Intermediate"),
                           Patch(facecolor="#B8741A", label="SecB")],
                  fontsize=6, frameon=False)
        plt.tight_layout()
        for fmt in ["svg", "pdf"]:
            fig.savefig(os.path.join(FIG_DIR, f"19_cnv_within_clone_barplot.{fmt}"),
                        bbox_inches="tight")
        plt.close()
        print("  Saved within_clone_barplot")

    # ── Panel C: Within-clone vs sample-level SecB scatter ────
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(clone_df["sample_secB_frac"] * 100, clone_df["frac_secB"] * 100,
               s=8, c="#7D4E4E", alpha=0.5, edgecolors="none")
    ax.plot([0, 100], [0, 100], color="red", ls="--", lw=1, label="Perfect agreement")
    ax.set_xlabel("Sample-level SecB fraction (%)")
    ax.set_ylabel("Within-clone SecB fraction (%)")
    ax.set_title("Clone-Level vs Sample-Level SecB Proportion", fontweight="bold")
    ax.legend(fontsize=6)
    ax.set_xlim(0, max(50, clone_df["sample_secB_frac"].max() * 110))
    ax.set_ylim(0, max(50, clone_df["frac_secB"].max() * 110))
    r = np.corrcoef(clone_df["sample_secB_frac"], clone_df["frac_secB"])[0, 1]
    ax.text(0.05, 0.95, f"r = {r:.2f}", transform=ax.transAxes, fontsize=7, va="top")
    plt.tight_layout()
    for fmt in ["svg", "pdf"]:
        fig.savefig(os.path.join(FIG_DIR, f"19_cnv_clone_vs_sample_secB.{fmt}"),
                    bbox_inches="tight")
    plt.close()
    print("  Saved clone_vs_sample_secB scatter")

    # ── Panel D: Entropy distribution ─────────────────────────
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.hist(clone_df["entropy"], bins=20, color="#E6A141", alpha=0.8,
            edgecolor="#333333", linewidth=0.4)
    ax.set_xlabel("Within-clone epitype entropy (bits)")
    ax.set_ylabel("Number of clones")
    ax.set_title("Epitype Diversity Within CNV Subclones", fontweight="bold")
    med = clone_df["entropy"].median()
    ax.axvline(med, color="red", ls="--", lw=1, label=f"Median = {med:.2f}")
    ax.legend(fontsize=6)
    plt.tight_layout()
    for fmt in ["svg", "pdf"]:
        fig.savefig(os.path.join(FIG_DIR, f"19_cnv_entropy.{fmt}"), bbox_inches="tight")
    plt.close()
    print("  Saved entropy histogram")

    # ── Per-sample grid (one panel per testable sample) ──────
    testable_ids = [r["sample_id"] for r in verdict_rows if r["verdict"] != "monoclonal"]
    if testable_ids:
        n = len(testable_ids)
        ncols = 8
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 1.6, nrows * 1.6))
        axes = np.atleast_2d(axes).ravel()
        verdict_lookup = {r["sample_id"]: r["verdict"] for r in verdict_rows}
        for i, sid in enumerate(testable_ids):
            ax = axes[i]
            g = clone_df[clone_df["sample_id"] == sid].sort_values("subclone")
            x = np.arange(len(g))
            bottom = np.zeros(len(g))
            for ep, color in [("frac_secA", "#E6A141"),
                              ("frac_intermediate", "#CF8C2E"),
                              ("frac_secB", "#B8741A")]:
                vals = g[ep].values * 100
                ax.bar(x, vals, bottom=bottom, color=color, width=0.85, edgecolor="none")
                bottom += vals
            border = "#CC3333" if verdict_lookup.get(sid) == "clonally_driven" else "#888888"
            for spine in ax.spines.values():
                spine.set_edgecolor(border)
                spine.set_linewidth(1.0 if verdict_lookup.get(sid) == "clonally_driven" else 0.4)
            ax.set_xticks([]); ax.set_yticks([]); ax.set_ylim(0, 100)
            ax.set_title(f"{sid}\nk={len(g)}", fontsize=5, color=border, pad=1)
        for j in range(len(testable_ids), len(axes)):
            axes[j].axis("off")
        fig.legend(handles=[Patch(facecolor="#E6A141", label="SecA"),
                            Patch(facecolor="#CF8C2E", label="Intermediate"),
                            Patch(facecolor="#B8741A", label="SecB"),
                            Patch(facecolor="none", edgecolor="#CC3333",
                                  linewidth=1.0, label="clonally_driven")],
                   fontsize=6, frameon=False, loc="lower center", ncol=4,
                   bbox_to_anchor=(0.5, -0.02))
        fig.suptitle("Per-Sample: Epitype Composition Within Each CNV Subclone",
                     fontsize=9, fontweight="bold", y=1.00)
        plt.tight_layout()
        for fmt in ["svg", "pdf"]:
            fig.savefig(os.path.join(FIG_DIR, f"19_cnv_per_sample_grid.{fmt}"),
                        bbox_inches="tight")
        plt.close()
        print(f"  Saved per_sample_grid ({n} samples)")

    print(f"\n{'=' * 65}")
    print("  DONE — Within-clone coexistence analysis")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
