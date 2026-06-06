#!/usr/bin/env python3
"""
fig_secretory_polarization / 00_prepare_data — build the shared schema_nmf cache
================================================================================
PURPOSE
    Prepare every intermediate DataFrame for the secretory-polarization figure
    family and write them as parquet/CSV so the rendering scripts can rebuild
    panels instantly without re-loading the epithelial h5ad. This is the
    upstream producer of the shared `meta.parquet` (schema_nmf 4-class) that
    feeds Fig 1F/G/H/I, Fig 2C/D/E/G, and SF5/6/9.

INPUTS
    - obj("atlas_epithelial")  (epithelial h5ad: X_umap_local, obs, X for DE)
        NOTE: original sourced celltype_h5ad/hgsc_atlas_final_epithelial.h5ad;
        here that object is config key "atlas_epithelial".
    - 11d NMF usage  : output_root/03_epithelial_nmf/11d_nmf_usage.csv (Factor_2)
    - 11v2 score parquets (progeny/hallmark/dorothea/cell_cycle) and flux parquet,
        under output_root/04_functional/ (functional-characterization caches).

OUTPUTS (under fig_data_dir/)
    - meta.parquet                       (schema_nmf 4-class per epithelial cell)
    - panel_b_site_proportions.csv       (Fig 2D)
    - panel_c_paired_site.csv            (Fig 2E)
    - panel_d_treatment_proportions.csv
    - panel_e_paired_treatment.csv
    - panel_f_marker_expression.parquet
    - panel_g_{progeny,hallmark,dorothea,flux}_radar.csv
    - panel_h_phase_proportions.csv
    - panel_i_deg_results.csv            (Fig 2C volcano DE)

MANUSCRIPT PANEL(S): upstream cache for Fig 1F/G/H/I, Fig 2C/D/E/G, SF5/6/9.

RUNTIME TIER: heavy (loads epithelial h5ad twice; runs Wilcoxon DE).

NOTE: epithelial polarization label standardized to "Intermediate" (was
"Transitioning"). schema_nmf is written with the standardized labels, so all
downstream figure scripts select on "Intermediate".
"""
from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import CFG, obj, path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

H5AD = obj("atlas_epithelial")
USAGE_CSV = path("output_root", "03_epithelial_nmf", "11d_nmf_usage.csv")
FUNC_DIR = path("output_root", "04_functional")  # progeny/hallmark/dorothea/cell_cycle/flux parquets
OUT = path("output_root", "fig_secretory_polarization", "data")
os.makedirs(OUT, exist_ok=True)

# Canonical polarization labels (config = single source of truth).
LBL = CFG["polarization"]["labels"]          # [SecA, Intermediate, SecB, Ciliated]
ORDER = ["SecA", "Intermediate", "SecB", "Ciliated"]

MIN_CELLS = 100  # minimum epithelial cells per patient-site or patient-treatment

# ============================================================================
# 1. LOAD & BUILD META
# ============================================================================

print("[1] Loading epithelial h5ad (backed)...")
adata = sc.read_h5ad(H5AD, backed="r")
obs = adata.obs.copy()

umap = pd.DataFrame(
    adata.obsm["X_umap_local"], index=obs.index, columns=["UMAP1", "UMAP2"],
)

nmf_usage = pd.read_csv(USAGE_CSV, index_col=0)
f2 = nmf_usage[CFG["polarization"]["factor"]]   # Factor_2
shared = obs.index.intersection(f2.index)
obs = obs.loc[shared]
umap = umap.loc[shared]
f2 = f2.loc[shared]
print(f"  {len(shared):,} cells")

meta = pd.DataFrame(index=shared)
meta["UMAP1"] = umap["UMAP1"].values
meta["UMAP2"] = umap["UMAP2"].values
meta["celltype_level2"] = obs["celltype_level2"].values
meta["treatment_status"] = obs["treatment_status"].values
meta["anatomic_site"] = obs["anatomic_site"].values
meta["metastatic_site"] = obs["metastatic_site"].values
meta["treatment_response"] = obs["treatment_response"].values
meta["patient_id"] = obs["patient_id"].values
meta["Factor_2"] = f2.values

# NMF percentile partition (canonical schema; see config polarization.partition).
is_ciliated = meta["celltype_level2"] == "Ciliated epithelial cell"
non_cil_f2 = meta.loc[~is_ciliated, "Factor_2"]
p50 = np.percentile(non_cil_f2, 50)
p75 = np.percentile(non_cil_f2, 75)
print(f"  NMF thresholds: p50={p50:.4f}, p75={p75:.4f}")

labels = pd.Series("SecA", index=meta.index)
labels[meta["Factor_2"] >= p50] = "Intermediate"
labels[meta["Factor_2"] >= p75] = "SecB"
labels[is_ciliated] = "Ciliated"
meta["schema_nmf"] = labels

for g in ORDER:
    n = (labels == g).sum()
    print(f"  {g}: {n:,} ({n / len(meta) * 100:.1f}%)")

# ============================================================================
# 2. PANEL B — Pre-treatment proportions by metastatic site
# ============================================================================

print("\n[2] Panel B: pre-treatment proportions by metastatic site...")
pre = meta[meta["treatment_status"] == "pre-treatment"].copy()
pre = pre[pre["metastatic_site"].isin(["primary", "ascites", "metastasis"])]

site_counts = pre.groupby(["metastatic_site", "schema_nmf"]).size().unstack(fill_value=0)
site_props = site_counts.div(site_counts.sum(axis=1), axis=0) * 100
site_props.to_csv(os.path.join(OUT, "panel_b_site_proportions.csv"))
print(site_props.round(1))

# ============================================================================
# 3. PANEL C — Paired patient SecB proportion across metastatic sites
# ============================================================================

print("\n[3] Panel C: paired patient SecB across metastatic sites...")
pre_sites = pre.copy()
pat_site = pre_sites.groupby(["patient_id", "metastatic_site"]).apply(
    lambda df: pd.Series({
        "n_cells": len(df),
        "secb_prop": (df["schema_nmf"] == "SecB").mean(),
    })
).reset_index()
pat_site = pat_site[pat_site["n_cells"] >= MIN_CELLS]

pat_counts = pat_site.groupby("patient_id")["metastatic_site"].nunique()
multi_site_pats = pat_counts[pat_counts >= 2].index
paired_c = pat_site[pat_site["patient_id"].isin(multi_site_pats)].copy()
print(f"  {len(multi_site_pats)} patients with >= 2 sites (>= {MIN_CELLS} cells each)")
print(f"  {len(paired_c)} data points")
paired_c.to_csv(os.path.join(OUT, "panel_c_paired_site.csv"), index=False)

# ============================================================================
# 4. PANEL D — Adnexa pre vs post-chemo proportions
# ============================================================================

print("\n[4] Panel D: adnexa pre vs post-chemo proportions...")
adnexa = meta[meta["anatomic_site"] == "adnexa"].copy()
adnexa_tx = adnexa[adnexa["treatment_status"].isin(["pre-treatment", "post-chemotherapy"])]

tx_counts = adnexa_tx.groupby(["treatment_status", "schema_nmf"]).size().unstack(fill_value=0)
tx_props = tx_counts.div(tx_counts.sum(axis=1), axis=0) * 100
tx_props.to_csv(os.path.join(OUT, "panel_d_treatment_proportions.csv"))
print(tx_props.round(1))

# ============================================================================
# 5. PANEL E — Paired patient SecB pre vs post-chemo
# ============================================================================

print("\n[5] Panel E: paired patient SecB pre vs post-chemo...")
all_tx = meta[meta["treatment_status"].isin(["pre-treatment", "post-chemotherapy"])].copy()
pat_tx = all_tx.groupby(["patient_id", "treatment_status"]).apply(
    lambda df: pd.Series({
        "n_cells": len(df),
        "secb_prop": (df["schema_nmf"] == "SecB").mean(),
    })
).reset_index()
pat_tx = pat_tx[pat_tx["n_cells"] >= MIN_CELLS]

pat_both = pat_tx.groupby("patient_id")["treatment_status"].nunique()
paired_pats = pat_both[pat_both == 2].index
paired_e = pat_tx[pat_tx["patient_id"].isin(paired_pats)].copy()
print(f"  {len(paired_pats)} patients with both pre + post (>= {MIN_CELLS} cells each)")
paired_e.to_csv(os.path.join(OUT, "panel_e_paired_treatment.csv"), index=False)

# ============================================================================
# 6. PANEL F — Marker gene expression
# ============================================================================

print("\n[6] Panel F: extracting marker gene expression...")
MARKERS_F = ["KRT7", "KRT19", "TACSTD2", "BCAM", "MECOM", "SOX17"]
marker_idx = [adata.var_names.get_loc(g) for g in MARKERS_F]
expr_raw = adata.X[:, marker_idx]
if sparse.issparse(expr_raw):
    expr_raw = expr_raw.toarray()
expr_df = pd.DataFrame(expr_raw, index=adata.obs.index, columns=MARKERS_F)
expr_df = expr_df.loc[shared]
expr_df["UMAP1"] = meta["UMAP1"].values
expr_df["UMAP2"] = meta["UMAP2"].values
expr_df.to_parquet(os.path.join(OUT, "panel_f_marker_expression.parquet"))
print(f"  Saved {len(MARKERS_F)} genes for {len(expr_df):,} cells")

# ============================================================================
# 7. PANEL G — Radar data (select most differential features)
# ============================================================================

print("\n[7] Panel G: computing radar data (most differential features)...")


def compute_group_means(scores_df, group_labels, feature_cols):
    numeric_cols = [c for c in feature_cols
                    if scores_df[c].dtype in ["float64", "float32", "int64", "int32"]]
    df = scores_df[numeric_cols].copy()
    df["group"] = group_labels
    return df.groupby("group").mean(numeric_only=True).reindex(ORDER)


def select_top_differential(group_means, n=8):
    ranges = group_means.max(axis=0) - group_means.min(axis=0)
    return ranges.nlargest(n).index.tolist()


def zscore_df(df):
    mu = df.mean(axis=0)
    sd = df.std(axis=0).replace(0, 1)
    return (df - mu) / sd


group_labels = meta["schema_nmf"].values

progeny = pd.read_parquet(os.path.join(FUNC_DIR, "progeny_scores.parquet"))
progeny = progeny.loc[shared].drop(columns=["celltype_v2"], errors="ignore")
prog_means = compute_group_means(progeny, group_labels, progeny.columns.tolist())
prog_top = select_top_differential(prog_means, n=8)
zscore_df(prog_means[prog_top]).to_csv(os.path.join(OUT, "panel_g_progeny_radar.csv"))
print(f"  PROGENy top 8: {prog_top}")

hallmark = pd.read_parquet(os.path.join(FUNC_DIR, "hallmark_scores.parquet"))
hallmark = hallmark.loc[shared].drop(columns=["celltype_v2"], errors="ignore")
hall_means = compute_group_means(hallmark, group_labels, hallmark.columns.tolist())
hall_top = select_top_differential(hall_means, n=8)
zscore_df(hall_means[hall_top]).to_csv(os.path.join(OUT, "panel_g_hallmark_radar.csv"))
print(f"  Hallmark top 8: {hall_top}")

dorothea = pd.read_parquet(os.path.join(FUNC_DIR, "dorothea_scores.parquet"))
dorothea = dorothea.loc[shared].drop(columns=["celltype_v2", "__index_level_0__"], errors="ignore")
doro_means = compute_group_means(dorothea, group_labels, dorothea.columns.tolist())
doro_top = select_top_differential(doro_means, n=8)
zscore_df(doro_means[doro_top]).to_csv(os.path.join(OUT, "panel_g_dorothea_radar.csv"))
print(f"  DoRoThEA top 8: {doro_top}")

flux = pd.read_parquet(os.path.join(FUNC_DIR, "flux_per_cell.parquet"))
flux_shared = flux.index.intersection(shared)
flux = flux.loc[flux_shared]
flux_labels = meta.loc[flux_shared, "schema_nmf"].values
flux_means = compute_group_means(flux, flux_labels, flux.columns.tolist())
flux_top = select_top_differential(flux_means, n=8)
zscore_df(flux_means[flux_top]).to_csv(os.path.join(OUT, "panel_g_flux_radar.csv"))
print(f"  Flux top 8: {flux_top}")

# ============================================================================
# 8. PANEL H — G2M cell cycle data
# ============================================================================

print("\n[8] Panel H: cell cycle phase proportions...")
cell_cycle = pd.read_parquet(os.path.join(FUNC_DIR, "cell_cycle_scores.parquet"))
cell_cycle = cell_cycle.loc[shared]
phase_ct = pd.crosstab(
    pd.Series(group_labels, name="group"),
    cell_cycle["phase"],
    normalize="index",
) * 100
phase_ct = phase_ct.reindex(ORDER).reindex(columns=["G1", "S", "G2M"], fill_value=0)
phase_ct.to_csv(os.path.join(OUT, "panel_h_phase_proportions.csv"))
print(phase_ct.round(1))

# ============================================================================
# 9. PANEL I — DEG volcano (SecA vs SecB)
# ============================================================================

print("\n[9] Panel I: running DEG SecA vs SecB...")
seca_idx = meta.index[meta["schema_nmf"] == "SecA"]
secb_idx = meta.index[meta["schema_nmf"] == "SecB"]
print(f"  SecA: {len(seca_idx):,}, SecB: {len(secb_idx):,}")
print("  Loading expression data (this may take a while)...")

np.random.seed(SEED)
rng = np.random.default_rng(SEED)
MAX_CELLS = 50000
seca_sub = (rng.choice(seca_idx, MAX_CELLS, replace=False)
            if len(seca_idx) > MAX_CELLS else seca_idx.values)
secb_sub = (rng.choice(secb_idx, MAX_CELLS, replace=False)
            if len(secb_idx) > MAX_CELLS else secb_idx.values)

all_sub = np.concatenate([seca_sub, secb_sub])
print(f"  Subsampled: {len(seca_sub):,} SecA + {len(secb_sub):,} SecB = {len(all_sub):,}")

adata_deg = sc.read_h5ad(H5AD)
adata_deg = adata_deg[all_sub].copy()
adata_deg.obs["deg_group"] = "SecA"
adata_deg.obs.loc[secb_sub, "deg_group"] = "SecB"

print("  Running Wilcoxon rank-sum...")
sc.tl.rank_genes_groups(
    adata_deg, groupby="deg_group", groups=["SecB"], reference="SecA",
    method="wilcoxon", n_genes=adata_deg.n_vars, use_raw=False,
)

result = adata_deg.uns["rank_genes_groups"]
deg_df = pd.DataFrame({
    "gene": result["names"]["SecB"],
    "log2fc": result["logfoldchanges"]["SecB"],
    "score": result["scores"]["SecB"],
    "pval": result["pvals"]["SecB"],
    "pval_adj": result["pvals_adj"]["SecB"],
})
deg_df.to_csv(os.path.join(OUT, "panel_i_deg_results.csv"), index=False)
print(f"  Saved {len(deg_df):,} genes")
sig = deg_df[(deg_df["pval_adj"] < 0.05) & (deg_df["log2fc"].abs() > 0.25)]
print(f"  Significant: {len(sig):,} (|log2FC|>0.25, padj<0.05)")

# ============================================================================
# 10. SAVE META
# ============================================================================

print("\n[10] Saving meta...")
meta.to_parquet(os.path.join(OUT, "meta.parquet"))
print(f"  Saved meta: {meta.shape}")

del adata, adata_deg
print("\nDone! All data saved to:", OUT)
