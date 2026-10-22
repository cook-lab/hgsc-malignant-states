#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/fig_secretory_polarization/00_prepare_data.py, the script that ran.
# Kept: section 1 and the save of meta.parquet (UMAP1/UMAP2 of the epithelial UMAP; the schema_nmf states read by Table S6
# and Fig 6h). Removed: the panel B-I data (site and treatment proportions, markers, radar, cell cycle, Wilcoxon DEG).
# The untrimmed script is in the lab project folder.
"""
Prepare all data for the Secretory Polarization manuscript figure.
Saves intermediate DataFrames as parquet/csv so the rendering script
can regenerate the PDF instantly without re-loading h5ad files.

Output → output/fig_secretory_polarization/data/
"""

import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")

# ============================================================================
# PATHS
# ============================================================================

BASE     = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD     = os.path.join(BASE, "celltype_h5ad", "hgsc_atlas_final_epithelial.h5ad")
DIR_11D  = os.path.join(BASE, "output", "11d_epithelial_nmf")
OUT      = os.path.join(BASE, "output", "fig_secretory_polarization", "data")
os.makedirs(OUT, exist_ok=True)

# ============================================================================
# 1. LOAD & BUILD META
# ============================================================================

print("[1] Loading h5ad (backed)...")
adata = sc.read_h5ad(H5AD, backed="r")
obs = adata.obs.copy()

# UMAP
umap = pd.DataFrame(
    adata.obsm["X_umap_local"], index=obs.index, columns=["UMAP1", "UMAP2"],
)

# NMF Factor 2
nmf_usage = pd.read_csv(os.path.join(DIR_11D, "11d_nmf_usage.csv"), index_col=0)
f2 = nmf_usage["Factor_2"]
shared = obs.index.intersection(f2.index)
obs = obs.loc[shared]
umap = umap.loc[shared]
f2 = f2.loc[shared]
print(f"  {len(shared):,} cells")

# Build meta
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

# NMF schema labels
is_ciliated = meta["celltype_level2"] == "Ciliated epithelial cell"
non_cil_f2 = meta.loc[~is_ciliated, "Factor_2"]
p50 = np.percentile(non_cil_f2, 50)
p75 = np.percentile(non_cil_f2, 75)
print(f"  NMF thresholds: p50={p50:.4f}, p75={p75:.4f}")

labels = pd.Series("SecA", index=meta.index)
labels[meta["Factor_2"] >= p50] = "Transitioning"
labels[meta["Factor_2"] >= p75] = "SecB"
labels[is_ciliated] = "Ciliated"
meta["schema_nmf"] = labels

for g in ["SecA", "Transitioning", "SecB", "Ciliated"]:
    n = (labels == g).sum()
    print(f"  {g}: {n:,} ({n/len(meta)*100:.1f}%)")

# ============================================================================
# 10. SAVE META (for UMAP panel A)
# ============================================================================

print("\n[10] Saving meta...")
meta.to_parquet(os.path.join(OUT, "meta.parquet"))
print(f"  Saved meta: {meta.shape}")

print("\nDone! All data saved to:", OUT)
