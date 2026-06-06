#!/usr/bin/env python3
"""
Atlas 01 — Step 02: concatenate per-study matrices, QC filter, Scrublet doublets

PURPOSE
    Concatenate the harmonised per-study count matrices (from step 01) on their
    gene intersection, enforce a common obs schema, apply per-cell QC cutoffs,
    drop low-cell samples, run Scrublet per sample, and filter doublets
    (score < 0.3). Produces the concatenated counts object that is the INPUT to
    the scVI/scANVI integration trust boundary (see this stage's README).

INPUTS
    DATA_ROOT/2026_final_atlas/processed/<study>.h5ad   (13 per-study objects from step 01)

OUTPUTS
    DATA_ROOT/2026_final_atlas/processed/atlas_concat_counts_only_X.h5ad   (concatenated counts)
    DATA_ROOT/2026_final_atlas/processed/atlas_concatenated_filtered.h5ad   (QC + Scrublet<0.3)
    output_root/01_preprocess_qc/02_qc/*.svg  (QC violins)

MANUSCRIPT PANEL(S)
    Pre-integration provenance; no panel rendered directly. The integration that
    consumes atlas_concatenated_filtered.h5ad is NOT reproduced here.

RUNTIME TIER
    heavy (concatenation of ~13 study matrices; per-sample Scrublet).

NOTE
    Migrated from 02_concat_qc_doublets.ipynb. Logic preserved exactly; the
    duplicated diagnostic cells and interactive plt.show() panels were removed,
    and QC violins are written to file. Scrublet seed pinned to config SEED.
"""

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# PATHS (resolved via central config)
# ============================================================================

PROCESSED_DIR = path("data_root", "2026_final_atlas", "processed")
FIG_DIR       = path("output_root", "01_preprocess_qc", "02_qc")
os.makedirs(FIG_DIR, exist_ok=True)

DATASET_NAMES = [
    "denisenko_2024", "geistlinger_2020", "loret_2022", "luo_2024", "nath_2021",
    "olalekan_2021", "olbrecht_2021", "regner_2021", "vazquez_garcia_2022",
    "xu_2022", "zhang_2022", "zheng_2023", "hornburg_2021",
]
dataset_paths = [os.path.join(PROCESSED_DIR, f"{n}.h5ad") for n in DATASET_NAMES]

CONCAT_PATH   = os.path.join(PROCESSED_DIR, "atlas_concat_counts_only_X.h5ad")
FILTERED_PATH = os.path.join(PROCESSED_DIR, "atlas_concatenated_filtered.h5ad")

# ============================================================================
# LOAD DATASETS
# ============================================================================

adata_objects = []
dataset_names = []
for p, name in zip(dataset_paths, DATASET_NAMES):
    print(f"Loading: {name}")
    adata_objects.append(sc.read_h5ad(p))
    dataset_names.append(name)
print(f"\nLoaded {len(adata_objects)} datasets")

# ============================================================================
# ENFORCE OBS SCHEMA (string-only, categorical-safe)
# ============================================================================

required_obs = [
    "barcode", "sample_id", "study", "patient_id", "sample_num", "treatment_status",
    "histological_subtype", "stage", "anatomic_site", "metastatic_site", "age",
    "treatment_response", "BRCA_status", "HRD_status", "TP53_status", "ref",
]
for i, a in enumerate(adata_objects):
    for col in required_obs:
        if col not in a.obs.columns:
            a.obs[col] = "NA"
    a.obs = a.obs.loc[:, required_obs]
    for col in required_obs:
        s = a.obs[col].astype("string").fillna("NA").replace("", "NA")
        a.obs[col] = s.astype(str)
    adata_objects[i] = a
print(f"Obs schema enforced ({len(required_obs)} columns)")

# ============================================================================
# HARMONIZE GENES (INTERSECTION)
# ============================================================================

common_genes = set(adata_objects[0].var_names)
for a in adata_objects[1:]:
    common_genes &= set(a.var_names)
common_genes = np.array(sorted(common_genes))
print(f"Number of common genes across all datasets: {len(common_genes)}")

for i, a in enumerate(adata_objects):
    a_sub = a[:, common_genes].copy()
    a_sub.X = a_sub.X.tocsr() if sp.issparse(a_sub.X) else sp.csr_matrix(a_sub.X)
    if a_sub.X.dtype != np.int32:
        a_sub.X.data = np.round(a_sub.X.data).astype(np.int32)
    a_sub.raw = None
    if hasattr(a_sub, "layers"):
        a_sub.layers.clear()
    adata_objects[i] = a_sub

# ============================================================================
# CONCATENATE
# ============================================================================

print("Merging data ...", flush=True)
adata_concat = sc.concat(adata_objects, join="inner", label="dataset", keys=dataset_names)
if not sp.isspmatrix_csr(adata_concat.X):
    adata_concat.X = adata_concat.X.tocsr()
if adata_concat.X.dtype != np.int32:
    adata_concat.X.data = np.rint(adata_concat.X.data).astype(np.int32)
adata_concat.obs_names_make_unique()
if "counts" in adata_concat.layers:
    del adata_concat.layers["counts"]
adata_concat.raw = None

print("Writing concatenated counts h5ad ...", flush=True)
adata_concat.write_h5ad(CONCAT_PATH, compression="gzip")

# Reload (matches original notebook flow: full load before filtering)
adata_concat = sc.read_h5ad(CONCAT_PATH)
print(f"Number of cells: {adata_concat.n_obs}")

# ============================================================================
# SAMPLE-LEVEL FILTERS
# ============================================================================

adata_concat = adata_concat[
    adata_concat.obs["study"].notna() & (adata_concat.obs["study"] != ""), :
].copy()
print(f"Cells after removing missing study: {adata_concat.n_obs}")

adata_concat = adata_concat[adata_concat.obs["histological_subtype"] == "Serous"].copy()
print(f"Cells after filtering histological_subtype != 'Serous': {adata_concat.n_obs}")

# ============================================================================
# PER-CELL QC METRICS + CUTOFFS
# ============================================================================

adata_concat.var_names_make_unique()
sc.pp.calculate_qc_metrics(adata_concat, inplace=True, percent_top=None, log1p=False)

CUTS = {"total_counts": {"min": 500}, "n_genes_by_counts": {"min": 300}}

# QC violins (pre-filter) → file
keys = ["total_counts", "n_genes_by_counts"]
axs = sc.pl.violin(adata_concat, keys=keys, groupby="study", log=False, rotation=90,
                   stripplot=False, multi_panel=True, show=False, scale="width")
if not isinstance(axs, (list, tuple)):
    axs = [axs]
for ax, key in zip(axs, keys):
    if "min" in CUTS[key]:
        ax.axhline(CUTS[key]["min"], linestyle="--", linewidth=1)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "02_qc_violins_prefilter.svg"), bbox_inches="tight")
plt.close("all")

qc_mask = (
    (adata_concat.obs["total_counts"] >= CUTS["total_counts"]["min"]) &
    (adata_concat.obs["n_genes_by_counts"] >= CUTS["n_genes_by_counts"]["min"])
)
adata_concat = adata_concat[qc_mask].copy()
print(f"Cells after QC filtering: {adata_concat.n_obs}")

# Drop samples with < 500 cells
sample_counts = adata_concat.obs["sample_num"].value_counts()
valid = sample_counts[sample_counts >= 500].index
adata_concat = adata_concat[adata_concat.obs["sample_num"].isin(valid)].copy()
print(f"Cells after filtering sample_num < 500: {adata_concat.n_obs}")

# ============================================================================
# SCRUBLET (per sample)
# ============================================================================

batch_key = "sample_id"
samples = adata_concat.obs[batch_key].unique()
doublet_scores = pd.Series(index=adata_concat.obs_names, dtype=float)
predicted_doublets = pd.Series(index=adata_concat.obs_names, dtype=bool)

for i, sample in enumerate(samples):
    mask = adata_concat.obs[batch_key] == sample
    adata_sub = adata_concat[mask].copy()
    if adata_sub.n_obs < 1000:
        doublet_scores.loc[mask] = np.nan
        predicted_doublets.loc[mask] = False
        print(f"[{i+1}/{len(samples)}] {sample}: SKIPPED — only {adata_sub.n_obs} cells (min 1000)")
        continue
    expected_rate = 0.008 * (adata_sub.n_obs / 1000)
    sc.external.pp.scrublet(adata_sub, expected_doublet_rate=expected_rate, random_state=SEED)
    doublet_scores.loc[mask] = adata_sub.obs["doublet_score"].values
    predicted_doublets.loc[mask] = adata_sub.obs["predicted_doublet"].values
    n_doublets = adata_sub.obs["predicted_doublet"].sum()
    pct = 100 * n_doublets / adata_sub.n_obs
    print(f"[{i+1}/{len(samples)}] {sample}: {n_doublets}/{adata_sub.n_obs} "
          f"({pct:.1f}%) [exp_rate={expected_rate:.4f}]")

adata_concat.obs["doublet_score_scrublet"] = doublet_scores
adata_concat.obs["doublet_scrublet"] = predicted_doublets

# Doublet-score violin by dataset → file
df = adata_concat.obs[["dataset", "doublet_score_scrublet"]].dropna()
plt.figure(figsize=(0.5 * df["dataset"].nunique() + 4, 5))
sns.violinplot(data=df, x="dataset", y="doublet_score_scrublet", inner="quartile",
               cut=0, scale="width", linewidth=0.8)
plt.xticks(rotation=90)
plt.ylabel("Scrublet doublet score")
plt.xlabel("Study")
plt.title("Scrublet doublet scores by study")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "02_scrublet_by_study.svg"), bbox_inches="tight")
plt.close("all")

# Filter doublets (score < 0.3; keep NaN/skipped samples)
print("Cells before doublet filter:", adata_concat.n_obs)
mask = (adata_concat.obs["doublet_score_scrublet"] < 0.3) | (adata_concat.obs["doublet_score_scrublet"].isna())
adata_concat = adata_concat[mask]
print("Cells after doublet filter :", adata_concat.n_obs)

# Re-apply low-cell sample filter post-doublet removal
sample_counts = adata_concat.obs["sample_num"].value_counts()
valid = sample_counts[sample_counts >= 500].index
adata_concat = adata_concat[adata_concat.obs["sample_num"].isin(valid)].copy()
print(f"Cells after re-filtering sample_num < 500: {adata_concat.n_obs}")

# ============================================================================
# FIX sample_id / sample_num SWAP + SAVE
# ============================================================================

adata_concat.obs["sample_id"], adata_concat.obs["sample_num"] = (
    adata_concat.obs["sample_num"].copy(),
    adata_concat.obs["sample_id"].copy(),
)
adata_concat.obs["sample_id"] = adata_concat.obs["sample_id"].astype(str)
adata_concat.obs["sample_num"] = adata_concat.obs["sample_num"].astype(str)

adata_concat.write_h5ad(FILTERED_PATH)
print(f"Done — wrote {FILTERED_PATH} (input to the integration trust boundary)")
