#!/usr/bin/env python
"""
Atlas 01 — Step 02: aggregate per-study h5ads into the raw atlas

PURPOSE
    Read the per-study raw-count h5ads (Vazquez-Garcia, Luo, Zheng, and a
    pre-merged aggregate), harmonise each to the common gene set, force sparse
    int32 counts, concatenate (inner join), drop samples with <500 cells, and
    write the raw atlas object (full gene space, counts in .X).

INPUTS
    raw_datasets/vazquez_garcia_2022.h5ad
    raw_datasets/luo_2024.h5ad
    raw_datasets/zheng_2023.h5ad
    raw_datasets/ovarian_cancer_aggregate.h5ad
        (per-study raw inputs; REQUIRED but not deposited — see config key raw_datasets)

OUTPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_raw.h5ad   (raw counts, full gene space)

RUNTIME TIER
    heavy (loads + concatenates multi-study raw matrices; high-memory CPU).

MANUSCRIPT ROLE
    Pre-integration aggregation. Raw substrate for the scVI/scANVI integration;
    no panel rendered directly.

NOTE
    AUTHORITATIVE original cluster script (atlas_00_aggregate.py). Analytical logic
    preserved EXACTLY (common-gene intersection, sparse int32 coercion, <500-cell
    sample filter). Only the hardcoded ComputeCanada cluster paths
    (/project/6090753/dcook/...) were replaced with central-config roots.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

import scanpy as sc, anndata as ad, gc, os, tempfile
import seaborn as sns
import numpy as np
import torch
import scipy.sparse as sp

# Config
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

# Resolve per-study raw inputs via the central raw_datasets root (config key raw_datasets)
from config.config import CFG, DATA_ROOT  # noqa: E402


def raw_dataset(fname):
    return str(DATA_ROOT / CFG["paths"]["raw_datasets"] / fname)

OUT_RAW = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_raw.h5ad")

# Read Vazquez-Garcia atlas
print("Reading Vazquez-Garcia", flush=True)
adata_vg = sc.read_h5ad(raw_dataset('vazquez_garcia_2022.h5ad'))

# Read Luo atlas
print("Reading Luo et al", flush=True)
adata_luo = sc.read_h5ad(raw_dataset('luo_2024.h5ad'))

## Adding single-value metadata that got dropped in export
adata_luo.obs['hist_subtype'] = "High-grade serous"
adata_luo.obs['source'] = "Luo_2024"
adata_luo.obs['method'] = "10x"
adata_luo.obs['enrichment'] = "Whole tissue"
adata_luo.obs['accession'] = "GSE222556"
adata_luo.obs['PFI'] = "Unknown"
adata_luo.obs['mut_profile'] = "Unknown"
adata_luo.obs['orig.ident'] = "Luo_2024"

# Read Zheng 2023
print("Reading Zheng et al", flush=True)
adata_zheng = sc.read_h5ad(raw_dataset('zheng_2023.h5ad'))
adata_zheng.obs['hist_subtype'] = "High-grade serous"
adata_zheng.obs['source'] = "Zheng_2023"
adata_zheng.obs['method'] = "10x"
adata_zheng.obs['enrichment'] = "Whole tissue"
adata_zheng.obs['accession'] = "10.17632/rc47y6m9mp.2"
adata_zheng.obs['PFI'] = "Unknown"
adata_zheng.obs['mut_profile'] = "Unknown"
adata_zheng.obs['orig.ident'] = "Zheng_2023"

# Reading other pre-merged datasets
print("Reading aggregate data sets", flush=True)
adata_agg = sc.read_h5ad(raw_dataset('ovarian_cancer_aggregate.h5ad'))

# Make age categorical
adata_agg.obs['age'] = adata_agg.obs['age'].astype('str').astype('category')
adata_vg.obs['age'] = adata_vg.obs['age'].astype('str').astype('category')
adata_luo.obs['age'] = adata_luo.obs['age'].astype('str').astype('category')

print("Merging data", flush=True)
adata_list = [adata_vg, adata_luo, adata_zheng, adata_agg]

# 1️⃣ Slim down each AnnData before concat
# ---- 0. Determine the common gene set & a reference order ----
common_genes = set(adata_list[0].var_names)
for a in adata_list[1:]:
    common_genes &= set(a.var_names)
common_genes = np.array(sorted(common_genes))          # keep alphabetical order

# ---- 1. Harmonise each object BEFORE concat ----
for a in adata_list:
    # a) keep only common genes **and reorder to reference order**
    a = a[:, common_genes].copy()
    # b) force sparse int32 counts
    if not sp.issparse(a.X):
        a.X = sp.csr_matrix(a.X)
    if a.X.dtype != np.int32:
        a.X.data = a.X.data.astype(np.int32)
    # c) drop heavy slots
    a.raw = None
    a.layers.clear()

# ---- 2. Now concat is just sparse vstack ----
print("🔹  Merging data …", flush=True)
adata = sc.concat(
    adata_list,
    join="inner",     # columns already match identically
)

del adata_vg, adata_agg, adata_luo
gc.collect()

# Remove samples with <500 cells
sample_counts = adata.obs['sample_id'].value_counts()
valid_ids = sample_counts[sample_counts>500].index
adata = adata[adata.obs["sample_id"].isin(valid_ids)].copy()

# obs['age'] keeps getting converted to the dtype object for some reason
#adata.obs['age'] = adata.obs['age'].astype('category')
adata.obs["age"] = adata.obs["age"].astype(str)

# ---------- Save Raw ----------
# ensure X is sparse
if not sp.issparse(adata.X):
    adata.X = sp.csr_matrix(adata.X)

# write out full raw counts + obs/var
print("Saving raw...", flush=True)
sc.write(
    OUT_RAW,
    adata,
    compression="gzip"
)
