#!/usr/bin/env python3
"""
00b — Extract UMAP coordinates + metadata from integration h5ads
================================================================

Purpose
    Cache lightweight parquets (UMAP coords + study + celltype_pred + QC) for
    the integration-quality supplementary figures (SF2B/C). Run once.

INPUTS  (via config)
    Pre-integration comparison objects (cluster job outputs; trust boundary):
      data_root/2026_final_atlas/pre-processing and integration/
        20260213_integration/output/anndata/{integrated_harmony,integrated_scvi,
        integrated_scanvi}.h5ad
    Final filtered atlas: obj("atlas_scanvi")

OUTPUTS
    output_root/figures/data/integration_harmony_umap.parquet
    output_root/figures/data/integration_scvi_umap.parquet
    output_root/figures/data/integration_scanvi_umap.parquet
    output_root/figures/data/atlas_final_umap.parquet

MANUSCRIPT PANEL(S)
    Upstream cache for SF2B (study composition) and SF2C (per-study UMAPs).
    (The harmony/scvi parquets support the integration method comparison;
    SF2B/C consume atlas_final_umap.parquet.)

RUNTIME TIER
    moderate (backed reads of several integration h5ads).
"""

import gc
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path  # noqa: E402

# ============================================================================
# PATHS (central config)
# ============================================================================

INTEG_DIR = path(
    "data_root", "2026_final_atlas", "pre-processing and integration",
    "20260213_integration", "output", "anndata",
)

EXTRACTIONS = [
    {
        "h5ad": os.path.join(INTEG_DIR, "integrated_harmony.h5ad"),
        "cache": path("output_root", "figures", "data", "integration_harmony_umap.parquet"),
        "label": "Harmony",
    },
    {
        "h5ad": os.path.join(INTEG_DIR, "integrated_scvi.h5ad"),
        "cache": path("output_root", "figures", "data", "integration_scvi_umap.parquet"),
        "label": "scVI",
    },
    {
        "h5ad": os.path.join(INTEG_DIR, "integrated_scanvi.h5ad"),
        "cache": path("output_root", "figures", "data", "integration_scanvi_umap.parquet"),
        "label": "scANVI (pre-filter)",
    },
    {
        "h5ad": obj("atlas_scanvi"),
        "cache": path("output_root", "figures", "data", "atlas_final_umap.parquet"),
        "label": "Final atlas",
    },
]

OBS_COLS = ["study", "celltype_pred", "n_genes_by_counts", "total_counts"]


def extract_umap(h5ad_path, cache_path, label):
    if os.path.exists(cache_path):
        print(f"  SKIP {label} — cache exists: {cache_path}")
        return

    import scanpy as sc

    print(f"  Loading {label} (backed): {h5ad_path}", flush=True)
    adata = sc.read_h5ad(h5ad_path, backed="r")
    print(f"    Shape: {adata.shape[0]:,} x {adata.shape[1]:,}")

    cols = [c for c in OBS_COLS if c in adata.obs.columns]
    df = adata.obs[cols].copy()

    umap = np.array(adata.obsm["X_umap"])
    df["UMAP1"] = umap[:, 0]
    df["UMAP2"] = umap[:, 1]

    for col in ["n_genes_by_counts", "total_counts"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    adata.file.close()
    del adata, umap
    gc.collect()

    df.to_parquet(cache_path)
    print(f"    Saved: {cache_path} ({len(df):,} cells)")


print("=" * 60)
print("00b — Extract integration UMAPs to parquet")
print("=" * 60)

for spec in EXTRACTIONS:
    extract_umap(spec["h5ad"], spec["cache"], spec["label"])

print("\nDone. SF2B/C can now read from output_root/figures/data/*.parquet")
