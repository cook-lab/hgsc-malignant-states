#!/usr/bin/env python3
"""
00b_extract_integration_umaps — cache integration UMAP coords + metadata to parquet
====================================================================================
PURPOSE
    Extract UMAP embeddings + obs from the integration-stage h5ads and the final
    atlas into lightweight parquet caches for the integration-quality
    supplemental figures. Run once.

INPUTS
    - <data_root>/2026_final_atlas/pre-processing and integration/
        20260213_integration/output/anndata/{integrated_harmony,
        integrated_scvi,integrated_scanvi}.h5ad   (integration trust boundary)
    - obj("atlas_scanvi")                          final atlas (UMAP + obs)

OUTPUTS (under output_root/_prep_caches/)
    - integration_harmony_umap.parquet
    - integration_scvi_umap.parquet
    - integration_scanvi_umap.parquet
    - atlas_final_umap.parquet

MANUSCRIPT PANEL(S): SF2B (study x cell-type bar), SF2C (per-study UMAPs).

RUNTIME TIER: moderate (backed reads of integration h5ads).
"""
from __future__ import annotations

import gc
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path  # noqa: E402

# ============================================================================
# PATHS
# ============================================================================

INTEG_DIR = obj("atlas_scanvi").replace(
    "hgsc_atlas_scanvi.h5ad",
    os.path.join("pre-processing and integration", "20260213_integration",
                 "output", "anndata"),
)

EXTRACTIONS = [
    {"h5ad": os.path.join(INTEG_DIR, "integrated_harmony.h5ad"),
     "cache": path("output_root", "_prep_caches", "integration_harmony_umap.parquet"),
     "label": "Harmony"},
    {"h5ad": os.path.join(INTEG_DIR, "integrated_scvi.h5ad"),
     "cache": path("output_root", "_prep_caches", "integration_scvi_umap.parquet"),
     "label": "scVI"},
    {"h5ad": os.path.join(INTEG_DIR, "integrated_scanvi.h5ad"),
     "cache": path("output_root", "_prep_caches", "integration_scanvi_umap.parquet"),
     "label": "scANVI (pre-filter)"},
    {"h5ad": obj("atlas_scanvi"),
     "cache": path("output_root", "_prep_caches", "atlas_final_umap.parquet"),
     "label": "Final atlas"},
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


if __name__ == "__main__":
    print("=" * 60)
    print("00b — Extract integration UMAPs to parquet")
    print("=" * 60)
    for spec in EXTRACTIONS:
        extract_umap(spec["h5ad"], spec["cache"], spec["label"])
    print("\nDone. Integration panel scripts can now read from _prep_caches parquets.")
