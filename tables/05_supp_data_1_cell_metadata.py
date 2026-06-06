#!/usr/bin/env python3
"""
Supp Data 1 — Cell-level metadata export
========================================
HGSC malignant-states atlas · supplemental table generator (authored in refactor).

PURPOSE
    Export the per-cell metadata table (obs) for the full integrated atlas as a
    committed, schema-stable generator (the canonical analysis only ever produced
    this via ad-hoc interactive export). Reads the config entry-point object so
    the export tracks whatever the deposited final atlas contains. Must run AFTER
    the Zhang-2022 treatment_status fix (canonical step 12a) is baked into
    hgsc_atlas_final.h5ad.

INPUTS
    - obj("atlas_final")  -> 2026_final_atlas/hgsc_atlas_final.h5ad (obs only; backed read)

OUTPUTS
    - supplemental/T1_atlas_metadata.csv    (one row per cell; n ~ 1,980,703)
      Fixes legacy filename: no double extension.

MANUSCRIPT PANEL(S)
    Supp Data 1 (cell-level metadata).

RUNTIME TIER
    moderate (obs is large; counts matrix is never materialized — backed read).
"""

import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# --- central config (tables/ is 1 level below repo root) ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.config import obj, path  # noqa: E402

# Columns to export, in order. Anything absent in obs is skipped (warned).
EXPORT_COLS = [
    "celltype_level1", "celltype_level2", "celltype_nmf",
    "study", "patient_id", "sample_id",
    "anatomic_site", "metastatic_site", "treatment_status",
    "tp53_status", "hrd_status", "brca_status",
    "total_counts", "n_genes_by_counts", "pct_counts_mt",
    "doublet_score", "scrublet_prediction",
]


def main():
    src = obj("atlas_final")
    print(f"[SD1] Reading obs (backed) from {src}")
    adata = ad.read_h5ad(src, backed="r")
    obs = adata.obs.copy()
    adata.file.close()

    # Standardize legacy epithelial label everywhere in the metadata.
    obs = obs.apply(lambda c: c.astype(str).str.replace(
        "Transitioning", "Intermediate", regex=False) if c.dtype == object or
        str(c.dtype) == "category" else c)

    present = [c for c in EXPORT_COLS if c in obs.columns]
    missing = [c for c in EXPORT_COLS if c not in obs.columns]
    if missing:
        print(f"[SD1] WARNING: columns absent in obs, skipped: {missing}")

    out_df = obs[present].copy()
    out_df.insert(0, "cell_barcode", obs.index.values)

    out = path("output_root", "supplemental", "T1_atlas_metadata.csv")
    out_df.to_csv(out, index=False)
    print(f"[SD1] Wrote {len(out_df):,} cells x {out_df.shape[1]} cols -> {out}")


if __name__ == "__main__":
    main()
