#!/usr/bin/env python3
"""
01_export_tma_barcode_patient_map — regenerate the TMA barcode->patient/type map
================================================================================
PURPOSE
    Regenerates the TMA barcode->patient/sample_type map consumed by
    figures/figure4/01 (Fig 4A/B composition dendrogram + donuts). Joins the
    deposited barcode->core map with the deposited core->patient map to attach
    each TMA cell's patient_id and sample_type. Run once.

INPUTS (deposited under data_root)
    - data_root/2026_final_xenium_analysis/output/figures/tma_core_patient_map.csv
        (columns: core_id, patient_id, sample_type)
    - data_root/2026_final_xenium_analysis/output/tma_barcode_core_map.csv
        (columns: barcode_orig, core_id, patient_id)

OUTPUTS
    - output_root/metadata/tma_barcode_patient_map.csv
        (columns: barcode_orig, patient_id, sample_type)

MANUSCRIPT PANEL(S): Fig 4A/B (TMA barcode->patient map) — consumed by
    figures/figure4/01_xenium_composition_by_tissue_dendro.py.

RUNTIME TIER: fast.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pandas as pd

# --- Central config (this file is 2 levels under the repo root) -------------
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# ============================================================================
# PATHS
# ============================================================================

CORE_MAP = path("data_root", "2026_final_xenium_analysis", "output",
                "figures", "tma_core_patient_map.csv")
BARCODE_CORE_MAP = path("data_root", "2026_final_xenium_analysis", "output",
                        "tma_barcode_core_map.csv")
OUT_CSV = path("output_root", "metadata", "tma_barcode_patient_map.csv")


def main():
    print("=" * 60)
    print("01 — Export TMA barcode->patient/sample_type map")
    print("=" * 60)

    print(f"  Reading core map: {CORE_MAP}", flush=True)
    core_map = pd.read_csv(CORE_MAP, usecols=["core_id", "patient_id", "sample_type"])
    core_map["core_id"] = core_map["core_id"].astype(str)

    print(f"  Reading barcode->core map: {BARCODE_CORE_MAP}", flush=True)
    bc_map = pd.read_csv(BARCODE_CORE_MAP, usecols=["barcode_orig", "core_id", "patient_id"])
    bc_map["core_id"] = bc_map["core_id"].astype(str)

    merged = bc_map.merge(core_map[["core_id", "sample_type"]], on="core_id", how="left")
    out = merged[["barcode_orig", "patient_id", "sample_type"]].copy()
    # Nullable integer so patient_id serializes as "2349" (not "2349.0"); any
    # unmatched barcode writes an empty cell that the consumer drops via notna().
    out["patient_id"] = out["patient_id"].astype("Int64")

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"    Saved: {OUT_CSV} ({len(out):,} barcodes)")


if __name__ == "__main__":
    main()
