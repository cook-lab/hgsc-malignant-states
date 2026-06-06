#!/usr/bin/env python3
"""
Atlas 02 — Step 12a: fix zhang_2022 treatment_status metadata

PURPOSE
    Correct a metadata annotation error: zhang_2022 cells labelled "post-treatment"
    are confirmed platinum post-chemotherapy. Updates treatment_status
    "post-treatment" -> "post-chemotherapy" in all atlas h5ad files. MUST run
    before any metadata export (Supp Data 1) and before downstream proportions.

INPUTS / OUTPUTS (modified in place)
    obj("atlas_celltype_l2")  = hgsc_atlas_celltype_level2.h5ad
    obj("atlas_epithelial")   = hgsc_atlas_epithelial.h5ad
    DATA_ROOT/2026_final_atlas/hgsc_atlas_celltype_level1.h5ad

MANUSCRIPT PANEL(S)
    Metadata correctness; underpins Fig 2F/2G (treatment) and Supp Data 1.

RUNTIME TIER
    heavy (loads + rewrites each atlas file).
"""

import gc
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import obj, path, SEED  # noqa: E402

np.random.seed(SEED)

# ============================================================================
# CONFIG
# ============================================================================

H5AD_FILES = [
    obj("atlas_celltype_l2"),
    obj("atlas_epithelial"),
    path("data_root", "2026_final_atlas", "hgsc_atlas_celltype_level1.h5ad"),
]

OLD_VALUE = "post-treatment"
NEW_VALUE = "post-chemotherapy"
TARGET_STUDY = "zhang_2022"


def main():
    print("=" * 70)
    print("Step 12a — Fix zhang_2022 treatment_status")
    print(f"  {OLD_VALUE} → {NEW_VALUE}")
    print("=" * 70)

    total_changed = 0
    for fpath in H5AD_FILES:
        if not Path(fpath).exists():
            print(f"\n  SKIP: {fpath} (file not found)")
            continue

        print(f"\n  Processing: {fpath}")
        t0 = time.time()
        adata = ad.read_h5ad(fpath)
        print(f"    {adata.n_obs:,} cells loaded in {time.time()-t0:.0f}s")

        ts = adata.obs["treatment_status"].astype(str).copy()
        mask = (ts == OLD_VALUE) & (adata.obs["study"].astype(str) == TARGET_STUDY)
        n_affected = int(mask.sum())

        n_all_pt = int((ts == OLD_VALUE).sum())
        if n_all_pt != n_affected:
            print(f"    WARNING: {n_all_pt - n_affected} post-treatment cells are NOT "
                  f"from {TARGET_STUDY}. Only fixing {TARGET_STUDY} cells.")

        if n_affected == 0:
            print("    No cells to fix — skipping write")
            del adata
            gc.collect()
            continue

        ts.loc[mask] = NEW_VALUE
        adata.obs["treatment_status"] = pd.Categorical(ts, categories=sorted(ts.unique()))
        n_remaining = int((adata.obs["treatment_status"].astype(str) == OLD_VALUE).sum())
        n_new = int((adata.obs["treatment_status"].astype(str) == NEW_VALUE).sum())
        print(f"    Fixed: {n_affected:,} cells ({OLD_VALUE} → {NEW_VALUE})")
        print(f"    Remaining '{OLD_VALUE}': {n_remaining} | Total '{NEW_VALUE}': {n_new:,}")

        t1 = time.time()
        adata.write_h5ad(fpath)
        print(f"    Saved in {time.time()-t1:.0f}s")

        total_changed += n_affected
        del adata
        gc.collect()

    print(f"\n{'=' * 70}")
    print(f"DONE — {total_changed:,} cells updated across {len(H5AD_FILES)} files")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
