#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/12a_fix_zhang2022_treatment_status.py, the script that ran.
# Kept: the in-place Zhang 2022 treatment_status fix (post-treatment -> post-chemotherapy). Removed: the closing pointer to
# follow-up scripts that are not part of the paper. The untrimmed script is in the lab project folder.
"""
Step 12a — Fix zhang_2022 treatment_status metadata
====================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Corrects a metadata annotation error: all cells from zhang_2022 that were
labelled "post-treatment" are confirmed as post-chemotherapy. This script
updates treatment_status from "post-treatment" → "post-chemotherapy" in
all atlas h5ad files.

Background:
    In the original metadata harmonisation, zhang_2022 samples were assigned
    "post-treatment" because the publication described them generically.
    Cross-referencing with the Methods section and GEO metadata confirms
    these patients received platinum-based chemotherapy (carboplatin/paclitaxel),
    making "post-chemotherapy" the correct label.

Files updated:
    1. hgsc_atlas_celltype_level2.h5ad    (29,934 cells affected)
    2. hgsc_atlas_epithelial.h5ad         (4,219 cells affected)
    3. hgsc_atlas_celltype_level1.h5ad    (29,934 cells affected)

NOTE: hgsc_atlas_celltype_level2_filtered.h5ad does not exist on disk
      (only a macOS ._ghost file remains); skipped.

After running this script, all downstream analyses (step 12 onward)
must be re-run.

Usage:
    python scripts/12a_fix_zhang2022_treatment_status.py
"""

import gc
import os
import time

import anndata as ad
import pandas as pd

# ============================================================================
# CONFIG
# ============================================================================

PROJECT_DIR = "/Volumes/CookLab/Sarah/2026_final_atlas"

H5AD_FILES = [
    "hgsc_atlas_celltype_level2.h5ad",
    "hgsc_atlas_epithelial.h5ad",
    "hgsc_atlas_celltype_level1.h5ad",
]

OLD_VALUE = "post-treatment"
NEW_VALUE = "post-chemotherapy"
TARGET_STUDY = "zhang_2022"

# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Step 12a — Fix zhang_2022 treatment_status")
    print(f"  {OLD_VALUE} → {NEW_VALUE}")
    print("=" * 70)

    total_changed = 0

    for fname in H5AD_FILES:
        path = os.path.join(PROJECT_DIR, fname)

        if not os.path.exists(path):
            print(f"\n  SKIP: {fname} (file not found)")
            continue

        print(f"\n  Processing: {fname}")
        t0 = time.time()

        # Read full object (need to modify and write back)
        print(f"    Loading...", flush=True)
        adata = ad.read_h5ad(path)
        print(f"    {adata.n_obs:,} cells loaded in {time.time()-t0:.0f}s")

        # Ensure treatment_status is a mutable string column
        ts = adata.obs["treatment_status"].astype(str).copy()

        # Identify cells to fix
        mask = (ts == OLD_VALUE) & (adata.obs["study"].astype(str) == TARGET_STUDY)
        n_affected = int(mask.sum())

        # Safety: verify ALL post-treatment cells are from zhang_2022
        mask_all_pt = ts == OLD_VALUE
        n_all_pt = int(mask_all_pt.sum())
        if n_all_pt != n_affected:
            other = n_all_pt - n_affected
            print(f"    WARNING: {other} post-treatment cells are NOT from "
                  f"{TARGET_STUDY}. Only fixing {TARGET_STUDY} cells.")

        if n_affected == 0:
            print(f"    No cells to fix — skipping write")
            del adata
            gc.collect()
            continue

        # Apply fix
        ts.loc[mask] = NEW_VALUE

        # Rebuild as categorical with updated categories
        new_cats = sorted(ts.unique())
        adata.obs["treatment_status"] = pd.Categorical(ts, categories=new_cats)

        # Verify
        n_remaining = int((adata.obs["treatment_status"].astype(str) == OLD_VALUE).sum())
        n_new = int((adata.obs["treatment_status"].astype(str) == NEW_VALUE).sum())
        print(f"    Fixed: {n_affected:,} cells  ({OLD_VALUE} → {NEW_VALUE})")
        print(f"    Remaining '{OLD_VALUE}': {n_remaining}")
        print(f"    Total '{NEW_VALUE}': {n_new:,}")

        # Write back
        print(f"    Writing...", flush=True)
        t1 = time.time()
        adata.write_h5ad(path)
        print(f"    Saved in {time.time()-t1:.0f}s")

        total_changed += n_affected
        del adata
        gc.collect()

    print(f"\n{'=' * 70}")
    print(f"DONE — {total_changed:,} cells updated across {len(H5AD_FILES)} files")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
