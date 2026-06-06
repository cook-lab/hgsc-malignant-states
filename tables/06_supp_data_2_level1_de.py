#!/usr/bin/env python3
"""
Supp Data 2 — Level-1 cell-type differential expression markers
===============================================================
HGSC malignant-states atlas · supplemental table generator (authored in refactor).

PURPOSE
    Committed, schema-stable generator for the one-vs-rest Wilcoxon marker table
    across level-1 cell types (the canonical analysis produced this via ad-hoc
    interactive export with no committed script). Reads the config entry-point
    object and runs scanpy rank_genes_groups, then tidies to a long-form table.

INPUTS
    - obj("atlas_final")  -> hgsc_atlas_final.h5ad
      (uses log-normalized X; expects obs column 'celltype_level1')

OUTPUTS
    - supplemental/T2_celltypelevel1_markers.csv
      Fixes legacy double-extension filename (was T2_*.csv.csv).

MANUSCRIPT PANEL(S)
    Supp Data 2 (level-1 DE markers).

RUNTIME TIER
    heavy (Wilcoxon DE over ~2M cells x level-1 groups).
"""

import sys
from pathlib import Path

import numpy as np
import scanpy as sc

# --- central config (tables/ is 1 level below repo root) ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.config import obj, path, SEED  # noqa: E402

GROUPBY = "celltype_level1"
METHOD = "wilcoxon"


def main():
    np.random.seed(SEED)
    src = obj("atlas_final")
    print(f"[SD2] Reading {src}")
    adata = sc.read_h5ad(src)

    # Standardize legacy epithelial label in the grouping variable.
    adata.obs[GROUPBY] = (adata.obs[GROUPBY].astype(str)
                          .str.replace("Transitioning", "Intermediate", regex=False))

    print(f"[SD2] rank_genes_groups ({METHOD}, one-vs-rest) on '{GROUPBY}'...")
    sc.tl.rank_genes_groups(adata, groupby=GROUPBY, method=METHOD, pts=True)

    de = sc.get.rank_genes_groups_df(adata, group=None)
    de = de.rename(columns={"group": "celltype_level1", "names": "gene"})
    de = de.sort_values(["celltype_level1", "pvals_adj", "scores"],
                        ascending=[True, True, False]).reset_index(drop=True)

    out = path("output_root", "supplemental", "T2_celltypelevel1_markers.csv")
    de.to_csv(out, index=False)
    print(f"[SD2] Wrote {len(de):,} marker rows -> {out}")


if __name__ == "__main__":
    main()
