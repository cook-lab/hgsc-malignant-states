#!/usr/bin/env python3
"""
Supp Data 3 — Ciliated vs secretory differential expression
===========================================================
HGSC malignant-states atlas · supplemental table generator (authored in refactor).

PURPOSE
    Committed, schema-stable generator for the ciliated-vs-secretory Wilcoxon DE
    on the epithelial subset (the canonical logic lived only inside the 11v4
    recharacterization sandbox; this is the clean committed export). Reads the
    config entry-point epithelial subset, collapses level-2 epithelial labels to
    a Ciliated vs Secretory contrast, and runs scanpy rank_genes_groups.

INPUTS
    - obj("atlas_celltype_dir")/hgsc_atlas_final_epithelial.h5ad
      (uses log-normalized X; expects obs column 'celltype_level2')

OUTPUTS
    - supplemental/T3_epithelial_markers.csv
      Fixes legacy spelling typo (was T3_epithelial_markeres.csv).

MANUSCRIPT PANEL(S)
    Supp Data 3 (ciliated-vs-secretory DE).

RUNTIME TIER
    moderate (Wilcoxon DE over ~575k epithelial cells; two groups).
"""

import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

# --- central config (tables/ is 1 level below repo root) ---
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.config import obj, path, SEED  # noqa: E402

GROUPBY = "epithelial_class"   # derived Ciliated vs Secretory contrast
METHOD = "wilcoxon"


def main():
    np.random.seed(SEED)
    src = os.path.join(obj("atlas_celltype_dir"), "hgsc_atlas_final_epithelial.h5ad")
    print(f"[SD3] Reading {src}")
    adata = sc.read_h5ad(src)

    # Standardize legacy label, then collapse to the Ciliated-vs-Secretory contrast.
    lvl2 = (adata.obs["celltype_level2"].astype(str)
            .str.replace("Transitioning", "Intermediate", regex=False))
    adata.obs[GROUPBY] = np.where(
        lvl2.str.contains("Ciliated", case=False), "Ciliated", "Secretory")

    print(f"[SD3] rank_genes_groups ({METHOD}) Ciliated vs Secretory...")
    sc.tl.rank_genes_groups(adata, groupby=GROUPBY, method=METHOD, pts=True)

    de = sc.get.rank_genes_groups_df(adata, group=None)
    de = de.rename(columns={"group": "epithelial_class", "names": "gene"})
    de = de.sort_values(["epithelial_class", "pvals_adj", "scores"],
                        ascending=[True, True, False]).reset_index(drop=True)

    out = path("output_root", "supplemental", "T3_epithelial_markers.csv")
    de.to_csv(out, index=False)
    print(f"[SD3] Wrote {len(de):,} marker rows -> {out}")


if __name__ == "__main__":
    main()
