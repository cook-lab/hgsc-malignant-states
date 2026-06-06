#!/usr/bin/env python3
"""
Prepare NMF-based epithelial labels (canonical epitype schema)
==============================================================
HGSC malignant-states atlas backend.

Writes a `celltype_nmf` obs column into hgsc_atlas_final.h5ad that replaces the
6 Leiden-based epithelial celltype_level2 labels with the NMF-derived
SecA / Intermediate / SecB / Ciliated epitypes (the canonical manuscript
schema). Non-epithelial cells keep their existing celltype_level2 labels.

Schema (NMF Factor_2 percentile partition; thresholds on non-ciliated epi):
  - SecA:          Factor_2 < p50
  - Intermediate:  p50 <= Factor_2 < p75   (was "Transitioning" in originals)
  - SecB:          Factor_2 >= p75
  - Ciliated:      celltype_level2 == "Ciliated epithelial cell"

INPUTS:
  - hgsc_atlas_final.h5ad (obs: celltype_level1, celltype_level2)   [config: atlas_final]
  - output_root/03_epithelial_nmf/11d_nmf_usage.csv (Factor_2)

OUTPUTS:
  - hgsc_atlas_final.h5ad: new `celltype_nmf` obs column (written in place)
  - output_root/03_epithelial_nmf/celltype_nmf_mapping.csv  (canonical schema map)
  - output_root/03_epithelial_nmf/celltype_nmf_summary.csv

MANUSCRIPT PANELS: produces the canonical epitype schema consumed by Fig 5F,
  Fig 3H, Fig 7E/F/G, Supp Data 7 (via 17b/20/22 chains).

RUNTIME TIER: moderate (backed obs read + single obs column write).

SEEDING: deterministic (percentile thresholds; no stochastic step).

Usage:
    python 02_prepare_nmf_labels.py
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import obj, path  # noqa: E402

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────
ATLAS_H5AD = obj("atlas_final")
NMF_USAGE  = path("output_root", "03_epithelial_nmf", "11d_nmf_usage.csv")
OUT_DIR    = path("output_root", "03_epithelial_nmf")
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("  Prepare NMF-based epithelial labels")
print("=" * 60)

# ── 1. Load metadata ─────────────────────────────────────────
print("\n[1] Loading atlas metadata (backed mode)...")
adata = ad.read_h5ad(ATLAS_H5AD, backed="r")
obs = adata.obs[["celltype_level1", "celltype_level2"]].copy()
print(f"    Atlas: {len(obs):,} cells")

# ── 2. Load NMF Factor_2 ─────────────────────────────────────
print("\n[2] Loading NMF Factor_2 usage...")
nmf_usage = pd.read_csv(NMF_USAGE, index_col=0)
f2 = nmf_usage["Factor_2"]
print(f"    NMF data: {len(f2):,} epithelial cells")

# ── 3. Compute thresholds & assign labels ────────────────────
print("\n[3] Computing NMF thresholds and assigning labels...")

is_epithelial = obs["celltype_level1"] == "Epithelial"
epi_barcodes = obs.index[is_epithelial]
print(f"    Epithelial cells in atlas: {is_epithelial.sum():,}")

shared = epi_barcodes.intersection(f2.index)
print(f"    Matched with NMF data: {len(shared):,}")

f2_matched = f2.loc[shared]
is_ciliated = obs.loc[shared, "celltype_level2"] == "Ciliated epithelial cell"
non_cil_f2 = f2_matched[~is_ciliated]

p50 = np.percentile(non_cil_f2, 50)
p75 = np.percentile(non_cil_f2, 75)
print(f"    NMF thresholds: p50={p50:.4f}, p75={p75:.4f}")

epi_labels = pd.Series("SecA epithelium", index=shared)
epi_labels[f2_matched >= p50] = "Intermediate epithelium"
epi_labels[f2_matched >= p75] = "SecB epithelium"
epi_labels[is_ciliated] = "Ciliated epithelial cell"

for g in ["SecA epithelium", "Intermediate epithelium", "SecB epithelium",
          "Ciliated epithelial cell"]:
    n = (epi_labels == g).sum()
    print(f"    {g}: {n:,} ({n/len(epi_labels)*100:.1f}%)")

# ── 4. Build full celltype_nmf column ────────────────────────
print("\n[4] Building full celltype_nmf column...")

celltype_nmf = obs["celltype_level2"].astype(str).copy()
celltype_nmf.loc[shared] = epi_labels

epi_not_in_nmf = epi_barcodes.difference(shared)
if len(epi_not_in_nmf) > 0:
    print(f"    WARNING: {len(epi_not_in_nmf)} epithelial cells not in NMF data; "
          f"retaining their celltype_level2 labels")

print(f"    Total unique labels: {celltype_nmf.nunique()}")

# ── 5. Save mapping CSV ──────────────────────────────────────
print("\n[5] Saving mapping CSV...")
mapping = pd.DataFrame({
    "celltype_level2": obs["celltype_level2"],
    "celltype_nmf": celltype_nmf,
    "celltype_level1": obs["celltype_level1"],
})
mapping.to_csv(os.path.join(OUT_DIR, "celltype_nmf_mapping.csv"))
print(f"    Saved: celltype_nmf_mapping.csv ({len(mapping):,} rows)")

summary = celltype_nmf.value_counts().reset_index()
summary.columns = ["celltype_nmf", "n_cells"]
summary["pct"] = summary["n_cells"] / summary["n_cells"].sum() * 100
summary = summary.sort_values("n_cells", ascending=False)
summary.to_csv(os.path.join(OUT_DIR, "celltype_nmf_summary.csv"), index=False)
print("    Saved: celltype_nmf_summary.csv")

# ── 6. Write to h5ad ─────────────────────────────────────────
print("\n[6] Writing celltype_nmf to atlas h5ad...")
del adata  # close backed reader
adata = ad.read_h5ad(ATLAS_H5AD, backed="r+")
adata.obs["celltype_nmf"] = celltype_nmf.values
adata.file.close()
print("    Done — celltype_nmf column added to hgsc_atlas_final.h5ad")

print(f"\n{'='*60}")
print("  Step complete!")
print(f"  Output: {OUT_DIR}")
print(f"{'='*60}")
