#!/usr/bin/env python3
"""
Add xenium_celltype metacolumn to hgsc_atlas_celltype_level2.h5ad
=================================================================
Starts from celltype_level1 with three targeted splits:
  - Epithelial  → Ciliated epithelium  /  Secretory epithelium
  - DC          → Plasmacytoid dendritic cell  /  Conventional dendritic cell
  - T/NK cell   → NK cell  /  T cell
All other compartments keep their level-1 label unchanged.
"""

import os
import pandas as pd
import scanpy as sc

BASE     = "/Volumes/CookLab/Sarah/2026_final_atlas"
H5AD_IN  = os.path.join(BASE, "hgsc_atlas_celltype_level2.h5ad")
H5AD_OUT = os.path.join(BASE, "hgsc_atlas_xenium.h5ad")

# ── load ───────────────────────────────────────────────────────────────
print("Loading atlas...", flush=True)
adata = sc.read_h5ad(H5AD_IN)
print(f"  Shape: {adata.shape[0]:,} × {adata.shape[1]:,}")

l1 = adata.obs["celltype_level1"].astype(str)
l2 = adata.obs["celltype_level2"].astype(str)

# ── start from level 1 ────────────────────────────────────────────────
xenium = l1.copy()

# ── Epithelial: ciliated vs secretory ─────────────────────────────────
epi_mask = (l1 == "Epithelial")
cil_mask = epi_mask & (l2 == "Ciliated epithelial cell")
sec_mask = epi_mask & ~cil_mask

xenium.loc[cil_mask] = "Ciliated epithelium"
xenium.loc[sec_mask] = "Secretory epithelium"

print(f"\nEpithelial split:")
print(f"  Ciliated epithelium:   {cil_mask.sum():>9,}")
print(f"  Secretory epithelium:  {sec_mask.sum():>9,}")

# ── DC: pDC vs conventional ──────────────────────────────────────────
dc_mask  = (l1 == "DC")
pdc_mask = dc_mask & (l2 == "Plasmacytoid dendritic cell")
cdc_mask = dc_mask & ~pdc_mask

xenium.loc[pdc_mask] = "Plasmacytoid dendritic cell"
xenium.loc[cdc_mask] = "Conventional dendritic cell"

print(f"\nDC split:")
print(f"  Plasmacytoid DC:       {pdc_mask.sum():>9,}")
print(f"  Conventional DC:       {cdc_mask.sum():>9,}")

# ── T/NK: NK vs T ────────────────────────────────────────────────────
tnk_mask = (l1 == "T/NK cell")
nk_mask  = tnk_mask & l2.str.contains("NK", na=False)
t_mask   = tnk_mask & ~nk_mask

xenium.loc[nk_mask] = "NK cell"
xenium.loc[t_mask]  = "T cell"

print(f"\nT/NK split:")
print(f"  NK cell:               {nk_mask.sum():>9,}")
print(f"  T cell:                {t_mask.sum():>9,}")

# ── store ──────────────────────────────────────────────────────────────
LABEL_ORDER = [
    "Secretory epithelium",
    "Ciliated epithelium",
    "Mesothelial",
    "Fibroblast",
    "Smooth muscle",
    "Pericyte",
    "Endothelial",
    "T cell",
    "NK cell",
    "B cell",
    "Plasma cell",
    "Macrophage",
    "Conventional dendritic cell",
    "Plasmacytoid dendritic cell",
    "Neutrophil",
    "Mast cell",
]

adata.obs["xenium_celltype"] = pd.Categorical(xenium, categories=LABEL_ORDER)

# ── summary ────────────────────────────────────────────────────────────
print(f"\n{'='*55}")
print(f"xenium_celltype summary ({adata.obs['xenium_celltype'].nunique()} types)")
print(f"{'='*55}")
counts = adata.obs["xenium_celltype"].value_counts()
for ct in LABEL_ORDER:
    n = counts.get(ct, 0)
    pct = 100 * n / len(adata)
    print(f"  {ct:<35} {n:>9,}  ({pct:>5.1f}%)")

# ── save ───────────────────────────────────────────────────────────────
print(f"\nSaving → {H5AD_OUT}...", flush=True)
adata.write_h5ad(H5AD_OUT)
fsize = os.path.getsize(H5AD_OUT) / (1024**3)
print(f"  File size: {fsize:.1f} GB")
print(f"\n  Original {H5AD_IN} is UNTOUCHED.")
print("Done.")
