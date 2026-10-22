#!/usr/bin/env python3
# Trimmed copy (2026-10-21) of 202601_final_atlas/scripts/22a_cibersort_reference_v2.py, the script that ran.
# Kept: the 300-cells-per-group atlas reference, its phenotype labels and the TCGA mixture link, which the BayesPrism step reads.
# Removed: the gene-check prints, the fallback signature matrix for the CIBERSORTx portal and the closing upload instructions.
# The untrimmed script is in the lab project folder.
"""
Step 22a v2 — Lean CIBERSORTx reference for TCGA deconvolution
================================================================
HGSC Single-Cell Atlas  ·  Cook Lab  ·  2026

Slimmed version: 300 cells/group, genes filtered to those expressed
in >=10% of cells in at least one group. Produces a much smaller
reference file that CIBERSORTx web portal can handle.

Also exports a pre-built signature matrix (group means) as fallback
for CIBERSORTx "custom" mode if scRNA-seq mode still fails.

Output (output/22_tcga_deconvolution/):
    - cibersortx_sc_reference_v2.txt     (lean scRNA-seq reference)
    - cibersortx_phenotypes_v2.txt       (phenotype labels)
    - cibersortx_sig_matrix.txt          (pre-built signature matrix, fallback)
    - cibersortx_mixture.txt             (TCGA, symlinked)

Usage:
    python scripts/22a_cibersort_reference_v2.py
"""

import os
import warnings
import gc

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
PROJECT_DIR  = "/Volumes/CookLab/Sarah/2026_final_atlas"
NMF_MAP      = os.path.join(PROJECT_DIR, "output", "17_nmf_labels",
                             "celltype_nmf_mapping.csv")
CELLTYPE_DIR = os.path.join(PROJECT_DIR, "celltype_h5ad")
TCGA_EXPR    = os.path.join(PROJECT_DIR, "data", "cibersort_data_prev",
                             "tcga_ecotyper.txt")
OUT_DIR      = os.path.join(PROJECT_DIR, "output", "22_tcga_deconvolution")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────
CELLS_PER_GROUP   = 300      # Reduced from 500
MIN_PCT_EXPRESSED = 0.10     # Gene must be expressed in >=10% of cells
                             # in at least one group
SEED = 42

# ── Functional grouping ──────────────────────────────────────
LEVEL1_TO_GROUP = {
    "Epithelial":    None,
    "T/NK cell":     "T_NK",
    "B cell":        "B_Plasma",
    "Plasma cell":   "B_Plasma",
    "Macrophage":    "Macrophage",
    "DC":            "DC",
    "Fibroblast":    "Fibroblast_Stromal",
    "Smooth muscle": "Fibroblast_Stromal",
    "Pericyte":      "Fibroblast_Stromal",
    "Endothelial":   "Endothelial",
    "Mesothelial":   "Mesothelial",
    "Neutrophil":    "Other_immune",
    "Mast cell":     "Other_immune",
}

NMF_EPI_TO_GROUP = {
    "SecA epithelium":          "SecA_epithelium",
    "Transitioning epithelium": "Transitioning_epithelium",
    "SecB epithelium":          "SecB_epithelium",
    "Ciliated epithelial cell": "Ciliated_epithelium",
}

H5AD_FILES = {
    "Epithelial":    "hgsc_atlas_final_epithelial.h5ad",
    "T/NK cell":     "hgsc_atlas_final_t_nk_cell.h5ad",
    "B cell":        "hgsc_atlas_final_b_cell.h5ad",
    "Plasma cell":   "hgsc_atlas_final_plasma_cell.h5ad",
    "Macrophage":    "hgsc_atlas_final_macrophage.h5ad",
    "DC":            "hgsc_atlas_final_dc.h5ad",
    "Fibroblast":    "hgsc_atlas_final_fibroblast.h5ad",
    "Smooth muscle": "hgsc_atlas_final_smooth_muscle.h5ad",
    "Pericyte":      "hgsc_atlas_final_pericyte.h5ad",
    "Endothelial":   "hgsc_atlas_final_endothelial.h5ad",
    "Mesothelial":   "hgsc_atlas_final_mesothelial.h5ad",
    "Neutrophil":    "hgsc_atlas_final_neutrophil.h5ad",
    "Mast cell":     "hgsc_atlas_final_mast_cell.h5ad",
}

# ============================================================================
print("=" * 65)
print("  Step 22a v2 — Lean CIBERSORTx reference")
print("=" * 65)

# ============================================================================
# 1. LOAD NMF MAPPING & BUILD GROUP ASSIGNMENTS
# ============================================================================
print("\n[1] Loading NMF mapping...")

nmf_map = pd.read_csv(NMF_MAP, index_col=0)

def assign_group(row):
    level1 = row["celltype_level1"]
    nmf    = row["celltype_nmf"]
    if level1 == "Epithelial":
        return NMF_EPI_TO_GROUP.get(nmf, None)
    return LEVEL1_TO_GROUP.get(level1, None)

nmf_map["deconv_group"] = nmf_map.apply(assign_group, axis=1)
nmf_map = nmf_map.dropna(subset=["deconv_group"])

group_counts = nmf_map["deconv_group"].value_counts().sort_index()
print(f"    {len(nmf_map):,} cells across {len(group_counts)} groups")

# ============================================================================
# 2. DOWNSAMPLE
# ============================================================================
print(f"\n[2] Downsampling to {CELLS_PER_GROUP} cells/group...")

selected_barcodes = {}
for group in sorted(group_counts.index):
    group_cells = nmf_map[nmf_map["deconv_group"] == group]
    n = min(CELLS_PER_GROUP, len(group_cells))
    selected = group_cells.sample(n=n, random_state=SEED).index.tolist()
    selected_barcodes[group] = selected
    print(f"    {group:<30s} {n:>5}")

total_selected = sum(len(v) for v in selected_barcodes.values())
print(f"    Total: {total_selected:,} cells")

# ============================================================================
# 3. LOAD RAW COUNTS
# ============================================================================
print("\n[3] Loading raw counts...")

barcode_to_group = {}
for group, barcodes in selected_barcodes.items():
    for bc in barcodes:
        barcode_to_group[bc] = group

selected_nmf = nmf_map.loc[list(barcode_to_group.keys())]
compartments_needed = selected_nmf["celltype_level1"].unique()

all_counts = []
all_barcodes = []
all_groups = []
gene_names = None

for level1 in sorted(compartments_needed):
    h5ad_file = os.path.join(CELLTYPE_DIR, H5AD_FILES[level1])
    if not os.path.exists(h5ad_file):
        continue

    comp_barcodes = selected_nmf[
        selected_nmf["celltype_level1"] == level1
    ].index.tolist()
    if not comp_barcodes:
        continue

    print(f"    {level1} ({len(comp_barcodes)} cells)...", end="")

    adata = ad.read_h5ad(h5ad_file, backed="r")
    obs_set = set(adata.obs.index)
    matching = [bc for bc in comp_barcodes if bc in obs_set]

    if not matching:
        print(" no matches")
        adata.file.close()
        continue

    idx_map = {bc: i for i, bc in enumerate(adata.obs.index)}
    indices = sorted([idx_map[bc] for bc in matching])

    adata_sub = adata[indices, :].to_memory()
    adata.file.close()

    counts = adata_sub.layers.get("counts", adata_sub.X)
    if sp.issparse(counts):
        counts = counts.tocsr()

    if gene_names is None:
        gene_names = adata_sub.var_names.tolist()

    sub_barcodes = adata_sub.obs.index.tolist()
    sub_groups = [barcode_to_group[bc] for bc in sub_barcodes]

    all_counts.append(counts)
    all_barcodes.extend(sub_barcodes)
    all_groups.extend(sub_groups)

    print(f" {len(sub_barcodes)} cells")
    del adata_sub; gc.collect()

# Stack
print(f"\n    Stacking...")
if all(sp.issparse(c) for c in all_counts):
    combined_counts = sp.vstack(all_counts, format="csr")
else:
    combined_counts = np.vstack([
        c.toarray() if sp.issparse(c) else c for c in all_counts
    ])

if sp.issparse(combined_counts):
    combined_dense = combined_counts.toarray()
else:
    combined_dense = combined_counts

combined_dense = np.rint(combined_dense).astype(int)
print(f"    Raw matrix: {combined_dense.shape[0]} cells x {combined_dense.shape[1]} genes")

# ============================================================================
# 4. FILTER GENES — keep those expressed in >=10% of cells in any group
# ============================================================================
print(f"\n[4] Filtering genes (>={MIN_PCT_EXPRESSED*100:.0f}% expressed in any group)...")

groups_array = np.array(all_groups)
unique_groups = sorted(set(all_groups))

keep_gene = np.zeros(len(gene_names), dtype=bool)
for group in unique_groups:
    mask = groups_array == group
    group_counts_mat = combined_dense[mask, :]
    pct_expressed = np.mean(group_counts_mat > 0, axis=0)
    keep_gene |= (pct_expressed >= MIN_PCT_EXPRESSED)

n_kept = keep_gene.sum()
n_total = len(gene_names)
print(f"    Genes kept: {n_kept:,} / {n_total:,} ({100*n_kept/n_total:.1f}%)")

# Filter
filtered_genes = [g for g, k in zip(gene_names, keep_gene) if k]
filtered_counts = combined_dense[:, keep_gene]

# ============================================================================
# 5. EXPORT scRNA-seq REFERENCE (lean)
# ============================================================================
print(f"\n[5] Exporting lean scRNA-seq reference...")
print(f"    Matrix: {len(filtered_genes)} genes x {len(all_barcodes)} cells")

ref_path = os.path.join(OUT_DIR, "cibersortx_sc_reference_v2.txt")
with open(ref_path, "w") as f:
    f.write("Gene\t" + "\t".join(all_barcodes) + "\n")
    for i, gene in enumerate(filtered_genes):
        gene_idx = gene_names.index(gene)
        row_vals = combined_dense[:, gene_idx]
        f.write(gene + "\t" + "\t".join(str(v) for v in row_vals) + "\n")

file_size_mb = os.path.getsize(ref_path) / 1e6
print(f"    Saved: cibersortx_sc_reference_v2.txt ({file_size_mb:.1f} MB)")

# ============================================================================
# 6. EXPORT PHENOTYPE FILE
# ============================================================================
print("\n[6] Exporting phenotype file...")

pheno_path = os.path.join(OUT_DIR, "cibersortx_phenotypes_v2.txt")
with open(pheno_path, "w") as f:
    for bc, group in zip(all_barcodes, all_groups):
        f.write(f"{bc}\t{group}\n")

print(f"    Saved: cibersortx_phenotypes_v2.txt")

# ============================================================================
# 8. ENSURE MIXTURE FILE
# ============================================================================
mixture_path = os.path.join(OUT_DIR, "cibersortx_mixture.txt")
if not os.path.exists(mixture_path):
    os.symlink(os.path.abspath(TCGA_EXPR), mixture_path)
    print(f"\n[8] Symlinked mixture file")
else:
    print(f"\n[8] Mixture file already exists")
