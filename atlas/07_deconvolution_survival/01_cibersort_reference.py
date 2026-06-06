#!/usr/bin/env python3
"""
Lean CIBERSORTx / BayesPrism scRNA-seq reference for TCGA deconvolution
======================================================================
HGSC malignant-states atlas backend.

Builds a lean single-cell reference (300 cells/group; genes expressed in >=10%
of cells in at least one group) plus a pre-built signature matrix, for TCGA-OV
bulk deconvolution. Epithelial cells are split into NMF epitypes
(SecA / Intermediate / SecB / Ciliated); TME compartments are functionally
grouped. The reference is consumed by the BayesPrism runs (02/03).

INPUTS:
  - output_root/03_epithelial_nmf/celltype_nmf_mapping.csv  (from 02_prepare_nmf_labels.py)
  - <data_root>/2026_final_atlas/celltype_h5ad/*.h5ad  (per-compartment subsets)
  - <data_root>/2026_final_atlas/data/cibersort_data_prev/tcga_ecotyper.txt (gene overlap + mixture)

OUTPUTS (output_root/07_deconvolution_survival/):
  - cibersortx_sc_reference_v2.txt   (lean scRNA-seq reference)
  - cibersortx_phenotypes_v2.txt
  - cibersortx_sig_matrix.txt        (pre-built signature matrix)
  - cibersortx_mixture.txt           (TCGA mixture, symlink)

MANUSCRIPT PANELS: upstream of Fig 7E/F/G (TCGA deconvolution + survival).

RUNTIME TIER: moderate (loads per-compartment subsets, downsampled).

SEEDING: downsample uses random_state=SEED for determinism.

Usage:
    python 01_cibersort_reference.py
"""

import os
import sys
import warnings
import gc

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from config.config import path, SEED  # noqa: E402

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────
NMF_MAP      = path("output_root", "03_epithelial_nmf", "celltype_nmf_mapping.csv")
CELLTYPE_DIR = path("data_root", "2026_final_atlas", "celltype_h5ad")
TCGA_EXPR    = path("data_root", "2026_final_atlas", "data", "cibersort_data_prev",
                    "tcga_ecotyper.txt")
OUT_DIR      = path("output_root", "07_deconvolution_survival")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Parameters (validated) ────────────────────────────────────
CELLS_PER_GROUP   = 300
MIN_PCT_EXPRESSED = 0.10

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
    "Intermediate epithelium":  "Intermediate_epithelium",
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

print("=" * 65)
print("  Lean CIBERSORTx / BayesPrism reference")
print("=" * 65)

# ── 1. Load NMF mapping & build group assignments ────────────
print("\n[1] Loading NMF mapping...")
nmf_map = pd.read_csv(NMF_MAP, index_col=0)


def assign_group(row):
    if row["celltype_level1"] == "Epithelial":
        return NMF_EPI_TO_GROUP.get(row["celltype_nmf"], None)
    return LEVEL1_TO_GROUP.get(row["celltype_level1"], None)


nmf_map["deconv_group"] = nmf_map.apply(assign_group, axis=1)
nmf_map = nmf_map.dropna(subset=["deconv_group"])
group_counts = nmf_map["deconv_group"].value_counts().sort_index()
print(f"    {len(nmf_map):,} cells across {len(group_counts)} groups")

# ── 2. Downsample ────────────────────────────────────────────
print(f"\n[2] Downsampling to {CELLS_PER_GROUP} cells/group...")
selected_barcodes = {}
for group in sorted(group_counts.index):
    group_cells = nmf_map[nmf_map["deconv_group"] == group]
    n = min(CELLS_PER_GROUP, len(group_cells))
    selected_barcodes[group] = group_cells.sample(n=n, random_state=SEED).index.tolist()
    print(f"    {group:<30s} {n:>5}")
total_selected = sum(len(v) for v in selected_barcodes.values())
print(f"    Total: {total_selected:,} cells")

# ── 3. Load raw counts ───────────────────────────────────────
print("\n[3] Loading raw counts...")
barcode_to_group = {bc: g for g, bcs in selected_barcodes.items() for bc in bcs}
selected_nmf = nmf_map.loc[list(barcode_to_group.keys())]
compartments_needed = selected_nmf["celltype_level1"].unique()

all_counts, all_barcodes, all_groups = [], [], []
gene_names = None
for level1 in sorted(compartments_needed):
    h5ad_file = os.path.join(CELLTYPE_DIR, H5AD_FILES[level1])
    if not os.path.exists(h5ad_file):
        continue
    comp_barcodes = selected_nmf[selected_nmf["celltype_level1"] == level1].index.tolist()
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
    all_counts.append(counts)
    all_barcodes.extend(sub_barcodes)
    all_groups.extend([barcode_to_group[bc] for bc in sub_barcodes])
    print(f" {len(sub_barcodes)} cells")
    del adata_sub
    gc.collect()

print("\n    Stacking...")
if all(sp.issparse(c) for c in all_counts):
    combined_counts = sp.vstack(all_counts, format="csr")
else:
    combined_counts = np.vstack([c.toarray() if sp.issparse(c) else c for c in all_counts])
combined_dense = combined_counts.toarray() if sp.issparse(combined_counts) else combined_counts
combined_dense = np.rint(combined_dense).astype(int)
print(f"    Raw matrix: {combined_dense.shape[0]} cells x {combined_dense.shape[1]} genes")

# ── 4. Filter genes ──────────────────────────────────────────
print(f"\n[4] Filtering genes (>={MIN_PCT_EXPRESSED*100:.0f}% expressed in any group)...")
groups_array = np.array(all_groups)
unique_groups = sorted(set(all_groups))
keep_gene = np.zeros(len(gene_names), dtype=bool)
for group in unique_groups:
    mask = groups_array == group
    pct_expressed = np.mean(combined_dense[mask, :] > 0, axis=0)
    keep_gene |= (pct_expressed >= MIN_PCT_EXPRESSED)
n_kept = keep_gene.sum()
print(f"    Genes kept: {n_kept:,} / {len(gene_names):,} ({100*n_kept/len(gene_names):.1f}%)")

filtered_genes = [g for g, k in zip(gene_names, keep_gene) if k]

key_genes = ["SOX17", "MECOM", "FBXO21", "PBX1", "WT1", "PAX8",
             "KRT17", "KRT19", "KRT7", "TACSTD2", "SLPI", "LCN2", "MMP7"]
for g in key_genes:
    print(f"    {g}: {'OK' if g in filtered_genes else 'MISSING'}")

tcga_genes = set()
with open(TCGA_EXPR, "r") as f:
    f.readline()
    for line in f:
        tcga_genes.add(line.split("\t")[0].strip())
overlap = set(filtered_genes) & tcga_genes
print(f"    Overlap with TCGA: {len(overlap):,} / {len(filtered_genes):,}")

# ── 5. Export scRNA-seq reference ────────────────────────────
print("\n[5] Exporting lean scRNA-seq reference...")
ref_path = os.path.join(OUT_DIR, "cibersortx_sc_reference_v2.txt")
with open(ref_path, "w") as f:
    f.write("Gene\t" + "\t".join(all_barcodes) + "\n")
    for gene in filtered_genes:
        row_vals = combined_dense[:, gene_names.index(gene)]
        f.write(gene + "\t" + "\t".join(str(v) for v in row_vals) + "\n")
file_size_mb = os.path.getsize(ref_path) / 1e6
print(f"    Saved: cibersortx_sc_reference_v2.txt ({file_size_mb:.1f} MB)")

# ── 6. Phenotype file ────────────────────────────────────────
print("\n[6] Exporting phenotype file...")
pheno_path = os.path.join(OUT_DIR, "cibersortx_phenotypes_v2.txt")
with open(pheno_path, "w") as f:
    for bc, group in zip(all_barcodes, all_groups):
        f.write(f"{bc}\t{group}\n")
print("    Saved: cibersortx_phenotypes_v2.txt")

# ── 7. Pre-built signature matrix ────────────────────────────
print("\n[7] Building pre-built signature matrix...")
sig_rows = []
for gene in filtered_genes:
    gene_idx = gene_names.index(gene)
    row = {"Gene": gene}
    for group in unique_groups:
        mask = groups_array == group
        group_raw = combined_dense[mask, gene_idx].astype(float)
        group_total = combined_dense[mask, :].sum(axis=1).astype(float)
        group_total[group_total == 0] = 1
        row[group] = round((group_raw / group_total * 1e6).mean(), 4)
    sig_rows.append(row)
sig_df = pd.DataFrame(sig_rows)
sig_path = os.path.join(OUT_DIR, "cibersortx_sig_matrix.txt")
sig_df.to_csv(sig_path, sep="\t", index=False)
sig_size_mb = os.path.getsize(sig_path) / 1e6
print(f"    Saved: cibersortx_sig_matrix.txt ({sig_size_mb:.1f} MB; "
      f"{len(filtered_genes)} genes x {len(unique_groups)} cell types)")

# ── 8. Ensure mixture file ───────────────────────────────────
mixture_path = os.path.join(OUT_DIR, "cibersortx_mixture.txt")
if not os.path.exists(mixture_path):
    os.symlink(os.path.abspath(TCGA_EXPR), mixture_path)
    print("\n[8] Symlinked mixture file")
else:
    print("\n[8] Mixture file already exists")

print(f"\n{'=' * 65}")
print("  Step complete!")
print(f"  Output: {OUT_DIR}")
print(f"{'=' * 65}")
