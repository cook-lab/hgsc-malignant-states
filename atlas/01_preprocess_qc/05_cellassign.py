#!/usr/bin/env python
"""
Atlas 01 — Step 05: CellAssign cell-type priors

PURPOSE
    Run CellAssign on the full-gene-space raw atlas using a curated binary
    marker matrix, to produce hard cell-type predictions (celltype_pred) that
    seed the scANVI semi-supervised refinement (step 06). Library-size factors
    are computed on the full count matrix; only marker-gene columns are passed
    to CellAssign. Labels are written back onto the raw object.

INPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_raw.h5ad   (from step 02)
    shared/cellassign_markers.csv   (binary marker matrix; rows=genes, cols=cell types)

OUTPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_raw.h5ad   (in place: + obs["celltype_pred"])

RUNTIME TIER
    heavy — GPU. CellAssign training.

MANUSCRIPT ROLE
    Cell-type priors for scANVI label propagation. The final cell-type labels
    used in the manuscript derive from these priors refined by scANVI (step 06).

NOTE
    AUTHORITATIVE original cluster script (atlas_04a_cellassign.py). Analytical
    logic preserved EXACTLY (size factors on full counts, marker-column subset,
    signatures = (marker_mat > 0), idxmax hard labels). Only the cluster paths
    were centralised; the marker matrix is the original cellassign_markers.csv
    copied into shared/.
    MARKER-MATRIX FLAG: this original matrix is the 53-gene × 11-type version
    (Epithelial/Mesothelial/Fibroblasts/Endothelial/T_NK/B/Macrophage/DC/Plasma/
    Mast/Other). The repo also carries shared/cellassign_markers_v3.csv (81-gene ×
    16-type) used by the later 20260213 production integration that produced the
    deposited hgsc_atlas_scanvi.h5ad. The two are NOT interchangeable — see README.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path  # noqa: E402

# =========================
# PATHS (central config)
# =========================
IN_RAW      = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_raw.h5ad")   # full gene-space AnnData (counts in .X)
MARKERS     = str(Path(__file__).resolve().parents[2] / "shared" / "cellassign_markers.csv")  # wide binary: rows=genes (index), cols=cell types
OUT_RAW     = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_raw.h5ad")
LABEL_COL   = "celltype_pred"
# =========================

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from scvi.external import CellAssign

print(f"🔹  Load RAW: {IN_RAW}", flush=True)
adata = sc.read_h5ad(IN_RAW)

# Keep counts sparse; avoid extra copies
if not sp.isspmatrix_csr(adata.X):
    adata.X = sp.csr_matrix(adata.X)


# ----------  CELLASSIGN  ----------
print("🔹  Prepare full raw-count matrix for size-factor calculation", flush=True)
full_counts = adata.X            # CSR sparse matrix
# Use var_names directly to avoid creating a huge .raw copy
genes       = adata.var_names.to_numpy()

# Compute true library sizes on the full matrix
lib_size     = np.ravel(full_counts.sum(axis=1))
size_factors = lib_size / lib_size.mean()

print("🔹  Load marker‐gene definitions and intersect with genes", flush=True)
# Expect a wide binary matrix with genes as index; keep it simple per your original
marker_gene_mat = pd.read_csv(MARKERS, index_col=0)
have = marker_gene_mat.index.intersection(genes)
if have.empty:
    raise ValueError("No marker genes found in the dataset!")

# Compute the integer column indices of those markers
idx = np.where(np.isin(genes, have))[0]
idx.sort()    # ensure sorted → sparse slice stays sparse

print(f"🔹  Subsetting {full_counts.shape[1]:,} → {idx.size:,} marker columns", flush=True)
marker_counts = full_counts[:, idx]              # still sparse

# Build the lean AnnData for CellAssign
bdata = ad.AnnData(
    X   = marker_counts,
    obs = adata.obs.copy(),
    var = pd.DataFrame(index=genes[idx])
)
bdata.obs["size_factor"] = size_factors

# Build signatures aligned to the subset var index (columns = cell types)
signatures = (marker_gene_mat.loc[bdata.var_names] > 0).astype(int)

# Minimal CellAssign calls; pass the size_factor so your computation is used
CellAssign.setup_anndata(bdata, size_factor_key="size_factor")
ca = CellAssign(bdata, signatures)
ca.train()
# After training:
pred = ca.predict()  # cells x celltypes probabilities (but index may be 0..n-1)

# 🔧 Ensure the prediction index matches your cells:
pred.index = bdata.obs_names    # <-- add this line

# Hard labels back onto RAW
adata.obs[LABEL_COL] = (
    pred.idxmax(axis=1)
        .reindex(adata.obs_names)
        .astype("category")
)

print("Non-NA assigned:", adata.obs[LABEL_COL].notna().sum(), "of", adata.n_obs, flush=True)
print(adata.obs[LABEL_COL].value_counts(dropna=False).head(20), flush=True)

print(f"🔹  Write RAW with CellAssign labels → {OUT_RAW}", flush=True)
adata.write_h5ad(OUT_RAW, compression="gzip")
print("✅  CellAssign finished.", flush=True)
