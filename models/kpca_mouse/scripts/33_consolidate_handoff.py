#!/usr/bin/env python3
"""Step 33 — consolidate ALL final labels into one handoff object.
Merges: celltype (step 27) + epithelial subclusters/SecB pole/SecB-hi + UCell
scores + de-novo NMF program usage (step 31) into one h5ad for Seurat export."""
import os, sys
from importlib import import_module
import numpy as np, pandas as pd, scanpy as sc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
S = import_module("00_setup")
OUT = os.path.join(S.OUTPUT, "33_handoff"); os.makedirs(OUT, exist_ok=True)
SCRATCH = "/private/tmp/claude-501/-Volumes-CookLab-Sarah-202601-final-atlas/1b3fdb57-aec0-43ae-8eba-c2e8c7bb77b9/scratchpad"

A = sc.read_h5ad(os.path.join(S.OUTPUT, "26_final_annotation", "kpca_annotated.h5ad"))
E = sc.read_h5ad(os.path.join(SCRATCH, "epi_split.h5ad"))     # tumour cells, epi_leiden, SecB/SecA/pol, X_umap_epi
usage = pd.read_csv(os.path.join(S.OUTPUT, "31_fulltumour_nmf", "31_program_usage.csv"), index_col=0)

# epithelial-subset columns, aligned to full object (NA for non-tumour)
A.obs["epi_subcluster"] = "non-epithelial"
A.obs.loc[E.obs_names, "epi_subcluster"] = E.obs["epi_leiden"].astype(str).values
# SecB pole (epi cl5), SecB-high (top25% SecB UCell), hypoxic-high (top25%)
secb_pole = pd.Series(False, index=A.obs_names); secb_pole.loc[E.obs_names] = (E.obs["epi_leiden"] == "5").values
A.obs["SecB_pole"] = secb_pole.values
sb = E.obs["SecB"].values; thr = np.percentile(sb, 75)
hi = pd.Series(False, index=A.obs_names); hi.loc[E.obs_names] = sb >= thr
A.obs["SecB_high_top25"] = hi.values
# refined per-cell scores from the epithelial UCell (adjusted gene sets)
for col, src in [("SecB_UCell_epi", "SecB"), ("SecA_UCell_epi", "SecA"), ("polarization_epi", "pol")]:
    v = pd.Series(np.nan, index=A.obs_names); v.loc[E.obs_names] = E.obs[src].values; A.obs[col] = v.values
# NMF program usage (mP1-4)
for p in usage.columns:
    v = pd.Series(np.nan, index=A.obs_names); v.loc[usage.index] = usage[p].values; A.obs[f"NMF_{p}"] = v.values
# epithelial UMAP coords as obs columns (NA for non-tumour)
for i, nm in enumerate(["epi_umap_1", "epi_umap_2"]):
    v = pd.Series(np.nan, index=A.obs_names); v.loc[E.obs_names] = E.obsm["X_umap_epi"][:, i]; A.obs[nm] = v.values
# a single readable identity column
A.obs["final_identity"] = A.obs["celltype"].astype(str)
A.obs.loc[A.obs["SecB_pole"].values, "final_identity"] = "Tumour_SecB_pole"
A.obs.loc[(A.obs["is_tumour"].values) & (~A.obs["SecB_pole"].values), "final_identity"] = "Tumour_EMT/other"

print("final_identity:\n", A.obs["final_identity"].value_counts().to_string())
print("\nepi_subcluster:\n", A.obs["epi_subcluster"].value_counts().to_string())
A.write(os.path.join(OUT, "kpca_final_annotated.h5ad"))
A.obs.to_csv(os.path.join(OUT, "kpca_final_metadata.csv"))
print(f"\nsaved {OUT}/kpca_final_annotated.h5ad  ({A.n_obs} cells x {A.n_vars} genes)")
