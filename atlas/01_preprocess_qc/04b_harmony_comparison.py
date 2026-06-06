#!/usr/bin/env python
"""
Atlas 01 — Step 04b: Harmony integration (COMPARISON / OPTIONAL)

PURPOSE
    Run Harmony batch correction (batch_key=sample_id) on the PCA embedding as an
    alternative integration for benchmarking against scVI/scANVI. This is a
    COMPARISON method — it is NOT on the canonical path that produces the atlas
    used downstream (that is scVI → CellAssign → scANVI, steps 04→05→06→07).

INPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_preprocess.h5ad   (from step 03)
        Requires obsm["X_pca"] to be present (see NOTE).

OUTPUTS
    output_root/01_preprocess_qc/integration/ovca_atlas_harmony.h5ad   (obsm["Harmony"])

RUNTIME TIER
    moderate-heavy (Harmony on full-atlas PCA).

MANUSCRIPT ROLE
    Integration-method comparison only (e.g. SF integration benchmark). Optional.

NOTE
    AUTHORITATIVE original cluster script (atlas_02_harmony.py). Logic preserved
    EXACTLY; only the cluster paths were centralised and the random seed sourced
    from config SEED (was a hardcoded 13).
    DEPENDENCY GAP: this script reads obsm["X_pca"], but the authoritative
    preprocess step (03_preprocess_hvg / atlas_01_preprocess.py) does NOT compute
    PCA — the monolithic atlas_integration.py ran sc.tl.pca inline before Harmony.
    To run this comparison you must first add a PCA pass to the preprocess object
    (sc.tl.pca on the HVGs). FLAGGED — see this stage's README.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from config.config import path, SEED  # noqa: E402

import scanpy as sc, anndata as ad, gc
from harmony import harmonize

IN  = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_preprocess.h5ad")
OUT = path("output_root", "01_preprocess_qc", "integration", "ovca_atlas_harmony.h5ad")

print("Loading pre-processed atlas …", flush=True)
adata = sc.read_h5ad(IN)

print("Running Harmony batch correction (sample_id) …", flush=True)
adata.obsm["Harmony"] = harmonize(
    adata.obsm["X_pca"],  # already stored by the preprocess step
    adata.obs,
    batch_key="sample_id",
    random_state=SEED,  # random seed (was hardcoded 13)
)

adata.write_h5ad(OUT, compression="gzip")
del adata
gc.collect()
print("Harmony embedding written to", OUT, flush=True)
