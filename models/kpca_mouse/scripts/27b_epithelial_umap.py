#!/usr/bin/env python3
# Recovered copy (2026-10-22). This code ran inline (a script in the session's temporary folder, epi_secgenes.py)
# from the mouse project folder 202608_kpca_secb_mouse/ on 2026-09-09 19:05 UTC and was never saved as a file; it is
# recovered from the analysis session log. Paths are relative to the mouse project folder, so run it from there.
# Writes epi_umap.h5ad: the 12,429 tumour cells of kpca_annotated.h5ad (step 27) with the epithelial UMAP
# (X_umap_epi) of Fig 3j, 3k, S10b and S10d. The file went to that temporary folder (the path below), which is gone;
# step 27c reads it there. Removed: the figure of SecA / SecB genes on the epithelial UMAP (24_C_epi_secAB_genes,
# not in the paper) and the matplotlib import it needed; the rest is unchanged.
# The run printed: "epithelial/tumour cells: 12429".

import scanpy as sc, numpy as np, pandas as pd
import importlib.util
spec=importlib.util.spec_from_file_location("s","scripts/00_setup.py"); S=importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
A=sc.read_h5ad('output/26_final_annotation/kpca_annotated.h5ad')
E=A[A.obs['is_tumour']].copy()
print('epithelial/tumour cells:',E.n_obs)
# fresh tumour-only UMAP
Ee=E.copy(); sc.pp.highly_variable_genes(Ee,n_top_genes=3000); Ee=Ee[:,Ee.var['highly_variable']].copy()
sc.pp.scale(Ee,max_value=10); sc.tl.pca(Ee,n_comps=30,random_state=0)
sc.pp.neighbors(Ee,n_neighbors=15,n_pcs=30,random_state=0); sc.tl.umap(Ee,random_state=0)
xy=Ee.obsm['X_umap']; E.obsm['X_umap_epi']=xy
E.write('/private/tmp/claude-501/-Volumes-CookLab-Sarah-202601-final-atlas/1b3fdb57-aec0-43ae-8eba-c2e8c7bb77b9/scratchpad/epi_umap.h5ad')
