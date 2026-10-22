#!/usr/bin/env python3
# Recovered copy (2026-10-22). This code ran inline (a `python3 -` heredoc) from the mouse project folder
# 202608_kpca_secb_mouse/ on 2026-09-05 00:04 UTC and was never saved as a file; it is recovered from the analysis
# session log. Paths are relative to the mouse project folder, so run it from there.
# Writes output/22_secb_axis/: seca_genes.csv and secb_genes.csv (the mouse SecA / SecB genes that step 24 reads),
# kpca_tumour.h5ad and the tumour matrix that 22b_secb_axis_ucell.R scores. The shell command created the output
# folder first (mkdir -p output/22_secb_axis); the os.makedirs line does that here. The rest is unchanged.
# The run printed: "tumour compartment (clusters 8,9,11): 7,752 cells", "exported for UCell".

import os; os.makedirs('output/22_secb_axis', exist_ok=True)
import scanpy as sc, numpy as np, pandas as pd
import sys; sys.path.insert(0,'scripts'); from importlib import import_module
S=import_module("00_setup")
B=sc.read_h5ad('output/21_working/kpca_working.h5ad'); B.obs['leiden']=B.obs['leiden'].astype(str)
T=B[B.obs.leiden.isin(['8','9','11'])].copy()
print(f"tumour compartment (clusters 8,9,11): {T.shape[0]:,} cells")
# re-embed tumour locally for a clean tumour-only UMAP
sc.pp.highly_variable_genes(T,n_top_genes=2000)
sc.pp.pca(T,n_comps=30,svd_solver='arpack',random_state=S.SEED)
sc.pp.neighbors(T,n_neighbors=15,n_pcs=20,random_state=S.SEED)
sc.tl.umap(T,random_state=S.SEED); T.obsm['X_umap_tum']=T.obsm['X_umap'].copy()
sc.tl.leiden(T,resolution=0.6,key_added='tum_leiden',flavor='igraph',n_iterations=2,directed=False,random_state=S.SEED)
T.write_h5ad('output/22_secb_axis/kpca_tumour.h5ad')
# export raw counts for UCell (genes x cells) — use human sets mapped to mouse
ortho=S.load_ortholog_map()
seca_h=['BCAM','MECOM','FBXO21','LGR5','LPAR3','PBX1','SOX17','RCN2']
secb_h=['KRT17','KRT19','KRT7','LCN2','PRSS22','SLPI','TACSTD2']
sA,_=S.map_genes(seca_h,ortho,restrict_to=set(T.var_names)); sB,_=S.map_genes(secb_h,ortho,restrict_to=set(T.var_names))
print("SecA mouse:",list(sA.values())); print("SecB mouse:",list(sB.values()))
pd.Series(list(sA.values())).to_csv('output/22_secb_axis/seca_genes.csv',index=False,header=['gene'])
pd.Series(list(sB.values())).to_csv('output/22_secb_axis/secb_genes.csv',index=False,header=['gene'])
# lognorm matrix for UCell (UCell wants counts or lognorm ranks; use lognorm X)
import scipy.io, scipy.sparse as sp
X=T.X.tocsr() if sp.issparse(T.X) else sp.csr_matrix(T.X)
scipy.io.mmwrite('output/22_secb_axis/tumour_X.mtx', X.T)  # genes x cells
pd.Series(T.var_names).to_csv('output/22_secb_axis/genes.csv',index=False,header=['gene'])
pd.Series(T.obs_names).to_csv('output/22_secb_axis/cells.csv',index=False,header=['cell'])
print("exported for UCell")
