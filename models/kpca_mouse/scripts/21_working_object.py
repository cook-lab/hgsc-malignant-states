#!/usr/bin/env python3
# Recovered copy (2026-10-22). This code ran inline (a `python3 -` heredoc) from the mouse project folder
# 202608_kpca_secb_mouse/ on 2026-09-04 23:19 UTC and was never saved as a file; it is recovered from the analysis
# session log. Paths are relative to the mouse project folder, so run it from there.
# Writes output/21_working/kpca_working.h5ad (on disk; its file time matches this run). The code below is unchanged.
# The run printed: "excluding old 6/15: 22,437 -> 20,634 cells", "re-embedded: 16 clusters".

import scanpy as sc, numpy as np, pandas as pd
import sys; sys.path.insert(0,'scripts'); from importlib import import_module
S=import_module("00_setup")
A=sc.read_h5ad('output/15_doublet_cleanup/kpca_all_clean.h5ad'); A.obs['leiden']=A.obs['leiden'].astype(str)
B=A[~A.obs.leiden.isin(['6','15'])].copy()
print(f"excluding old 6/15: {A.shape[0]:,} -> {B.shape[0]:,} cells")
# clean transient cols
for c in list(B.obs.columns):
    if c.startswith(('l1.5','l2.0','_')): del B.obs[c]
sc.pp.highly_variable_genes(B,n_top_genes=S.N_HVG)
sc.pp.pca(B,n_comps=50,svd_solver='arpack',random_state=S.SEED)
sc.pp.neighbors(B,n_neighbors=15,n_pcs=40,random_state=S.SEED)
sc.tl.umap(B,random_state=S.SEED)
sc.tl.leiden(B,resolution=1.0,key_added='leiden',flavor='igraph',n_iterations=2,directed=False,random_state=S.SEED)
print(f"re-embedded: {B.obs.leiden.nunique()} clusters")
import os; os.makedirs('output/21_working',exist_ok=True)
B.write_h5ad('output/21_working/kpca_working.h5ad')
print("wrote output/21_working/kpca_working.h5ad")
# check SecA/SecB gene availability (7 each, mouse orthologs)
seca=['Bcam','Mecom','Sox17','Pbx1','Lgr5','Rcn2','Wt1']
secb=['Krt7','Krt19','Krt17','Lcn2','Mmp7','Slpi','Tacstd2']
def det(g): 
    v=B[:,g].X; v=np.asarray(v.todense()).ravel() if hasattr(v,'todense') else np.asarray(v).ravel(); return 100*(v>0).mean()
print("\nSecA genes:",{g:f'{det(g):.0f}%' for g in seca if g in B.var_names})
print("SecB genes:",{g:f'{det(g):.0f}%' for g in secb if g in B.var_names})
print("missing:",[g for g in seca+secb if g not in B.var_names])
