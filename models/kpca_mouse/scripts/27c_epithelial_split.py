#!/usr/bin/env python3
# Recovered copy (2026-10-22). This code ran inline (a script in the session's temporary folder, split.py) from the
# mouse project folder 202608_kpca_secb_mouse/ on 2026-09-09 19:16 UTC and was never saved as a file; it is
# recovered from the analysis session log. Paths are relative to the mouse project folder, so run it from there.
# Writes epi_split.h5ad, which 33_consolidate_handoff.py reads (its SCRATCH path): epithelial Leiden clusters
# (epi_leiden, resolution 0.5) and SecA / SecB UCell from 5 + 4 genes (Python UCell, maxRank 1500) with polarization
# = SecB - SecA. Removed: the split-definition figure (24_E_split_definition, not in the paper) and the matplotlib
# import it needed; the rest is unchanged.
# The run printed: the epithelial Leiden cluster table (cluster 5: 2,769 cells, polarization 0.440).

import scanpy as sc, numpy as np, pandas as pd
E=sc.read_h5ad('/private/tmp/claude-501/-Volumes-CookLab-Sarah-202601-final-atlas/1b3fdb57-aec0-43ae-8eba-c2e8c7bb77b9/scratchpad/epi_umap.h5ad')
xy=E.obsm['X_umap_epi']
# recompute graph on epi HVG for leiden
Ee=E.copy(); sc.pp.highly_variable_genes(Ee,n_top_genes=3000); Ee=Ee[:,Ee.var['highly_variable']].copy()
sc.pp.scale(Ee,max_value=10); sc.tl.pca(Ee,n_comps=30,random_state=0); sc.pp.neighbors(Ee,n_neighbors=15,n_pcs=30,random_state=0)
sc.tl.leiden(Ee,resolution=0.5,random_state=0,flavor='igraph',n_iterations=2,directed=False)
E.obs['epi_leiden']=Ee.obs['leiden'].values
SECA=["Bcam","Mecom","Pbx1","Sox17","Rcn2"]; SECB=["Krt19","Krt7","Lcn2","Prss22"]
def ucell(ad,genes,maxrank=1500):
    gi=[ad.var_names.get_loc(g) for g in genes]; n=len(genes); X=ad.X; s=np.empty(ad.n_obs)
    for st in range(0,ad.n_obs,2000):
        en=min(st+2000,ad.n_obs); b=X[st:en]; b=b.toarray() if hasattr(b,'toarray') else np.asarray(b)
        o=np.argsort(-b,axis=1,kind='stable'); r=np.empty_like(o); rows=np.arange(b.shape[0])[:,None]; r[rows,o]=np.arange(1,b.shape[1]+1)[None,:]
        sig=r[:,gi].astype(float); sig[sig>maxrank]=maxrank+1; U=sig.sum(1); s[st:en]=1-(U-n*(n+1)/2)/(n*maxrank)
    return np.clip(s,0,1)
E.obs['SecA']=ucell(E,SECA); E.obs['SecB']=ucell(E,SECB); E.obs['pol']=E.obs['SecB']-E.obs['SecA']
# per-cluster polarization
cl=E.obs.groupby('epi_leiden',observed=True).agg(n=('pol','size'),pol=('pol','mean'),SecB=('SecB','mean'),SecA=('SecA','mean')).sort_values('pol',ascending=False)
print(cl.round(3).to_string())
E.write('/private/tmp/claude-501/-Volumes-CookLab-Sarah-202601-final-atlas/1b3fdb57-aec0-43ae-8eba-c2e8c7bb77b9/scratchpad/epi_split.h5ad')
