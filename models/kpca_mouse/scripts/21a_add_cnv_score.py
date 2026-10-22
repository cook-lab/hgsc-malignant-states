#!/usr/bin/env python3
# Recovered copy (2026-10-22). This code ran inline (a `python3 -` heredoc) from the mouse project folder
# 202608_kpca_secb_mouse/ on 2026-09-04 23:15 UTC and was never saved as a file; it is recovered from the analysis
# session log. Paths are relative to the mouse project folder, so run it from there.
# Adds the per-cell CNV score of step 19 (cnv_score) to the step-15b object, in place; step 21 carries it into
# kpca_working.h5ad. It also stored Leiden clusterings at resolutions 1.5 and 2.0 (l1.5, l2.0), which step 21
# deletes. Removed: a diagnostic table of how the old clusters 6 and 15 split at resolution 2.0 (Krt / Cd3e / Nkg7
# co-expression); the rest is unchanged.
# The run printed: the Leiden cluster counts at resolutions 1.5 and 2.0.

import scanpy as sc, numpy as np, pandas as pd
A=sc.read_h5ad('output/15_doublet_cleanup/kpca_all_clean.h5ad'); A.obs['leiden']=A.obs['leiden'].astype(str)
print("neighbors present:", 'connectivities' in A.obsp)
# push resolution
for res in [1.5,2.0]:
    sc.tl.leiden(A,resolution=res,key_added=f'l{res}',flavor='igraph',n_iterations=2,directed=False,random_state=0)
    print(f"res {res}: {A.obs[f'l{res}'].nunique()} clusters")
pc=pd.read_csv('output/19_identity_confidence/19_percell.csv',index_col=0); A.obs=A.obs.join(pc[['cnv_score']])
A.write_h5ad('output/15_doublet_cleanup/kpca_all_clean.h5ad')  # save with new leiden keys
