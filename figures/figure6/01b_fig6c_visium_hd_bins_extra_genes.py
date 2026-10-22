#!/usr/bin/env python3
# PURPOSE: Fig 6c -- adds log1p CP10k of ITGB6, ITGAV, CA9, CTSL and TIMP1 (16 um bins of WT01B1) to the WT01 bin table
#   of 01_fig6abcd_visium_hd_bins.py, in place. Run it after 01. Panel 6c reads the CA9 column; no panel reads the
#   other four.
# Recovered (2026-10-22): this ran inline in the final run's fig6_singler output folder on 2026-09-25 16:33 UTC and was
#   never saved as a file; it is recovered from the analysis session log. Changed: the folder and file names, which
#   now go through OUTPUT_ROOT and the study codes; the code is otherwise as it ran.
# OUTPUTS: <fig6_singler>/visium_n2_bins_WT01.csv.gz, rewritten with the five columns (read by 13_fig6_draw_panels.py).
# MANUSCRIPT PANEL(S): Fig 6c. RUNTIME TIER: fast.
import os
O=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler")
import numpy as np, pandas as pd, scanpy as sc
L=f"{O}/visium_local/WT01B1"; B=pd.read_csv(f"{O}/visium_n2_bins_WT01.csv.gz")
a=sc.read_10x_h5(f"{L}/filtered_feature_bc_matrix.h5"); a.var_names_make_unique(); a=a[B.barcode.values].copy()
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a); Xc=a.X.tocsc(); rv={g:i for i,g in enumerate(a.var_names)}
for g in ["ITGB6","ITGAV","CA9","CTSL","TIMP1"]:
    if g in rv: B[g]=Xc[:,rv[g]].toarray().ravel()
B.to_csv(f"{O}/visium_n2_bins_WT01.csv.gz",index=False); print(B.columns.tolist())
