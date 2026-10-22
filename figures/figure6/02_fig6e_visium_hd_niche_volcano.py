#!/usr/bin/env python3
# Visium HD WT01: genome-wide niche mac bins (>=2 SecB-hi within 3 bins) vs distal mac bins (none). Descriptive (n=1).
# log2FC of mean CP10k (+0.01), bin-level Mann-Whitney + BH (bins not independent -> descriptive only).
import sys, numpy as np, pandas as pd, scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
BIN=sys.argv[1] if len(sys.argv)>1 else "16"   # "16" or "8" um bins
import os
L=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler/visium_local/WT01B1")+("_8um" if BIN=="8" else ""); O=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler")
SECB=["KRT17","KRT19","KRT7","LCN2","PRSS22","SLPI","TACSTD2"]; SECA=["MECOM","FBXO21","LGR5","LPAR3","PBX1","SOX17","RCN2"]
EPI=["EPCAM","KRT8","KRT18","PAX8","WFDC2","MUC16"]; MAC=["CD68","CD163","CSF1R","AIF1","C1QA","C1QB","C1QC","FCGR3A","CD14","LYZ"]
a=sc.read_10x_h5(f"{L}/filtered_feature_bc_matrix.h5"); a.var_names_make_unique()
pos=pd.read_parquet(f"{L}/tissue_positions.parquet" if BIN=="16" else f"{L}/spatial/tissue_positions.parquet").set_index("barcode")
a=a[a.obs_names.isin(pos.index)].copy(); sc.pp.filter_cells(a,min_counts=20 if BIN=="16" else 10); pos=pos.loc[a.obs_names]
raw=a.X.tocsc(); rv={g:i for i,g in enumerate(a.var_names)}
umi=lambda gl: np.asarray(raw[:,[rv[g] for g in gl if g in rv]].sum(1)).ravel()
epiU,macU=umi(EPI),umi(MAC)
lib=np.asarray(a.X.sum(1)).ravel(); CP=a.X.multiply(1e4/lib[:,None]).tocsr(); LN=CP.copy(); LN.data=np.log1p(LN.data); LNc=LN.tocsc()
def z(gl): M=LNc[:,[rv[g] for g in gl if g in rv]].toarray(); return ((M-M.mean(0))/(M.std(0)+1e-9)).mean(1)
pol=z(SECB)-z(SECA); rc=pos[["array_row","array_col"]].values.astype(float)
EMIN,MMIN,RAD=(5,3,3.0) if BIN=="16" else (3,2,6.0)   # 8 um: same ~48 um radius
epi=(epiU>=EMIN)&(macU==0); mac=(macU>=MMIN)&(macU>=epiU); pe=pol[epi]; secb=epi&(pol>=np.quantile(pe,2/3))
tree=cKDTree(rc); mi=np.where(mac)[0]; nB=np.array([secb[n].sum() for n in tree.query_ball_point(rc[mi],r=RAD,p=np.inf)])
ni,di=mi[nB>=2],mi[nB==0]
det=np.maximum(np.asarray((CP[ni]>0).mean(0)).ravel(),np.asarray((CP[di]>0).mean(0)).ravel()); keep=np.where(det>=0.05)[0]
mA=np.asarray(CP[ni].mean(0)).ravel(); mB=np.asarray(CP[di].mean(0)).ravel()
A=LNc[ni][:,keep].tocsc(); B=LNc[di][:,keep].tocsc()
P=[mannwhitneyu(A[:,j].toarray().ravel(),B[:,j].toarray().ravel()).pvalue for j in range(len(keep))]
D=pd.DataFrame({"gene":a.var_names[keep],"log2FC":np.log2(mA[keep]+0.01)-np.log2(mB[keep]+0.01),"p_bin":P,
                "det_niche":np.asarray((CP[ni]>0).mean(0)).ravel()[keep],"det_distal":np.asarray((CP[di]>0).mean(0)).ravel()[keep]})
sb=np.where(secb)[0]; D["mean_secb_epi"]=np.asarray(CP[sb].mean(0)).ravel()[keep]; D["mean_niche"]=mA[keep]
# carry-over rule: a transcript carried into a mixed bin from SecB epithelium cannot exceed its level in pure SecB-hi bins
D["mac_sourced"]=D.mean_niche>=D.mean_secb_epi
D["fdr_bin"]=multipletests(D.p_bin,method="fdr_bh")[1]; D=D.sort_values("log2FC",ascending=False)
D.to_csv(f"{O}/visiumWT01_niche_vs_distal_volcano"+("_8um" if BIN=="8" else "")+".csv",index=False)
print(len(ni),len(di),len(D)); print(D.head(30).round(3).to_string(index=False)); print(D.tail(15).round(3).to_string(index=False))
