#!/usr/bin/env python3
# Visium HD (16 um), WT01 + WT03 (descriptive n = 2): candidate genes, both partners.
#   mac side  : niche mac bins (>=2 SecB-hi within 3 bins) vs distal (none); log2FC of mean CP10k, adjusted for macrophage
#               content (minus median log2FC of pan-mac markers). Carry-over rule: a gene is attributable to macrophages only
#               if its mean in niche mac bins >= its mean in pure SecB-hi epithelial bins (else "carry-over", uninformative).
#   SecB side : SecB-hi vs SecA-hi pure epithelial bins; log2FC of mean CP10k.
# Bin-level Mann-Whitney p (descriptive only).
import numpy as np, pandas as pd, scanpy as sc
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
import os
L=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler/visium_local"); O=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler")
SAMPLES={"WT01B1":"WT01","WT03B1":"WT03"}
SECB=["KRT17","KRT19","KRT7","LCN2","PRSS22","SLPI","TACSTD2"]; SECA=["MECOM","FBXO21","LGR5","LPAR3","PBX1","SOX17","RCN2"]
EPI=["EPCAM","KRT8","KRT18","PAX8","WFDC2","MUC16"]; MAC=["CD68","CD163","CSF1R","AIF1","C1QA","C1QB","C1QC","FCGR3A","CD14","LYZ"]
PANMAC=["CD68","CSF1R","AIF1","CD14","FCGR3A","TYROBP","LAPTM5"]
CAND=["TNF","NFKB1","NFKBIA","BIRC3","BCL2L1","MCL1","VEGFA","ADM","EREG","CD55","CD47","LGALS9","LGALS3","CD274","C7","IL32","TGFB1",
      "INHBA","SLC16A3","LDHA","SLC2A1","NDRG1","PGK1","ENO1","ALDOA","PDK1","BNIP3","HK2","CA9","ICAM1","VCAM1","CD44","SPP1","CTSL",
      "CTSD","MMP7","MMP9","TIMP1","MIF","MDK","SAA1","GPNMB","MARCO","APOE","TREM2","CXCL2","CCL20","C1QC","MRC1","CD163","STAB1","A2M",
      "CD83","MXD1","TGM2","ITGB6","ITGAV","HAVCR2","LAG3"]
rows=[]
for s,pid in SAMPLES.items():
    a=sc.read_10x_h5(f"{L}/{s}/filtered_feature_bc_matrix.h5"); a.var_names_make_unique()
    pos=pd.read_parquet(f"{L}/{s}/tissue_positions.parquet").set_index("barcode"); a=a[a.obs_names.isin(pos.index)].copy()
    sc.pp.filter_cells(a,min_counts=20); pos=pos.loc[a.obs_names]; raw=a.X.tocsc(); rv={g:i for i,g in enumerate(a.var_names)}
    umi=lambda gl: np.asarray(raw[:,[rv[g] for g in gl if g in rv]].sum(1)).ravel(); epiU,macU=umi(EPI),umi(MAC)
    lib=np.asarray(a.X.sum(1)).ravel(); CP=a.X.multiply(1e4/lib[:,None]).tocsc(); LN=CP.copy(); LN.data=np.log1p(LN.data)
    def z(gl): M=LN[:,[rv[g] for g in gl if g in rv]].toarray(); return ((M-M.mean(0))/(M.std(0)+1e-9)).mean(1)
    pol=z(SECB)-z(SECA); rc=pos[["array_row","array_col"]].values.astype(float)
    epi=(epiU>=5)&(macU==0); mac=(macU>=3)&(macU>=epiU); pe=pol[epi]; secb=epi&(pol>=np.quantile(pe,2/3)); seca=epi&(pol<=np.quantile(pe,1/3))
    mi=np.where(mac)[0]; nB=np.array([secb[n].sum() for n in cKDTree(rc).query_ball_point(rc[mi],r=3.0,p=np.inf)])
    ni,di,bi,ai=mi[nB>=2],mi[nB==0],np.where(secb)[0],np.where(seca)[0]
    mean=lambda idx,g: CP[idx,rv[g]].mean()
    l2=lambda x,y: np.log2(x+0.01)-np.log2(y+0.01)
    adj=np.median([l2(mean(ni,g),mean(di,g)) for g in PANMAC if g in rv])
    for g in CAND:
        if g not in rv: continue
        col=LN[:,rv[g]].toarray().ravel()
        mn,md,mb,ma=mean(ni,g),mean(di,g),mean(bi,g),mean(ai,g)
        rows.append(dict(sample=pid,gene=g,side="mac",effect=l2(mn,md)-adj,mac_sourced=mn>=mb,det=(col[ni]>0).mean(),
                         p=mannwhitneyu(col[ni],col[di]).pvalue))
        rows.append(dict(sample=pid,gene=g,side="secb",effect=l2(mb,ma),mac_sourced=np.nan,det=(col[bi]>0).mean(),
                         p=mannwhitneyu(col[bi],col[ai]).pvalue))
    print(pid,"panmac adj %.2f"%adj,flush=True)
pd.DataFrame(rows).to_csv(f"{O}/visium_candidates_2sections.csv",index=False); print("VCAND_DONE")
