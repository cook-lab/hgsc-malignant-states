#!/usr/bin/env python3
# PURPOSE: Fig 6a-d -- Visium HD (16 um bins, whole transcriptome), per-bin classes, polarization and macrophage
#   program. Descriptive. Sections: WT01B1 (= WT01, the Fig 4 ROI tissue; the bins of panels a-d) and WT03B1 (= WT03,
#   a second Visium HD section processed alongside; its bin table is not read by any panel).
# Label-independent of Xenium: bins are classed from raw marker UMIs.
#   pure epithelial bin = >=5 epithelial UMIs and ZERO macrophage UMIs
#   macrophage bin      = >=3 macrophage UMIs and >= epithelial UMIs
#   SecB-hi / SecA-hi   = top / bottom tertile of (SecB - SecA score) among pure epithelial bins
#   mac (SecB niche)    = mac bin with >=2 SecB-hi bins within 3 bins Chebyshev (~48 um); mac (distal) = with none
# OUTPUTS: <fig6_singler>/visium_n2_bins_WT01.csv.gz, visium_n2_bins_WT03.csv.gz (barcode, row, col, cls, pol, mac_prog,
#   log1p CP10k of CD47, SPP1, MIF, CD55, SLC16A3, ICAM1). The WT01 table is read by 13_fig6_draw_panels.py (6a-d),
#   fig6_roi_picker/roi_server.py and tables/10_tableS10_secb_niche_macrophage_degs.py (bin-class counts).
#   Panel 6c also needs a CA9 column, which 01b_fig6c_visium_hd_bins_extra_genes.py adds (run it after this script).
# MANUSCRIPT PANEL(S): Fig 6a, 6b, 6c, 6d. RUNTIME TIER: moderate.
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, spearmanr, pearsonr
import os
L=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler/visium_local"); O=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "fig6_singler")
SAMPLES={"WT01B1":"WT01","WT03B1":"WT03"}
SECB=["KRT17","KRT19","KRT7","LCN2","PRSS22","SLPI","TACSTD2"]; SECA=["MECOM","FBXO21","LGR5","LPAR3","PBX1","SOX17","RCN2"]
EPI=["EPCAM","KRT8","KRT18","PAX8","WFDC2","MUC16"]; MAC=["CD68","CD163","CSF1R","AIF1","C1QA","C1QB","C1QC","FCGR3A","CD14","LYZ"]
PROG_UP=["LDHA","ALDOA","ENO1","PGK1","PDK1","SLC2A1","SLC16A3","NDRG1","VEGFA"]; PROG_DN=["C1QC","MRC1","A2M","C7"]
def zscore(X,var,gl):
    ix=[var[g] for g in gl if g in var]; M=X[:,ix]; M=M.toarray() if sp.issparse(M) else M
    return ((M-M.mean(0))/(M.std(0)+1e-9)).mean(1)
for s,pid in SAMPLES.items():
    a=sc.read_10x_h5(f"{L}/{s}/filtered_feature_bc_matrix.h5"); a.var_names_make_unique()
    pos=pd.read_parquet(f"{L}/{s}/tissue_positions.parquet").set_index("barcode")
    a=a[a.obs_names.isin(pos.index)].copy(); sc.pp.filter_cells(a,min_counts=20); pos=pos.loc[a.obs_names]
    raw=a.X.tocsc(); rvar={g:i for i,g in enumerate(a.var_names)}
    umi=lambda gl: np.asarray(raw[:,[rvar[g] for g in gl if g in rvar]].sum(1)).ravel()
    epiU,macU=umi(EPI),umi(MAC)
    sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
    X=a.X.tocsr(); Xc=X.tocsc(); var={g:i for i,g in enumerate(a.var_names)}
    pol=zscore(Xc,var,SECB)-zscore(Xc,var,SECA)
    rc=pos[["array_row","array_col"]].values.astype(float)
    epi_dom=(epiU>=5)&(macU==0); mac_bin=(macU>=3)&(macU>=epiU)
    pe=pol[epi_dom]; secb_hi=epi_dom&(pol>=np.quantile(pe,2/3)); seca_hi=epi_dom&(pol<=np.quantile(pe,1/3))
    tree=cKDTree(rc)
    def nb_count(q,t,r=3.0):
        tg=np.zeros(len(t),bool); tg[:]=t
        return np.array([tg[nb].sum() for nb in tree.query_ball_point(rc[q],r=r,p=np.inf)],int)
    m_nsecb=nb_count(mac_bin,secb_hi); b_nmac=nb_count(secb_hi,mac_bin)
    mi=np.where(mac_bin)[0]; bi=np.where(secb_hi)[0]
    CMP={"V_mac_niche":(mi[m_nsecb>=2],mi[m_nsecb==0]),
         "V_secb_contact":(bi[b_nmac>=1],bi[b_nmac==0]),
         "V_secb_vs_seca":(bi,np.where(seca_hi)[0])}
    # per-bin table for maps
    cls=np.full(a.n_obs,"other",object); cls[epi_dom]="epithelium"; cls[seca_hi]="SecA-hi"; cls[secb_hi]="SecB-hi"
    cls[mac_bin]="mac (other)"; cls[CMP["V_mac_niche"][1]]="mac (distal)"; cls[CMP["V_mac_niche"][0]]="mac (SecB niche)"
    prog=zscore(Xc,var,PROG_UP)-zscore(Xc,var,PROG_DN)
    B=pd.DataFrame({"barcode":a.obs_names,"row":pos.array_row.values,"col":pos.array_col.values,"cls":cls,"pol":pol,"mac_prog":prog})
    for g in ["CD47","SPP1","MIF","CD55","SLC16A3","ICAM1"]:
        if g in var: B[g]=Xc[:,var[g]].toarray().ravel()
    B.to_csv(f"{O}/visium_n2_bins_{pid}.csv.gz",index=False)
    print(f"{pid}: bins {a.n_obs}  pure-epi {epi_dom.sum()}  SecB-hi {secb_hi.sum()}  SecA-hi {seca_hi.sum()}  mac {mac_bin.sum()}  "
          f"niche/distal {len(CMP['V_mac_niche'][0])}/{len(CMP['V_mac_niche'][1])}  "
          f"SecB contact/none {len(CMP['V_secb_contact'][0])}/{len(CMP['V_secb_contact'][1])}",flush=True)
print("VISIUM_N2_DONE")
