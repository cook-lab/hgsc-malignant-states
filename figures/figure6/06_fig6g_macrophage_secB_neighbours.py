#!/usr/bin/env python3
# PURPOSE: Fig 6g co-localisation (Xenium, SingleR labels). Macrophage program = z(glycolytic/HIF) -
# z(resident/complement), z within each sample's macrophages (whole tissue or TMA array file). Within each unit
# (WT section; TMA patient) macrophages are split into program tertiles: hypoxic-high (top) vs resident-high
# (bottom). Per macrophage:
#   secb_frac50 = SecB / all secretory epithelium within 50 um (>=3 epithelial neighbours)
#   dist_vasc   = distance to nearest endothelial cell or pericyte (not plotted)
# A whole tissue or TMA core enters with >= 30 macrophages, >= 50 secretory and >= 5 vascular cells.
# TMA: the final TMA core set (spatial/02_qc/06b_final_tma_cores.py: per-array core QC on SingleR labels,
# primary HGSC, treatment-naive, >= 500 cells per patient), neighbours searched within core (TMA_1/TMA_2
# coordinates overlap).
# OUTPUTS (<fig6_singler>/): fig6_colocal_percell.csv.gz (per-macrophage record), fig6_colocal_units.csv (one row per
#   whole tissue or TMA patient; drawn and tested by 13_fig6_draw_panels.py); prints the unit medians of the text.
# MANUSCRIPT PANEL(S): Fig 6g ("twice the proportion of SecB neighbors"). RUNTIME TIER: moderate.
import gzip, io, os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from figures._prep.final_tma_cores import final_core_ids
from scipy.io import mmread; from scipy import sparse; from scipy.spatial import cKDTree
R=os.environ.get("OUTPUT_ROOT", "output"); DEG=f"{R}/fig6_viewer/degdata"; SG=f"{R}/singler_recompute"; O=f"{R}/fig6_singler"
WT=["WT01","WT02","WT03","WT04","WT05","WT06","WT07","WT08"]; TMA=["tma","COEUR_TMA_G","COEUR_TMA_H"]
GLYCO=["LDHA","ALDOA","ENO1","PGK1","PDK1","SLC2A1","SLC16A3","NDRG1","VEGFA"]; RESID=["C1QC","MRC1","A2M","C7"]
MAPG=["VEGFA","CD55","SLC16A3","ICAM1","BIRC3","LGALS9","SLC2A1","CD44","TNF"]
EPI=["SecA epithelium","Intermediate epithelium","SecB epithelium"]; VASC=["Endothelial","Pericyte"]
def load(s):
    d=pd.read_csv(f"{SG}/{s}_singler_cells.csv.gz"); d=d[d.cell_label_singler.notna()&(d.cell_label_singler!="")]
    genes=np.array([l.strip() for l in open(f"{DEG}/{s}_genes.txt")]); gi={g:i for i,g in enumerate(genes)}
    cb=pd.read_csv(f"{DEG}/{s}_cells.csv.gz",usecols=["barcode"]).barcode.values
    with gzip.open(f"{DEG}/{s}_counts.mtx.gz","rb") as f: M=sparse.csr_matrix(mmread(io.BytesIO(f.read()))).T.tocsr()
    tot=np.asarray(M.sum(1)).ravel(); tot[tot==0]=1
    need=[g for g in dict.fromkeys(GLYCO+RESID+MAPG) if g in gi]
    X=pd.DataFrame(np.log1p(M[:,[gi[g] for g in need]].toarray()/tot[:,None]*1e4),columns=need,index=cb)
    d=d[d.barcode.isin(cb)].reset_index(drop=True); return d, X.loc[d.barcode.values].reset_index(drop=True)
def program(d,X):
    mac=(d.cell_label_singler=="Macrophage").values; z=lambda c:(X.loc[mac,c]-X.loc[mac,c].mean())/(X.loc[mac,c].std()+1e-9)
    g=z([c for c in GLYCO if c in X]).mean(1); r=z([c for c in RESID if c in X]).mean(1)
    p=np.full(len(d),np.nan); p[mac]=((g-g.mean())/g.std()-(r-r.mean())/r.std()).values; return p
rows=[]
def per_unit(d,prog,unit_name,cohort,sample):
    lab=d.cell_label_singler.values; xy=d[["x","y"]].values
    epi=np.isin(lab,EPI); sb=lab=="SecB epithelium"; vs=np.isin(lab,VASC); mac=lab=="Macrophage"
    if mac.sum()<30 or epi.sum()<50 or vs.sum()<5: return None
    te=cKDTree(xy[epi]); nb=te.query_ball_point(xy[mac],r=50.0); sbe=sb[epi]
    ne=np.array([len(n) for n in nb]); ns=np.array([sbe[n].sum() for n in nb])
    dv,_=cKDTree(xy[vs]).query(xy[mac]); frac=np.where(ne>=3,ns/np.maximum(ne,1),np.nan)
    return pd.DataFrame({"cohort":cohort,"sample":sample,"unit":unit_name,"barcode":d.barcode.values[mac],"program":prog[mac],
                         "secb_frac50":frac,"n_epi50":ne,"dist_vasc":dv})
for s in WT+TMA:
    d,X=load(s); prog=program(d,X)
    if s in WT:
        t=per_unit(d,prog,s,"WT",s); rows.append(t)
    else:
        d=d[d.core_id.astype(str).isin(final_core_ids(s))]   # final TMA cores
        parts=[per_unit(g,prog[g.index.values],None,"TMA",s) for c,g in d.groupby("core_id")]
        for (c,g),t in zip(d.groupby("core_id"),parts):
            if t is not None: t["unit"]=f"{s}|{g.patient_id.iloc[0]}"; rows.append(t)
    print(s,"done",flush=True)
T=pd.concat([r for r in rows if r is not None]); T.to_csv(f"{O}/fig6_colocal_percell.csv.gz",index=False)
U=[]
for u,g in T.groupby("unit"):
    g=g[np.isfinite(g.program)]
    if len(g)<30: continue
    q1,q2=np.quantile(g.program,[1/3,2/3]); lo,hi=g[g.program<=q1],g[g.program>=q2]
    # Fig 6g summary per unit = MEAN % SecB among epithelial neighbours (x100) per tertile. The per-macrophage
    # fractions are zero-inflated (in 4/8 tissues >50% of resident-high macrophages have none), so a median is
    # mostly 0 and drops the signal. Distance to vasculature (not plotted): median.
    U.append(dict(cohort=g.cohort.iloc[0],unit=u,n_mac=len(g),secb_resident=100*lo.secb_frac50.mean(),secb_hypoxic=100*hi.secb_frac50.mean(),
                  vasc_resident=lo.dist_vasc.median(),vasc_hypoxic=hi.dist_vasc.median()))
U=pd.DataFrame(U).dropna(); U.to_csv(f"{O}/fig6_colocal_units.csv",index=False)
from scipy.stats import wilcoxon
for co,g in U.groupby("cohort"):
    print(co,len(g),"SecB %% (unit mean; median over units) resident %.3f hypoxic %.3f p=%.2g"%(g.secb_resident.median(),g.secb_hypoxic.median(),wilcoxon(g.secb_resident,g.secb_hypoxic).pvalue),
          "| vasc dist resident %.1f hypoxic %.1f p=%.2g"%(g.vasc_resident.median(),g.vasc_hypoxic.median(),wilcoxon(g.vasc_resident,g.vasc_hypoxic).pvalue),
          "| n up secb",int((g.secb_hypoxic>g.secb_resident).sum()),"vasc",int((g.vasc_hypoxic>g.vasc_resident).sum()))
print("COLOCAL_DONE")
