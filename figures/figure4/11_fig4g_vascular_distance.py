#!/usr/bin/env python3
# PURPOSE: Fig 4g distance to the nearest vascular cell (pericyte or endothelial) of SecA, Intermediate and SecB cells,
#   from the SingleR labels; median per unit. WT: whole section. TMA: the final TMA core set
#   (spatial/02_qc/06b_final_tma_cores.py: per-array core QC on SingleR labels, primary HGSC, treatment-naive,
#   >= 500 cells per patient); nearest vascular cell searched within the core, then the median over the patient's cells.
# INPUTS: singler_recompute/<sample>_singler_cells.csv.gz (8 WT, tma, COEUR_TMA_G, COEUR_TMA_H); final TMA core set.
# OUTPUTS: singler_recompute/vascular_distance_summary_SINGLER.csv (cell_label, sample_id, tissue, median_dist, n),
#   drawn and tested by 12_fig4g_vascular_distance_plot.py (VASC_SUMMARY_CSV).
# MANUSCRIPT PANEL(S): Fig 4g. RUNTIME TIER: moderate.
import os, sys
import pandas as pd, numpy as np
from scipy.spatial import ConvexHull, cKDTree
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from figures._prep.final_tma_cores import final_core_ids
R=os.environ.get("OUTPUT_ROOT", "output"); S=f"{R}/singler_recompute"
OUTV=f"{S}/vascular_distance_summary_SINGLER.csv"
WT=["WT01","WT02","WT03","WT04","WT05","WT06","WT07","WT08"]
EPI=["SecA epithelium","Intermediate epithelium","SecB epithelium"]
VASC=["Pericyte","Endothelial"]
def load(s): return pd.read_csv(f"{S}/{s}_singler_cells.csv.gz")

vasc=[]
# ---- WT ----
for s in WT:
    d=load(s); d=d[d.cell_label_singler.notna()]
    vc=d[d.cell_label_singler.isin(VASC)][["x","y"]].values
    if len(vc):
        tree=cKDTree(vc)
        for e in EPI:
            ec=d[d.cell_label_singler==e][["x","y"]].values
            if len(ec):
                dd,_=tree.query(ec,k=1)
                vasc.append({"cell_label":e,"sample_id":s,"tissue":"whole_tissue",
                             "median_dist":float(np.median(dd)),"n":len(ec)})
# ---- TMA (CHUM tma + COEUR G/H) ----
specs=[("tma",""),("COEUR_TMA_G","COEUR-"),("COEUR_TMA_H","COEUR-")]
tma_vasc_core=[]
for s,prefix in specs:
    d=load(s); d=d[d.cell_label_singler.notna()].copy()
    d["core_id"]=d.core_id.astype(str)
    d=d[d.core_id.isin(final_core_ids(s))]          # final cores only (tumour, on-core, assigned)
    d["patient_id"]=prefix+d.patient_id.astype(str).str.replace(r"\.0$","",regex=True)
    for core,g in d.groupby("core_id"):
        pid=g.patient_id.iloc[0]
        vc=g[g.cell_label_singler.isin(VASC)][["x","y"]].values
        if len(vc):
            tree=cKDTree(vc)
            for e in EPI:
                ec=g[g.cell_label_singler==e][["x","y"]].values
                if len(ec):
                    dd,_=tree.query(ec,k=1)
                    tma_vasc_core.append({"cell_label":e,"patient_id":pid,"dist":dd})
# aggregate TMA vascular per patient (median over that patient's epithelial cells)
vv=pd.DataFrame(tma_vasc_core)
if len(vv):
    vv=vv.explode("dist"); vv["dist"]=vv["dist"].astype(float)
    agg=vv.groupby(["patient_id","cell_label"]).agg(median_dist=("dist","median"),n=("dist","size")).reset_index()
    for _,r in agg.iterrows():
        vasc.append({"cell_label":r.cell_label,"sample_id":r.patient_id,"tissue":"TMA",
                     "median_dist":r.median_dist,"n":int(r.n)})
pd.DataFrame(vasc).to_csv(OUTV,index=False)
vd=pd.DataFrame(vasc); print("vascular rows:",len(vd)," TMA patients:",vd[vd.tissue=='TMA'].sample_id.nunique())
print("DONE")
