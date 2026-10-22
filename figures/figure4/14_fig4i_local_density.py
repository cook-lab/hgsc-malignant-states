#!/usr/bin/env python3
# PURPOSE: Fig 4i local cell density around SecA, Intermediate and SecB cells per sample: the per-sample test of the
# Fig 4f density curve, with the TMA as validation.
# Density = the Fig 4f definition (08_fig4f_gam_features.R): all SingleR-labelled cells within 50 um of a
# secretory cell, the cell itself excluded, per mm^2 of the 50 um circle. Per unit = median over that class's cells.
#   WT : whole section (core 8).
#   TMA: final TMA core set (figures/_prep/final_tma_cores.py); neighbours searched within the core (core coordinates
#        overlap across cores); a patient's cells are pooled over its retained cores before the median.
#   Sensitivity (TMA, not plotted): only cells >= 50 um inside their core's convex hull, so no neighbourhood is cut
#        by the core edge (columns *_core50).
# INPUTS: singler_recompute/<sample>_singler_cells.csv.gz (8 WT, tma, COEUR_TMA_G, COEUR_TMA_H); final TMA core set.
# OUTPUT: singler_recompute/per_patient_local_density_SINGLER.csv -> figure4/15_fig4i_local_density_plot.py
# MANUSCRIPT PANEL(S): Fig 4i. RUNTIME TIER: moderate.
import os, sys
import numpy as np, pandas as pd
from scipy.spatial import ConvexHull, cKDTree
from scipy.stats import wilcoxon
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from figures._prep.final_tma_cores import final_core_ids
R=os.environ.get("OUTPUT_ROOT", "output"); S=f"{R}/singler_recompute"
OUT=f"{S}/per_patient_local_density_SINGLER.csv"
WT=["WT01","WT02","WT03","WT04","WT05","WT06","WT07","WT08"]
KEY={"SecA epithelium":"SecA","Intermediate epithelium":"Int","SecB epithelium":"SecB"}
RADIUS=50.0; CIRC=np.pi*(RADIUS/1000)**2          # mm^2 of the 50 um circle
def load(s):
    d=pd.read_csv(f"{S}/{s}_singler_cells.csv.gz"); return d[d.cell_label_singler.notna()].reset_index(drop=True)
def neighbours(xy, focal):
    n=cKDTree(xy).query_ball_point(xy[focal], r=RADIUS, return_length=True)
    return np.asarray(n)-1                         # minus the focal cell itself
def hull_depth(xy):
    # distance from each point to the convex-hull boundary (unit normals; n.x + b <= 0 inside)
    if len(xy)<3: return np.zeros(len(xy))
    eq=ConvexHull(xy).equations
    return np.min(-(xy@eq[:,:2].T+eq[:,2]),axis=1)
def medians(g, suffix=""):
    # med_nb = median neighbour count (exact .0 / .5); ldens = med_nb per mm^2. Tests use med_nb: same signs and ranks,
    # and its ties and zero differences are exact (float densities read back from CSV can differ in the last bit).
    row={}
    for k in KEY.values():
        v=g.loc[g.cls==k,"nb"]; m=float(np.median(v)) if len(v) else np.nan
        row[f"med_nb_{k}{suffix}"]=m; row[f"ldens_{k}{suffix}"]=m/CIRC; row[f"n_{k}{suffix}"]=int(len(v))
    return row

rows=[]
# ---- WT ----
for s in WT:
    d=load(s); focal=d.cell_label_singler.isin(KEY).values
    c=pd.DataFrame({"cls":d.cell_label_singler[focal].map(KEY).values,
                    "nb":neighbours(d[["x","y"]].values, focal)})
    rows.append({"patient_id":s,"source":"WT",**medians(c)})
    print(s,"done",flush=True)
# ---- TMA (CHUM tma + COEUR G/H), final cores ----
cells=[]
# patient ids as '<TMA>|<number>' (TMA B blocks G and H pooled under 'COEUR')
for s,prefix in [("tma","tma|"),("COEUR_TMA_G","COEUR|"),("COEUR_TMA_H","COEUR|")]:
    d=load(s); d["core_id"]=d.core_id.astype(str)
    d=d[d.core_id.isin(final_core_ids(s))]
    d["patient_id"]=prefix+d.patient_id.astype(str).str.replace(r"\.0$","",regex=True)
    for core,g in d.groupby("core_id"):
        focal=g.cell_label_singler.isin(KEY).values
        if not focal.any(): continue
        xy=g[["x","y"]].values
        cells.append(pd.DataFrame({"patient_id":g.patient_id.values[focal],
                                   "cls":g.cell_label_singler.values[focal],
                                   "nb":neighbours(xy, focal),
                                   "depth":hull_depth(xy)[focal]}))
    print(s,"done",flush=True)
C=pd.concat(cells); C["cls"]=C.cls.map(KEY)
for p,g in C.groupby("patient_id"):
    rows.append({"patient_id":p,"source":"TMA",**medians(g),**medians(g[g.depth>=RADIUS],"_core50")})
T=pd.DataFrame(rows); T.to_csv(OUT,index=False)

for src,sfx in [("WT",""),("TMA",""),("TMA","_core50")]:
    g=T[T.source==src][[f"med_nb_SecA{sfx}",f"med_nb_SecB{sfx}"]].dropna()
    a,b=g.iloc[:,0],g.iloc[:,1]
    print(f"{src}{sfx}: n {len(g)}  median SecA {a.median()/CIRC:.0f}  SecB {b.median()/CIRC:.0f} cells/mm2  "
          f"lower around SecB {int((b<a).sum())}/{len(g)}  paired Wilcoxon p {wilcoxon(a,b).pvalue:.2g}")
print("wrote",OUT); print("DONE")
