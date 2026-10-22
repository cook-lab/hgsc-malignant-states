#!/usr/bin/env python3
# Recovered copy (2026-10-22). TMA B (COEUR blocks G and H) de-array: each cell -> core (inside the core circle, nearest of
# the 3 closest centres) -> COEUR patient (core rank -> patient through rank punch 1 / 2 of the clinical workbook). It ran
# inline (a `python3 -` heredoc) from 202609_final_hgsc_spatial/ on 2026-09-13 18:46 UTC and was never saved as a file;
# it is recovered from the analysis session log. Removed: a third TMA (not in the paper) from the loop; the rest is
# unchanged.
# Inputs: the raw Xenium cells.parquet of each block; data/clinical data/coeur_clinical_tmagh.xlsx; and, in the session's
# temporary folder SP (gone), meta_<tag>.json (the block's xmin / ymin) and final/state/<tag>.json, the core circles
# placed in dearray_editor.html. The same circles (centre, radius, rank) are written to <tag>_core_summary.csv.
# Outputs (data/dearray/): <tag>_cell_core_assignment.parquet (read by 03_build_sfe_new_coeur.R), <tag>_core_summary.csv.
import pandas as pd,numpy as np,json,warnings; warnings.filterwarnings("ignore")
from scipy.spatial import cKDTree
SP="/private/tmp/claude-501/-Volumes-CookLab-Sarah-202609-final-hgsc-spatial/b60ba64b-4537-4ba1-b9f7-c540a53e9c8c/scratchpad/dearray";OUT="data/dearray"
paths={"G":"data/HGSC xenium /COEUR_TMA_G/cells.parquet","H":"data/HGSC xenium /COEUR_TMA_H/cells.parquet"}
def toint(s):
    try:return int(str(s).strip())
    except:return None
clin=pd.read_excel("data/clinical data/coeur_clinical_tmagh.xlsx",sheet_name="clinical data");clin.columns=[str(c).strip() for c in clin.columns]
pcol=[c for c in clin.columns if "COEUR" in c.upper()][0]
r1=[c for c in clin.columns if c.lower().replace(" ","")=="rankpunch1"][0];r2=[c for c in clin.columns if c.lower().replace(" ","")=="rankpunch2"][0]
def lut_for(b):
    g=clin[clin['Bloc']==b];L=pd.concat([g[[pcol,r1]].rename(columns={r1:"rank"}),g[[pcol,r2]].rename(columns={r2:"rank"})]).dropna(subset=["rank"])
    L["rank"]=L["rank"].astype(int);return L.drop_duplicates("rank").set_index("rank")[pcol]
for tag in ["G","H"]:
    meta=json.load(open(f"{SP}/meta_{tag}.json"));doc=json.load(open(f"{SP}/final/state/{tag}.json"))
    cores=[c for c in doc["cores"] if c.get("status")!="removed" and not str(c.get('label','')).lower().startswith('ori')]
    cen=np.array([[c["x"],c["y"]] for c in cores]);rad=np.array([c["r"] for c in cores]);cid=[c["id"] for c in cores];labs=[str(c.get("label","")) for c in cores]
    cells=pd.read_parquet(paths[tag]);xy=cells[["x_centroid","y_centroid"]].values-np.array([meta["xmin"],meta["ymin"]])
    dist,idx=cKDTree(cen).query(xy,k=min(3,len(cen)))
    if idx.ndim==1: idx=idx[:,None];dist=dist[:,None]
    assign=np.full(len(xy),-1)
    for kk in range(idx.shape[1]):
        un=assign<0;w=dist[un,kk]<=rad[idx[un,kk]];rr=np.where(un)[0][w];assign[rr]=idx[rr,kk]
    lut=lut_for(tag) if tag in ("G","H") else None
    out=pd.DataFrame({"cell_id":cells["cell_id"].values,"x_centroid":cells["x_centroid"].values,"y_centroid":cells["y_centroid"].values,
                      "core_id":pd.array([cid[i] if i>=0 else None for i in assign],dtype="string"),
                      "core_rank":pd.array([toint(labs[i]) if i>=0 else None for i in assign],dtype="Int64"),
                      "mapped":assign>=0})
    if lut is not None:
        out["patient_coeur_tfri"]=pd.array([int(lut.loc[toint(labs[i])]) if (i>=0 and toint(labs[i]) in lut.index) else None for i in assign],dtype="Int64")
    out.to_parquet(f"{OUT}/{tag}_cell_core_assignment.parquet",index=False)
    cc=np.bincount(assign[assign>=0],minlength=len(cen))
    summ=pd.DataFrame({"core_id":cid,"core_rank":pd.array([toint(l) for l in labs],dtype="Int64"),
                       "x_um":cen[:,0]+meta["xmin"],"y_um":cen[:,1]+meta["ymin"],"radius_um":rad,"n_cells":cc})
    if lut is not None:
        summ["patient_coeur_tfri"]=pd.array([int(lut.loc[toint(l)]) if toint(l) in lut.index else None for l in labs],dtype="Int64")
    summ.to_csv(f"{OUT}/{tag}_core_summary.csv",index=False)
    npat=summ["patient_coeur_tfri"].nunique() if lut is not None else 0
    print(f"{tag}: {len(cores)} cores, {int((assign>=0).sum())} cells mapped, {npat} patients")
