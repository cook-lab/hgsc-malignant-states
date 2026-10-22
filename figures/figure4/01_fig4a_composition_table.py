#!/usr/bin/env python3
# PURPOSE: Fig 4a per-patient cell-type composition from the SingleR labels (cell_label_singler), one rule for every
#   arm: whole tissues (8), FTE cores (TMA A slides, sample type "fallopian") and HGSC TMA patients of the final TMA
#   core set (spatial/02_qc/06b_final_tma_cores.py: per-array core QC on SingleR labels, primary HGSC,
#   treatment-naive, >= 500 cells per patient). Patient bars with < 500 cells are dropped.
# INPUTS: singler_recompute/<sample>_singler_cells.csv.gz (WT01-WT08, tma, COEUR_TMA_G, COEUR_TMA_H); final TMA core set.
# OUTPUTS: singler_recompute/fig4a_composition_by_patient_SINGLER.csv; prints the cells per whole tissue and the
#   patients per arm (whole-tissue and FTE cell counts of the text, Figures S11-S12 and Table S8).
# MANUSCRIPT PANEL(S): Fig 4a. Plot: COMP_CSV=<OUT> PANEL_A_ONLY=1 python figures/figure4/02_fig4a_composition_plot.py
# RUNTIME TIER: moderate.
import pandas as pd, numpy as np, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from figures._prep.final_tma_cores import final_core_ids
S=os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "singler_recompute")
OUT=f"{S}/fig4a_composition_by_patient_SINGLER.csv"
WT8=["WT01","WT02","WT03","WT04","WT05","WT06","WT07","WT08"]
rows=[]
print("whole tissue (SingleR):")
for s in WT8:
    d=pd.read_csv(f"{S}/{s}_singler_cells.csv.gz",usecols=["cell_label_singler"])
    c=d["cell_label_singler"].value_counts()
    print(f"  {s:12s} {len(d):>7,} cells  ciliated={100*(d.cell_label_singler=='Ciliated epithelium').mean():.1f}%")
    for lab,n in c.items():
        rows.append({"patient_id":s,"tissue_type":"HGSC whole tissue","cell_label":lab,"n":int(n)})

def tma_arm(fname, tissue_from_sampletype, pid_prefix=""):
    d=pd.read_csv(f"{S}/{fname}_singler_cells.csv.gz",
                  usecols=["patient_id","core_id","sample_type","cell_label_singler"],dtype={"core_id":str})
    d=d[d.patient_id.notna()].copy()
    # HGSC: final TMA cores only; FTE cores (sample_type "fallopian") are kept as before
    d=d[(d.sample_type=="fallopian") | d.core_id.isin(final_core_ids(fname))]
    d["patient_id"]=pid_prefix+d.patient_id.astype(str).str.replace(r"\.0$","",regex=True)
    if tissue_from_sampletype:
        d["tissue_type"]=np.where(d.sample_type=="fallopian","Fallopian tube","HGSC TMA")
    else:
        d["tissue_type"]="HGSC TMA"
    return d[["patient_id","tissue_type","cell_label_singler"]].rename(columns={"cell_label_singler":"cell_label"})

print("CHUM discovery TMA (FTE + HGSC):")
chum=tma_arm("tma", True, "")
print(f"  patients: {chum.patient_id.nunique()}  ({chum[chum.tissue_type=='Fallopian tube'].patient_id.nunique()} FTE / {chum[chum.tissue_type=='HGSC TMA'].patient_id.nunique()} HGSC)")
print("COEUR TMA G+H (HGSC):")
coeur=pd.concat([tma_arm("COEUR_TMA_G",False,"COEUR-"),tma_arm("COEUR_TMA_H",False,"COEUR-")],ignore_index=True)
print(f"  patients: {coeur.patient_id.nunique()}")

tma_all=pd.concat([chum,coeur],ignore_index=True)
g=tma_all.groupby(["patient_id","tissue_type","cell_label"],as_index=False).size().rename(columns={"size":"n"})
rows.extend(g.to_dict("records"))

out=pd.DataFrame(rows)
# min-500-cells/patient filter (spatial-study QC; applies per patient x tissue_type)
MIN=500
tot=out.groupby(["patient_id","tissue_type"]).n.transform("sum")
before=out[["patient_id","tissue_type"]].drop_duplicates().shape[0]
out=out[tot>=MIN].copy()
after=out[["patient_id","tissue_type"]].drop_duplicates().shape[0]
print(f"min-{MIN} filter: {before} -> {after} patient-bars")
out.to_csv(OUT,index=False)
print("\nper-arm patients:"); print(out.groupby("tissue_type")["patient_id"].nunique().to_string())
print("wrote",OUT)
