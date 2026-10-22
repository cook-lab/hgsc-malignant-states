#!/usr/bin/env python3
# Recovered copy (2026-10-22). TMA B (COEUR) per-patient clinical extract: one row per COEUR-TFRI patient (first
# occurrence) from the "clinical data" sheet of coeur_clinical_tmagh.xlsx: OS time and event, PFS time
# (Progression_Overall; the 9999 sentinel -> missing) and event, FIGO stage code, debulking, age. It ran inline
# (`python3 - "$SP"`) on 2026-09-15 18:02 UTC and was never saved as a file; it is recovered from the analysis session
# log. The code is unchanged. Usage: python coeur_clinical_extract.py <out_dir>
# Output: <out_dir>/coeur_clinical_perpatient.csv, copied (18:03 UTC) to
# 202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup/coeur_clinical_perpatient_extracted.csv,
# which figures/figure7/01_fig7abcd_tma_cox_km_tables.R reads.

import openpyxl, csv, sys
SP=sys.argv[1]
p="/Volumes/CookLab/Sarah/202609_final_hgsc_spatial/data/clinical data/coeur_clinical_tmagh.xlsx"
wb=openpyxl.load_workbook(p, read_only=True, data_only=True)
ws=wb["clinical data"]; rows=list(ws.iter_rows(values_only=True))
h=[str(x).strip() if x is not None else "" for x in rows[0]]; idx={c:i for i,c in enumerate(h)}
seen={}; out=[]
def g(r,c): return r[idx[c]]
for r in rows[1:]:
    pid=g(r,"COEUR-TFRI")
    if pid is None or pid in seen: continue
    seen[pid]=1
    po=g(r,"Progression_Overall")
    # 9999 sentinel -> missing time (only 1 case); else use as PFS time (months)
    if isinstance(po,(int,float)) and po>=9990: po=None
    out.append(dict(patient=pid, os_time=g(r,"Overall_Survival_Time"), os_event=g(r,"Vital_Status_Code"),
        pfs_time=po, pfs_event=g(r,"Progression_Status_Code"),
        figo_code=g(r,"FIGO_Stage_Code"), debulking=g(r,"Debulking"), age=g(r,"Age_at_Diagnosis")))
w=csv.DictWriter(open(f"{SP}/coeur_clinical_perpatient.csv","w"), fieldnames=list(out[0].keys()))
w.writeheader(); [w.writerow(o) for o in out]
print("re-extracted with Progression_Overall as PFS time:", len(out), "patients")
