#!/usr/bin/env python3
"""
Supplementary Table 8 -- Xenium spatial cohort: whole tissues, HGSC TMA patients, FTE cores (de-identified)
==========================================================================================================

Regenerates the table end-to-end from the files below. Counts use ALL analysed cells (no subsampling).

INPUTS (read-only)
    $OUTPUT_ROOT/singler_recompute/<sample>_singler_cells.csv.gz     analysed cells + SingleR labels
        whole tissues: WT01 WT02 WT03 WT04 WT05 WT06 WT07 WT08 (core 8)
        tma          : CHUM TMA (slides 1+2) -> "TMA A" (HGSC cores + 15 FTE cores)
        COEUR_TMA_G/H: COEUR TMA blocks G + H -> "TMA B"
    data/clinical data/20260831_maryam_verified_cohort.xlsx      whole-tissue stage + pre/post-chemo (verified sheet)
    data/clinical data/cooklab_hgsc_tissue_expansion.xlsx        whole-tissue pathology text (histology check only)
    data/clinical data/chum_tma_clinical_data_updated.xlsx       TMA A clinical (raw)
    data/clinical data/coeur_clinical_tmagh.xlsx (sheet "clinical data")   TMA B clinical (raw)
    2026_final_xenium_analysis/output/10_clinical/per_patient_features.csv            TMA A clinical as used by Fig 7
    202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup/
        coeur_clinical_perpatient_extracted.csv                                        TMA B clinical as used by Fig 7
    $OUTPUT_ROOT/final_panels/figure7/Fig7_ABCD_cox_2core_results.csv                       Fig 7a-d n (cross-check)
    $OUTPUT_ROOT/07_core_qc/final_tma_{cores,counts}.csv                      final TMA core set
        (spatial/02_qc/06b_final_tma_cores.py, read through figures/_prep/final_tma_cores.py): per-core QC flags,
        `final` and the reason a core is not final; every TMA analysis uses the final cores (2026-10-20)

The survival-cohort flags re-implement figures/figure7/01_fig7abcd_tma_cox_km_tables.R exactly (same inputs, same
filters) and are asserted to reproduce the n of Fig7_ABCD_cox_2core_results.csv.

PRIVACY
    Source identifiers (biobank / hospital / COEUR-TFRI / pathology accession numbers) are replaced by new study
    codes (WT01-WT08, TMA-P001..., FTE01-FTE15). The code -> source-ID key is written ONLY to PRIVATE_KEY_DIR
    (outside the repository and the output folder). No dates, names, free text or source IDs are written.
    Age is released as 10-year bands (<40 and >=80 lumped). OS/PFS are durations in months (0.1), not dates.

OUTPUTS
    $OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable8_xenium_clinical_cohort.xlsx
        sheets: README, Whole_tissues, TMA_patients, TMA_cores, TMA_core_selection, FTE_cores, Summary
        (+ one CSV per data sheet). Patient codes are unchanged by the final-core columns (same patients, same seed;
        asserted against the existing private key).
    The Whole_tissues CSV is read by figures/supplementary/17_figS11-13_xenium_celltype_maps.R and
    21_figS16-18_per_sample_gams.R (treatment status per WT code).

MANUSCRIPT ITEM(S)
    Table S8; Methods / Results cohort numbers (378 patients / 503 cores; Figure 7 OS n = 104, PFS n = 103).
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from supp_xlsx import write_table  # noqa: E402
from figures._prep.final_tma_cores import CORES_CSV, final_cores  # noqa: E402

SHARE = os.environ.get("DATA_ROOT", "/path/to/deposit")
WS = f"{SHARE}/202609_final_hgsc_spatial"
CLIN = f"{WS}/data/clinical data"
SG = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "singler_recompute")
CHUM_FIG7 = f"{SHARE}/2026_final_xenium_analysis/output/10_clinical/per_patient_features.csv"
COEUR_FIG7 = f"{WS}/exploratory/survival_polarization_by_coregroup/coeur_clinical_perpatient_extracted.csv"
FIG7_RES = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels/figure7/Fig7_ABCD_cox_2core_results.csv")
PRIVATE_KEY_DIR = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "private_keys")   # local only; never inside the repo or the publication folder
SEED = 20260928
MIN_CELLS = 500   # spatial QC rule (per patient for patient-level, per core for core-level statistics)

WT8 = ["WT01", "WT02", "WT03", "WT04", "WT05", "WT06", "WT07", "WT08"]
SEC = ["SecA epithelium", "Intermediate epithelium", "SecB epithelium"]
TMA_LABEL = {"tma": "A", "COEUR_TMA_G": "B", "COEUR_TMA_H": "B"}
ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV"}


def log(*a):
    print(*a, flush=True)


def pct(a, b):
    return round(100 * a / b, 1) if b else np.nan


def age_band(a):
    if pd.isna(a):
        return "Unknown"
    a = float(a)
    if a < 40:
        return "<40"
    if a >= 80:
        return ">=80"
    lo = int(a // 10 * 10)
    return f"{lo}-{lo + 9}"


def comp(d):
    """cell-composition columns from a frame of cell_label_singler"""
    lab = d.cell_label_singler
    n = len(lab); nsec = lab.isin(SEC).sum()
    return {"Cells profiled (all analysed cells)": n,
            "Cells with a SingleR label": int(lab.notna().sum()),
            "% secretory epithelium (of all cells)": pct(nsec, n),
            "% ciliated epithelium (of all cells)": pct((lab == "Ciliated epithelium").sum(), n),
            "% SecA (of secretory)": pct((lab == SEC[0]).sum(), nsec),
            "% Intermediate (of secretory)": pct((lab == SEC[1]).sum(), nsec),
            "% SecB (of secretory)": pct((lab == SEC[2]).sum(), nsec)}


def read_sg(nm, cols):
    return pd.read_csv(f"{SG}/{nm}_singler_cells.csv.gz", usecols=cols, dtype={"patient_id": str, "core_id": str})


key_rows = []   # private key (never written to the output folder)

# ====================================================================== 1. whole tissues
ver = pd.read_excel(f"{CLIN}/20260831_maryam_verified_cohort.xlsx")
ver.columns = [c.strip() for c in ver.columns]
tcol = [c for c in ver.columns if c.startswith("Unnamed")][0]            # the pre-/post-chemo column
ver["sample"] = ver["Tumour Bank ID"].astype(str).map(
    lambda s: "SP24_" + s.split("-")[-1] if s.startswith("SP-24") else f"OTB_{s}")
ver = ver.set_index("sample")
exp = pd.read_excel(f"{CLIN}/cooklab_hgsc_tissue_expansion.xlsx")
exp["sample"] = exp["Tumour Bank ID"].astype(str).map(
    lambda s: "SP24_" + s.split("-")[-1] if s.startswith("SP-24") else f"OTB_{s}")
exp = exp.set_index("sample")
path_cols = ["Histopathology / Pathology of diagnosis method (Ex. Hypoxia)", "Surgical pathology tumour morphology"]

wt_rows = []
for s in WT8:
    d = read_sg(s, ["cell_label_singler"])
    r = ver.loc[s]
    trt = {"pre-chemo": "Treatment-naive", "post-chemo": "Post-chemotherapy"}[str(r[tcol]).strip()]
    stage = ROMAN[int(r["Clinical Tumour stage at diagnosis"])] + str(r["Clinical Tumour sub stage"]).strip().upper()
    txt = " ".join(str(exp.at[s, c]) for c in path_cols).lower()
    hist = "HGSC" if ("high grade serous" in txt or "hgsc" in txt or "hgrso" in txt) else "to verify"
    wt_rows.append(dict(_src=s, **{"Treatment status at sampling": trt, "FIGO stage at diagnosis": stage,
                                   "Histology": hist, "Grade": "High grade" if hist == "HGSC" else "to verify",
                                   "Age band (years)": "Not available", "Tissue site of section": "Not available"},
                        **comp(d)))
wt = pd.DataFrame(wt_rows)
# codes: ordered by treatment status then stage (not by source ID)
wt["_o"] = wt["Treatment status at sampling"].map({"Treatment-naive": 0, "Post-chemotherapy": 1})
wt = wt.sort_values(["_o", "FIGO stage at diagnosis", "Cells profiled (all analysed cells)"]).reset_index(drop=True)
wt.insert(0, "Tissue code", [f"WT{i + 1:02d}" for i in range(len(wt))])
key_rows += [("Whole tissue", c, s) for c, s in zip(wt["Tissue code"], wt._src)]
wt = wt.drop(columns=["_src", "_o"])
assert (wt["Cells with a SingleR label"] == wt["Cells profiled (all analysed cells)"]).all()
log(wt[["Tissue code", "Treatment status at sampling", "FIGO stage at diagnosis", "Cells profiled (all analysed cells)"]])

# ====================================================================== 2. TMA cells -> cores -> patients
cols = ["cell_label_singler", "polarization_UCell", "patient_id", "core_id", "sample_type", "cell_area"]
cells = []
for nm in ["tma", "COEUR_TMA_G", "COEUR_TMA_H"]:
    d = read_sg(nm, cols)
    d["tma_obj"] = nm
    cells.append(d)
cells = pd.concat(cells, ignore_index=True)
cells["TMA"] = cells.tma_obj.map(TMA_LABEL)
n_tma_all = cells.groupby("tma_obj").size()
tum = cells[(cells.sample_type == "tumour") & cells.patient_id.notna() & cells.core_id.notna()
            & (cells.core_id != "Off core")].copy()
unassigned = cells[~cells.index.isin(tum.index) & (cells.sample_type != "fallopian")]
log("TMA cells not assigned to a patient core:", unassigned.groupby("tma_obj").size().to_dict())
tum["pkey"] = tum.TMA + ":" + tum.patient_id
core = tum.groupby(["TMA", "tma_obj", "pkey", "core_id"]).size().rename("n").reset_index()
# final TMA core set: per-array core QC (SingleR labels) + primary HGSC, treatment-naive, >= 500 cells per patient
fct = final_cores().rename(columns={"tma": "TMA", "core_id_raw": "core_id"})
core["block"] = core.tma_obj.map({"tma": "A", "COEUR_TMA_G": "G", "COEUR_TMA_H": "H"})
core = core.merge(fct[["TMA", "block", "core_id", "n_cells", "n_secretory", "n_flags", "recommendation", "final",
                       "reason_excluded"] + [c for c in fct.columns if c.startswith("flag_")]],
                  on=["TMA", "block", "core_id"], how="left", validate="1:1")
assert len(core) == len(fct) and core.final.notna().all() and (core.n == core.n_cells).all(), \
    "final core table and the cell tables describe different cores"
core["qc_pass"] = core.recommendation != "exclude"
core["final"] = core.final.astype(bool)
core["n_final_cells"] = core.n.where(core.final, 0)
pat = core.groupby(["TMA", "pkey"]).agg(n_cores=("n", "size"), n_cores_qc=("qc_pass", "sum"), n_final=("final", "sum"),
                                        cells=("n", "sum"), cells_final=("n_final_cells", "sum"),
                                        blocks=("tma_obj", lambda x: ",".join(sorted(set(x)))))
pat = pat.reset_index()
assert not pat.blocks.str.contains(",").any(), "a TMA B patient spans blocks G and H"
# composition on the final cores (the cells every TMA analysis uses); empty for patients outside the final cohort
tum = tum.merge(core[["tma_obj", "core_id", "final"]], on=["tma_obj", "core_id"], how="left")
lab = tum[tum.final].groupby("pkey").cell_label_singler
pat = pat.merge(pd.DataFrame({
    "nsec": lab.apply(lambda s: s.isin(SEC).sum()), "secA": lab.apply(lambda s: (s == SEC[0]).sum()),
    "secB": lab.apply(lambda s: (s == SEC[2]).sum())}).reset_index(), on="pkey", how="left")
pat[["nsec", "secA", "secB"]] = pat[["nsec", "secA", "secB"]].fillna(0).astype(int)
pat["src_id"] = pat.pkey.str.split(":").str[1]

# ---------------------------------------------------------------------- 2a. clinical (raw workbooks)
chum = pd.read_excel(f"{CLIN}/chum_tma_clinical_data_updated.xlsx")
chum["src_id"] = chum.patient_id.astype(str)
ppf = pd.read_csv(CHUM_FIG7, dtype={"patient_id": str})
res_map = ppf.dropna(subset=["residual_raw"]).drop_duplicates("residual_raw").set_index("residual_raw").residual_category
unmapped = set(chum.residual.dropna()) - set(res_map.index)
assert not unmapped, f"CHUM residual strings without a category: {unmapped}"
RES_LABEL = {"R0": "None (R0)", "optimal_le1cm": "<=1 cm", "suboptimal_1to2cm": ">1-2 cm",
             "suboptimal_gt2cm": ">2 cm", "suboptimal_miliary": "Miliary", "residual_unspecified": "Present, size unknown"}
DEBULK = {"None (R0)": "Optimal", "<=1 cm": "Optimal", ">1-2 cm": "Suboptimal", ">2 cm": "Suboptimal",
          "Miliary": "Suboptimal", "Suboptimal (size unknown)": "Suboptimal"}
a = pd.DataFrame({
    "src_id": chum.src_id, "TMA": "A", "age": chum.age,
    "stage": chum.stage.map(lambda x: ROMAN.get(int(x)) if pd.notna(x) else np.nan),
    "residual": chum.residual.map(res_map).map(RES_LABEL),
    "treat": chum.treatment_status.map({"naïve": "Treatment-naive", "post-chemotherapy": "Post-chemotherapy (NACT)"}),
    "hist": pd.Series("HGSC", index=chum.index).where(chum.stage.notna()), "grade": np.nan,
    "os_m": chum.survival_months, "os_e": chum.survival_outcome,
    "pfs_m": chum.progression_months, "pfs_e": chum.progression_outcome})
# cross-check the raw-derived values against the Fig 7 table where both exist
chk = a.merge(ppf, left_on="src_id", right_on="patient_id")
assert (chk.os_m == chk.survival_months).all() and (chk.age_x == chk.age_y).all()
assert (chk.residual.map(DEBULK).fillna("NA") == chk.residual_binary.str.capitalize().fillna("NA")).all()

co = pd.read_excel(f"{CLIN}/coeur_clinical_tmagh.xlsx", sheet_name="clinical data")
co = co[co.Bloc.isin(["G", "H"])].copy()
co["src_id"] = co["COEUR-TFRI"].astype(str)


def num(x):
    return pd.to_numeric(x, errors="coerce")


RES_CO = {0: "None (R0)", 1: "<=1 cm", 2: ">1-2 cm", 3: ">2 cm", 4: "Miliary", 5: "Suboptimal (size unknown)"}
HIST_CO = {1: "HGSC", 2: "Endometrioid", 3: "Low-grade serous", 4: "Clear cell", 5: "Mucinous", 6: "Mixed"}
b = pd.DataFrame({
    "src_id": co.src_id, "Bloc": co.Bloc, "TMA": "B", "age": num(co.Age_at_Diagnosis),
    "stage": num(co.FIGO_Stage_Code).map(ROMAN),
    "residual": num(co["Residual disease.1"]).map(RES_CO),
    "treat": num(co.Post_Chemo_Code).map({0: "Treatment-naive", 1: "Post-chemotherapy (NACT)"}),
    "hist": num(co.Histopathologyrevision).map(HIST_CO),
    "grade": num(co.Tumour_Grade).map({1: "G1", 2: "G2", 3: "G3"}),
    "os_m": num(co.Overall_Survival_Time), "os_e": num(co.Vital_Status_Code),
    "pfs_m": num(co.Progression_Overall), "pfs_e": num(co.Progression_Status_Code)})
assert not b.duplicated(["src_id", "Bloc"]).any()
ext = pd.read_csv(COEUR_FIG7, dtype={"patient": str})
chk = b.merge(ext, left_on="src_id", right_on="patient")
for x, y in [("os_m", "os_time"), ("pfs_m", "pfs_time"), ("age_x", "age_y")]:
    dd = (chk[x] - num(chk[y])).abs()
    assert (dd.fillna(0) < 1e-6).all(), f"COEUR {x} differs from the Fig 7 extract"
# Fig-7 debulking (0/1) == residual-derived dichotomy
assert (chk.residual.map(DEBULK).map({"Optimal": 0, "Suboptimal": 1}).fillna(-1)
        == num(chk.debulking).fillna(-1)).all()
b = b.drop(columns="Bloc")

clin = pd.concat([a, b], ignore_index=True)
clin["pkey"] = clin.TMA + ":" + clin.src_id
assert not clin.pkey.duplicated().any()
pat = pat.merge(clin.drop(columns=["TMA", "src_id"]), on="pkey", how="left")
missing_clin = ~pat.pkey.isin(set(clin.pkey))
log("TMA patients without a clinical row:", int(missing_clin.sum()))

# ---------------------------------------------------------------------- 2b. Fig 7 survival cohort (re-implemented)
fchum = ppf[["patient_id", "os_time_5y", "os_event_5y", "pfs_time_5y", "pfs_event_5y", "stage", "residual_binary", "age"]].copy()
fchum = pd.DataFrame({"pkey": "A:" + fchum.patient_id, "os_time": fchum.os_time_5y, "os_event": fchum.os_event_5y,
                      "pfs_time": fchum.pfs_time_5y, "pfs_event": fchum.pfs_event_5y,
                      "stage_late": (num(fchum.stage) >= 3).astype(float).where(num(fchum.stage).notna()),
                      "debulk": fchum.residual_binary.where(fchum.residual_binary.isin(["optimal", "suboptimal"])),
                      "age_f": num(fchum.age)})
e = ext[num(ext.os_event).isin([0, 1]) & ext.os_time.notna()].copy()   # R: character %in% c(0,1)
ot, pt = num(e.os_time), num(e.pfs_time)
fco = pd.DataFrame({"pkey": "B:" + e.patient,
                    "os_event": ((num(e.os_event) == 1) & (ot <= 60)).astype(int), "os_time": ot.clip(upper=60),
                    "pfs_event": ((num(e.pfs_event) == 1) & (pt <= 60)).astype(int), "pfs_time": pt.clip(upper=60),
                    "stage_late": (num(e.figo_code) >= 3).astype(float).where(num(e.figo_code).notna()),
                    "debulk": num(e.debulking).map({0: "optimal", 1: "suboptimal"}), "age_f": num(e.age)})
fclin = pd.concat([fchum, fco], ignore_index=True)
f7 = pat[["pkey", "n_final", "nsec"]].merge(fclin, on="pkey", how="inner")
f7 = f7[f7.n_final >= 2]                          # >= 2 final cores (01_fig7abcd_tma_cox_km_tables.R, 2026-10-20)
# primary HGSC, treatment-naive at sampling (as 01_fig7abcd_tma_cox_km_tables.R, 2026-10-17): drop neoadjuvant chemotherapy and
# reviewed non-HGSC histology; unreviewed histology is kept
excl7 = set(pat.pkey[(pat["treat"] == "Post-chemotherapy (NACT)") | (pat["hist"].notna() & (pat["hist"] != "HGSC"))])
f7 = f7[~f7.pkey.isin(excl7)]
cov = f7.stage_late.notna() & f7.debulk.notna() & f7.age_f.notna()
os_ok = cov & np.isfinite(f7.os_time) & f7.os_event.notna() & (f7.os_time > 0)
pfs_ok = cov & np.isfinite(f7.pfs_time) & f7.pfs_event.notna() & (f7.pfs_time > 0)
res7 = pd.read_csv(FIG7_RES)
n_os7 = int(res7[res7.endpoint == "OS_5yr"].n.max()); n_pfs7 = int(res7[res7.endpoint == "PFS_5yr"].n.max())
N_2CORE_MERGED = len(f7)
log(f"Fig 7: >=2-core merged {N_2CORE_MERGED}; OS {int(os_ok.sum())} (fig {n_os7}); PFS {int(pfs_ok.sum())} (fig {n_pfs7})")
assert int(os_ok.sum()) == n_os7 and int(pfs_ok.sum()) == n_pfs7, "survival cohort not reproduced"
pat["in_os"] = pat.pkey.isin(set(f7.pkey[os_ok]))
pat["in_pfs"] = pat.pkey.isin(set(f7.pkey[pfs_ok]))
N_OS_EV = int(f7.os_event[os_ok].sum()); N_PFS_EV = int(f7.pfs_event[pfs_ok].sum())
# patients with >=2 cores and complete clinical data in the raw workbook but absent from the Fig 7 CHUM table
pat["chum_missing_from_fig7"] = (pat.TMA == "A") & ~pat.pkey.isin(fchum.pkey)

# ---------------------------------------------------------------------- 2c. codes + public table
rng = np.random.RandomState(SEED)
pat = pat.sort_values(["TMA", "pkey"]).reset_index(drop=True)
order = np.concatenate([rng.permutation(np.where(pat.TMA == t)[0]) for t in ["A", "B"]])
pat = pat.iloc[order].reset_index(drop=True)
pat.insert(0, "Patient code", [f"TMA-P{i + 1:03d}" for i in range(len(pat))])
key_rows += [(f"TMA {t}", c, s) for t, c, s in zip(pat.TMA, pat["Patient code"], pat.src_id)]


def ev(x):
    return x.map({1: "Event", 0: "Censored"}) if x is not None else x


pat["in_final"] = pat.n_final > 0
pat["why_not_final"] = np.select(
    [pat.in_final, pat["treat"] == "Post-chemotherapy (NACT)", pat["hist"].notna() & (pat["hist"] != "HGSC"),
     pat.n_cores_qc == 0],
    ["", "Neoadjuvant chemotherapy before sampling", "Not HGSC on central review", "All cores excluded by core QC"],
    default=f"<{MIN_CELLS} cells on QC-passing cores")


tp = pd.DataFrame({
    "Patient code": pat["Patient code"],
    "TMA": "TMA " + pat.TMA,
    "Tumour cores profiled": pat.n_cores,
    "Cores passing core QC": pat.n_cores_qc,
    "Final cores (used in all TMA analyses)": pat.n_final,
    "Cells profiled (on-core, all analysed cells)": pat.cells,
    "Cells on final cores": pat.cells_final,
    "In final TMA cohort": np.where(pat.in_final, "Yes", "No"),
    "Reason not in final TMA cohort": pat.why_not_final,
    "Age band at diagnosis (years)": pat.age.map(age_band),
    "FIGO stage": pat.stage.fillna("Unknown"),
    "Residual disease after primary surgery": pat.residual.fillna("Unknown"),
    "Debulking (optimal <=1 cm)": pat.residual.map(DEBULK).fillna("Unknown"),
    "Treatment status of sampled tumour": pat.treat.fillna("Unknown"),
    "Histology (central review)": pat["hist"].fillna("Unknown"),
    "Tumour grade": pat.grade.fillna("Not recorded"),
    "OS (months)": pat.os_m.round(1),
    "OS status": pat.os_e.map({1: "Dead", 0: "Alive"}).fillna("Unknown"),
    "PFS (months)": pat.pfs_m.round(1),
    "PFS status": pat.pfs_e.map({1: "Progressed", 0: "Progression-free"}).fillna("Unknown"),
    "In Fig 7 OS analysis (n=%d)" % n_os7: np.where(pat.in_os, "Yes", "No"),
    "In Fig 7 PFS analysis (n=%d)" % n_pfs7: np.where(pat.in_pfs, "Yes", "No"),
    "% SecA (of secretory, final cores)": [pct(x, y) for x, y in zip(pat.secA, pat.nsec)],
    "% SecB (of secretory, final cores)": [pct(x, y) for x, y in zip(pat.secB, pat.nsec)],
})
tp.loc[pat.TMA == "A", "Tumour grade"] = "Not recorded"
tp.loc[(pat.TMA == "A") & pat.stage.isna(), "Histology (central review)"] = "Unknown"

# ---------------------------------------------------------------------- 2d. per-core sheet + core selection flow
code_of_pkey = dict(zip(pat.pkey, pat["Patient code"]))
core["Patient code"] = core.pkey.map(code_of_pkey)
core["_cnum"] = pd.to_numeric(core.core_id, errors="coerce")
core = core.sort_values(["Patient code", "_cnum", "core_id"]).reset_index(drop=True)
core["Core"] = core.groupby("Patient code").cumcount() + 1
yn = lambda x: np.where(x.astype(bool), "Yes", "No")   # noqa: E731
pat_why = dict(zip(pat.pkey, pat.why_not_final))
WHY_CORE = {"core QC (>= 2 flags)": "Core QC (>= 2 flags)"}
ct = pd.DataFrame({
    "Patient code": core["Patient code"], "TMA": "TMA " + core.TMA, "Block": core.block, "Core": core.Core,
    "Cells": core.n, "Secretory cells": core.n_secretory,
    "Flag: <200 cells": yn(core.flag_small), "Flag: <10% epithelial": yn(core.flag_low_epi),
    "Flag: low SingleR confidence": yn(core.flag_low_confidence),
    "Flag: composition outlier": yn(core.flag_composition_outlier),
    "Flag: discordant duplicate core": yn(core.flag_discordant_pair),
    "Flags (n)": core.n_flags, "Core QC (>=2 flags = exclude)": core.recommendation,
    "Final core (used in all TMA analyses)": yn(core.final),
    "Reason not final": np.where(core.final, "",
                                 np.where(~core.qc_pass, "Core QC (>= 2 flags)", core.pkey.map(pat_why)))})
assert ct["Patient code"].notna().all()
assert ct.loc[ct["Reason not final"] != "", "Final core (used in all TMA analyses)"].eq("No").all()
flow = pd.read_csv(os.path.join(os.path.dirname(CORES_CSV), "final_tma_counts.csv"))
fl = flow[flow.step.str.startswith("final (") & (flow.TMA == "A+B")].iloc[0]
assert (fl.patients, fl.cores, fl.cells) == (int(pat.in_final.sum()), int(core.final.sum()), int(pat.cells_final.sum())), \
    "final_tma_counts.csv and this table disagree"
flow = flow.rename(columns={"step": "Step", "patients": "Patients", "cores": "Cores", "cells": "Cells",
                            "secretory_cells": "Secretory cells"})
flow["TMA"] = flow.TMA.map({"A": "TMA A", "B": "TMA B", "A+B": "TMA A+B"})
STEP = {"profiled tumour cores": "1. Tumour cores profiled",
        "after core QC": "2. Cores passing core QC (< 2 flags)",
        "after core QC, HGSC, treatment-naive": "3. ... from patients with primary HGSC, treatment-naive at sampling",
        "final (patients with >= 500 cells on retained cores)": "4. Final: ... patients with >= 500 cells on these cores",
        "final patients with >= 2 final cores": "Final patients with 2 final cores (Fig 7a-d before covariates)"}
assert set(flow.Step) == set(STEP)
flow["Step"] = flow.Step.map(STEP)

# ====================================================================== 3. FTE cores
fte = cells[cells.sample_type == "fallopian"]
fc = pd.DataFrame([dict(patient_id=p_, core_id=k_, **comp(g)) for (p_, k_), g in fte.groupby(["patient_id", "core_id"])])
perm = np.random.RandomState(SEED + 1).permutation(len(fc))
fc = fc.iloc[perm].reset_index(drop=True)
fc.insert(0, "FTE core code", [f"FTE{i + 1:02d}" for i in range(len(fc))])
key_rows += [("FTE core", c, f"patient {p} / core {k}") for c, p, k in zip(fc["FTE core code"], fc.patient_id, fc.core_id)]
ftab = fc[["FTE core code", "Cells profiled (all analysed cells)", "% secretory epithelium (of all cells)",
           "% ciliated epithelium (of all cells)", "% SecA (of secretory)", "% SecB (of secretory)"]].copy()
ftab.insert(1, "TMA", "TMA A (CHUM)")
for c in ftab.columns[2:]:
    ftab[c] = pd.to_numeric(ftab[c])
ftab["Cells profiled (all analysed cells)"] = ftab["Cells profiled (all analysed cells)"].astype(int)

# ====================================================================== 4. Summary
def dist(s, order=None):
    vc = s.value_counts()
    keys = order if order else list(vc.index)
    return "; ".join(f"{k} {int(vc.get(k, 0))}" for k in keys if vc.get(k, 0))


def med_iqr(x):
    x = pd.Series(x).dropna()
    return f"{x.median():.0f} ({x.quantile(.25):.0f}-{x.quantile(.75):.0f}); n={len(x)}" if len(x) else "NA"


STG = ["I", "II", "III", "IV", "Unknown"]
coh = pat.in_final
groups = [("TMA A+B, all profiled patients", pat.cells > 0),
          ("TMA A (CHUM), final cohort", (pat.TMA == "A") & coh), ("TMA B (COEUR G+H), final cohort", (pat.TMA == "B") & coh),
          ("TMA A+B, final cohort (all TMA analyses)", coh),
          ("Fig 7 OS cohort (>=2 final cores, complete covariates)", pat.in_os), ("Fig 7 PFS cohort", pat.in_pfs)]
rows = []
for name, m in groups:
    p, t_ = pat[m], tp[m.values]
    rows.append({"Cohort": name, "Patients": int(m.sum()), "Tumour cores": int(p.n_cores.sum()),
                 "Cores passing core QC": int(p.n_cores_qc.sum()), "Final cores": int(p.n_final.sum()),
                 "Cells": int(p.cells.sum()), "Cells on final cores": int(p.cells_final.sum()),
                 "Patients with 2 final cores": int((p.n_final >= 2).sum()),
                 "Age at diagnosis, median (IQR)": med_iqr(p.age),
                 "FIGO stage": dist(t_["FIGO stage"], STG),
                 "Treatment status": dist(t_["Treatment status of sampled tumour"]),
                 "Residual disease": dist(t_["Residual disease after primary surgery"]),
                 "Histology": dist(t_["Histology (central review)"]),
                 "Deaths (any follow-up)": int((p.os_e == 1).sum()), "Progressions (any follow-up)": int((p.pfs_e == 1).sum()),
                 "5-yr OS events (as analysed)": N_OS_EV if name.startswith("Fig 7 OS") else "",
                 "5-yr PFS events (as analysed)": N_PFS_EV if name.startswith("Fig 7 PFS") else ""})
rows.insert(0, {"Cohort": "Whole tissues (core 8)", "Patients": len(wt), "Tumour cores": "",
                "Cells": int(wt["Cells profiled (all analysed cells)"].sum()),
                "FIGO stage": dist(wt["FIGO stage at diagnosis"].str.extract(r"^(IV|III|II|I)")[0], STG),
                "Treatment status": dist(wt["Treatment status at sampling"]), "Histology": dist(wt.Histology),
                "Age at diagnosis, median (IQR)": "Not available"})
rows.append({"Cohort": "FTE cores (TMA A)", "Patients": len(ftab), "Tumour cores": "",
             "Cells": int(ftab["Cells profiled (all analysed cells)"].sum())})
summ = pd.DataFrame(rows, columns=list(rows[1].keys()))   # column order of the TMA rows (the first row lacks some)

# ====================================================================== 5. re-identification screen (report only)
QI = ["TMA", "Age band at diagnosis (years)", "FIGO stage", "Residual disease after primary surgery",
      "Treatment status of sampled tumour", "Histology (central review)", "Tumour grade"]
k = tp.groupby(QI).transform("size") if False else tp.groupby(QI)["Patient code"].transform("size")
log(f"k-anonymity on {QI}: unique rows = {(k == 1).sum()}, k<3 rows = {(k < 3).sum()} of {len(tp)}")
for c in QI[1:]:
    vc = tp.groupby("TMA")[c].value_counts()
    small = vc[vc < 3]
    if len(small):
        log(f"  small categories in {c}: {small.to_dict()}")

# ====================================================================== 6. private key (outside repo/output)
os.makedirs(PRIVATE_KEY_DIR, exist_ok=True)
new_key = pd.DataFrame(key_rows, columns=["arm", "study_code", "source_id"])
kf = f"{PRIVATE_KEY_DIR}/SuppTable8_code_key_PRIVATE.csv"
if os.path.exists(kf):   # codes are used across figures and tables: they must never change between runs
    old_key = pd.read_csv(kf, dtype=str)
    m_ = old_key.merge(new_key.astype(str), on=["arm", "study_code"], how="outer", indicator=True)
    assert (m_._merge == "both").all() and (m_.source_id_x == m_.source_id_y).all(), "study codes changed"
new_key.to_csv(kf, index=False)
log("private key ->", PRIVATE_KEY_DIR)

# ====================================================================== 7. write
n_fin = int(coh.sum())
readme = [
    "De-identified description of the Xenium spatial transcriptomics cohort. All study codes (WT01-WT08, TMA-P001..., "
    "FTE01-FTE15) are new and do not encode any hospital, biobank, COEUR-TFRI or pathology accession number; the key is "
    "held privately by the authors. No dates, names or free text are included. Ages are 10-year bands (<40 and >=80 "
    "lumped). OS/PFS are durations in months from diagnosis, rounded to 0.1.",
    f"Whole_tissues: the 8 whole HGSC sections (4 treatment-naive, 4 post-chemotherapy at sampling, per the "
    f"PI-lab verified clinical sheet). Cells = all analysed cells of each section; composition from the "
    "SingleR labels with secretory cells split into SecA/Intermediate/SecB on the frozen polarization thresholds. "
    "Age and anatomical site of the sectioned block are not in the staged clinical files (reported 'Not available').",
    f"TMA_patients: TMA A = CHUM TMA (two slides; {int((pat.TMA == 'A').sum())} HGSC patients, designed as duplicate "
    f"cores); TMA B = COEUR TMA blocks G and H ({int((pat.TMA == 'B').sum())} patients; no patient is on both blocks). "
    "Cores = tumour cores assigned to the patient (off-core cells and three block-H cores with no patient excluded). "
    f"'In final TMA cohort' marks the patients whose cores are used in every TMA analysis (n = {n_fin}, "
    f"{int(pat.n_final.sum())} final cores, {int(pat.cells_final.sum()):,} cells): cores passing core QC, from patients "
    f"with primary HGSC (central review; unreviewed kept) who were treatment-naive at sampling and have >= {MIN_CELLS} "
    "cells on those cores. Secretory composition is computed on the final cores and left empty for other patients.",
    "TMA_cores: every profiled tumour core with its core QC result. Core QC was run once per array (TMA A; TMA B "
    "blocks G and H) with the same method on the final SingleR labels: five flags (< 200 cells; < 10% epithelial "
    "cells; low SingleR confidence, i.e. median score below the array's 5th percentile or > 10% pruned labels; "
    "composition outlier, Mahalanobis distance on centred log-ratio cell-type proportions, chi-square p < 0.01; "
    "duplicate cores of one patient discordant, Jensen-Shannon divergence above the array's 90th percentile or "
    "Spearman correlation < 0.5). Cores with >= 2 flags are excluded; cores with 1 flag are kept. 'Final core' marks "
    "the cores used in every TMA analysis; 'Reason not final' gives the core QC or patient-level reason. Core = "
    "1, 2 in core-map order within the patient.",
    "TMA_core_selection: patients, cores, cells and secretory cells after each selection step, by TMA.",
    "Clinical fields: TMA A from the CHUM clinical workbook, TMA B from the COEUR clinical workbook (G/H rows). "
    "Treatment status = neoadjuvant chemotherapy before the banked/sampled tumour (COEUR Post_Chemo_Code; CHUM "
    "treatment_status). Residual disease after primary surgery; debulking dichotomy optimal = none or <=1 cm. "
    "Histology = central pathology review (COEUR); TMA A patients are HGSC per cohort definition. Grade recorded "
    "for TMA B only. OS = months from diagnosis to death/last contact; PFS = months from diagnosis to first "
    "progression/last contact (COEUR 'Progression_Overall'). Figure 7 analyses censor both at 60 months.",
    f"Fig 7 flags: final-cohort patients with >= 2 final cores and complete stage, debulking and age, as in "
    f"figures/figure7/01_fig7abcd_tma_cox_km_tables.R; reproduces OS n = {n_os7} ({N_OS_EV} 5-yr deaths) and "
    f"PFS n = {n_pfs7} ({N_PFS_EV} 5-yr progressions).",
    "FTE_cores: the 15 fallopian-tube epithelium cores on TMA A (one per donor); cells and SingleR composition.",
    "Summary: cohort-level counts; age median (IQR) is computed from exact ages but only the aggregate is shown.",
    "Built by tables/08_tableS8_xenium_cohort.py.",
]
write_table(8, "xenium_clinical_cohort",
            title="Supplementary Table 8. Xenium spatial transcriptomics cohort: whole tissues, HGSC TMA patients and FTE cores",
            readme=readme,
            sheets=[("Whole_tissues", wt), ("TMA_patients", tp), ("TMA_cores", ct), ("TMA_core_selection", flow),
                    ("FTE_cores", ftab), ("Summary", summ)])
log(summ.T.to_string())
log("TMA A patients absent from the Fig 7 CHUM clinical table (per_patient_features.csv):",
    pat.loc[pat.chum_missing_from_fig7, ["Patient code", "n_cores", "cells", "stage", "residual"]].to_dict("records"))
