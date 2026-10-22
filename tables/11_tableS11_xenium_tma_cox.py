#!/usr/bin/env python3
"""Supplementary Table 11 -- Xenium HGSC TMA multivariable Cox models (Figure 7a-d; Results, TMA survival).

Regenerates the table end-to-end:
  1. Runs the Figure 7 model code (figures/figure7/01_fig7abcd_tma_cox_km_tables.R, 2026-09-25 log1p
     version) VERBATIM for feature construction, clinical merge, 2-core restriction, complete-case filter,
     scaling and the coxph() call -- no covariate / cohort / endpoint changes. The only additions are
     exporting what the generator computed but did not save: every model coefficient (stage, residual
     disease, age), the feature SD used for scaling, and a cox.zph() proportional-hazards diagnostic on the
     SAME fitted model objects (diagnostic only; does not alter any estimate).
     It also computes the Figure 7a,b Kaplan-Meier comparison (lowest vs highest polarization tertile of each
     endpoint's model cohort, middle tertile left out, as for TCGA-OV in Figure 7g,h; 5-year KM survival per
     group; log-rank test).
  2. Checks the exported feature HRs against the promoted figure source
     $OUTPUT_ROOT/final_panels/figure7/Fig7_ABCD_cox_2core_results.csv, the Figure 7c,d FDR (Benjamini-Hochberg
     across the 6 secretory-measure tests) against Fig7_ABCD_forest_fdr.csv and the KM rows against
     Fig7_ABCD_km_logrank.csv (all must match to 1e-8).
  3. Writes SuppTable11 via the shared writer supp_xlsx.write_table.

Intermediates: $OUTPUT_ROOT/supp_tables/ST11/ (read-only on all inputs; nothing in final_panels/figure7/ is touched).
"""
import os, sys, subprocess, textwrap
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from supp_xlsx import write_table  # noqa: E402

REPO = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, REPO)
from figures._prep.final_tma_cores import final_patients  # noqa: E402
FIG_CSV = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels", "figure7", "Fig7_ABCD_cox_2core_results.csv")
FDR_CSV = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels", "figure7", "Fig7_ABCD_forest_fdr.csv")
KM_CSV = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "final_panels", "figure7", "Fig7_ABCD_km_logrank.csv")
WORK = os.path.join(os.environ.get("OUTPUT_ROOT", "output"), "supp_tables/ST11")
os.makedirs(WORK, exist_ok=True)

# ---------------------------------------------------------------------------------------------------------
# R: sections 1-3 and the model loop are copied from 01_fig7abcd_tma_cox_km_tables.R; OUT redirected to WORK.
# ---------------------------------------------------------------------------------------------------------
R_CODE = r'''
suppressWarnings(suppressMessages({library(data.table); library(survival)}))
source(Sys.getenv("FINAL_TMA_CORES_R"))
SG <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "singler_recompute")
SP <- file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup")
OUT<- Sys.getenv("ST11_WORK")
SECA<-"SecA epithelium"; INT<-"Intermediate epithelium"; SECB<-"SecB epithelium"; SEC<-c(SECA,INT,SECB)

## ---- 1. Per-patient features on SingleR labels (densities pooled over the patient's cores) ----
feat_list<-list()
for (nm in c("tma","COEUR_TMA_G","COEUR_TMA_H")) {
  d <- fread(file.path(SG, paste0(nm,"_singler_cells.csv.gz")))
  d <- d[!is.na(cell_label_singler) & cell_label_singler!="" & sample_type=="tumour" &
         !is.na(core_id) & core_id!="Off core" & !is.na(patient_id)]
  d <- d[as.character(core_id) %in% final_core_ids(nm)]   # final TMA cores only (2026-10-20); n_cores counts them
  ca <- d[, .(area=sum(cell_area,na.rm=TRUE)/1e6), by=.(patient_id,core_id)]
  pa <- ca[, .(area=sum(area), n_cores=uniqueN(core_id)), by=patient_id]
  cnt<- d[, .N, by=.(patient_id,cell_label_singler)]
  cw <- dcast(cnt, patient_id~cell_label_singler, value.var="N", fill=0)
  f  <- merge(pa, cw, by="patient_id")
  ctcols <- setdiff(names(cw),"patient_id")
  for (ct in ctcols) set(f, j=paste0("dens_", gsub(" ","_",ct)), value=f[[ct]]/f$area)
  sec <- d[cell_label_singler %in% SEC]
  sp  <- sec[, .(n_sec=.N, secA=sum(cell_label_singler==SECA), secB=sum(cell_label_singler==SECB),
                 polar_mean=mean(polarization_UCell,na.rm=TRUE)), by=patient_id]
  sp[, `:=`(prop_SecA_of_sec=fifelse(n_sec>0,secA/n_sec,NA_real_),
            prop_SecB_of_sec=fifelse(n_sec>0,secB/n_sec,NA_real_))]
  f <- merge(f, sp[,.(patient_id,prop_SecA_of_sec,prop_SecB_of_sec,polar_mean)], by="patient_id", all.x=TRUE)
  f[, tma:=paste0("sfe_",nm)]
  feat_list[[nm]]<-f
}
feat<-rbindlist(feat_list, fill=TRUE); feat[, patient_id:=as.character(patient_id)]
for (c0 in grep("^dens_",names(feat),value=TRUE)) feat[is.na(get(c0)), (c0):=0]

## ---- 2. Clinical (identical to generator) ----
chum<-as.data.table(read.csv(file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "2026_final_xenium_analysis/output/10_clinical/per_patient_features.csv"), stringsAsFactors=FALSE))[, .(patient_id,os_time_5y,os_event_5y,pfs_time_5y,pfs_event_5y,stage,residual_binary,age)]
chum[, `:=`(patient_id=as.character(patient_id), tma="sfe_tma",
            os_time=os_time_5y, os_event=os_event_5y, pfs_time=pfs_time_5y, pfs_event=pfs_event_5y, age=as.numeric(age),
            stage_late=as.integer(suppressWarnings(as.numeric(stage))>=3),
            debulk=fifelse(residual_binary=="optimal","optimal",fifelse(residual_binary=="suboptimal","suboptimal",NA_character_)))]
coeur<-as.data.table(read.csv(file.path(SP,"coeur_clinical_perpatient_extracted.csv"), stringsAsFactors=FALSE))
coeur<-coeur[os_event %in% c(0,1) & !is.na(os_time)]
coeur[, patient_id:=as.character(patient)]
coeur[, `:=`(os_event=as.integer(os_event==1 & as.numeric(os_time)<=60),
             pfs_event=as.integer(pfs_event==1 & suppressWarnings(as.numeric(pfs_time))<=60),
             os_time=pmin(as.numeric(os_time),60), pfs_time=pmin(suppressWarnings(as.numeric(pfs_time)),60),
             age=as.numeric(age), stage_late=as.integer(suppressWarnings(as.numeric(figo_code))>=3),
             debulk=fifelse(as.character(debulking)=="0","optimal",fifelse(as.character(debulking)=="1","suboptimal",NA_character_)))]
mkc<-function(x) x[, .(patient_id,tma,os_time,os_event,pfs_time,pfs_event,stage_late,debulk,age)]
cG<-copy(coeur)[,tma:="sfe_COEUR_TMA_G"]; cH<-copy(coeur)[,tma:="sfe_COEUR_TMA_H"]
clin<-rbind(mkc(chum),mkc(cG),mkc(cH))

## ---- 2b. Primary HGSC, treatment-naive at sampling (as 01_fig7abcd_tma_cox_km_tables.R, 2026-10-17) ----
CLIN_DIR<-file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial/data/clinical data")
cx<-as.data.table(readxl::read_excel(file.path(CLIN_DIR,"chum_tma_clinical_data_updated.xlsx")))
ex<-cx[treatment_status=="post-chemotherapy", .(patient_id=as.character(patient_id), tma="sfe_tma", why="NACT")]
cb<-as.data.table(readxl::read_excel(file.path(CLIN_DIR,"coeur_clinical_tmagh.xlsx"), sheet="clinical data"))[Bloc %in% c("G","H")]
cb[, `:=`(patient_id=as.character(`COEUR-TFRI`), tma=paste0("sfe_COEUR_TMA_",Bloc),
          nact=suppressWarnings(as.numeric(Post_Chemo_Code))==1, hist=suppressWarnings(as.numeric(Histopathologyrevision)))]
ex<-rbind(ex, cb[nact %in% TRUE, .(patient_id, tma, why="NACT")],
              cb[!is.na(hist) & hist!=1, .(patient_id, tma, why="non-HGSC")])

## ---- 3. Merge + 2-core restriction ----
d<-merge(feat, clin, by=c("patient_id","tma")); d[, debulk:=factor(debulk,levels=c("optimal","suboptimal"))]
d2<-d[n_cores>=2]
d2<-d2[!ex, on=.(patient_id, tma)]
cat(sprintf("2-core cohort: %d (CHUM %d, COEUR_G %d, COEUR_H %d)\n", nrow(d2),
  sum(d2$tma=="sfe_tma"),sum(d2$tma=="sfe_COEUR_TMA_G"),sum(d2$tma=="sfe_COEUR_TMA_H")))
fwrite(data.table(n_two_core_clinical=nrow(d2)), file.path(OUT,"st11_two_core_n.csv"))   # for the README (not in the generator)

## ---- 4. Multivariable Cox (adj: +stage+debulk+age), HR per +1 SD -- generator loop + full export ----
FEATS<-c("prop_SecA_of_sec","prop_SecB_of_sec","polar_mean")   # pre-specified secretory measures only (author 2026-10-19)
zc<-function(x){s<-sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s}
res<-list(); full<-list(); zph<-list(); coh<-list(); km<-list()
for (EPi in c("OS","PFS")){
  tcol<-if(EPi=="OS")"os_time" else "pfs_time"; ecol<-if(EPi=="OS")"os_event" else "pfs_event"
  for (fe in FEATS){
    sub<-copy(d2); covs<-c("stage_late","debulk","age")
    keep<-is.finite(sub[[fe]])&is.finite(sub[[tcol]])&!is.na(sub[[ecol]])&sub[[tcol]]>0
    for(cv in covs) keep<-keep & !is.na(sub[[cv]])
    sub<-sub[keep]; if(nrow(sub)<15||sum(sub[[ecol]])<8) next
    val <- if(grepl("^dens_",fe)) log1p(sub[[fe]]) else sub[[fe]]
    sub[, z:=zc(val)]
    fit<-tryCatch(coxph(as.formula(paste0("Surv(",tcol,",",ecol,")~z+stage_late+debulk+age")),data=sub),error=function(e)NULL)
    if(is.null(fit)) next; s<-summary(fit)
    ep<-paste0(EPi,"_5yr")
    res[[length(res)+1]]<-data.table(endpoint=ep, feature=fe, n=s$n, n_events=s$nevent,
      HR=s$conf.int["z","exp(coef)"], HR_lower=s$conf.int["z","lower .95"], HR_upper=s$conf.int["z","upper .95"],
      p_value=s$coefficients["z","Pr(>|z|)"], feat_mean=mean(val), feat_sd=sd(val),
      raw_median=median(sub[[fe]]), raw_q1=quantile(sub[[fe]],.25), raw_q3=quantile(sub[[fe]],.75),
      frac_zero=mean(sub[[fe]]==0),
      stage_late_n=sum(sub$stage_late==1), subopt_n=sum(sub$debulk=="suboptimal"), age_median=median(sub$age),
      concordance=s$concordance[1], lr_test_p=s$logtest["pvalue"])
    cf<-s$coefficients; ci<-s$conf.int
    full[[length(full)+1]]<-data.table(endpoint=ep, feature=fe, term=rownames(cf), coef=cf[,"coef"],
      se=cf[,"se(coef)"], HR=ci[,"exp(coef)"], HR_lower=ci[,"lower .95"], HR_upper=ci[,"upper .95"],
      z_stat=cf[,"z"], p_value=cf[,"Pr(>|z|)"], n=s$n, n_events=s$nevent)
    zz<-tryCatch(cox.zph(fit, transform="km"), error=function(e) NULL)
    if(!is.null(zz)) zph[[length(zph)+1]]<-data.table(endpoint=ep, feature=fe, term=rownames(zz$table),
      chisq=zz$table[,"chisq"], df=zz$table[,"df"], p_value=zz$table[,"p"])
    if (fe=="polar_mean") coh[[length(coh)+1]]<-sub[, .(n=.N, events=sum(get(ecol)),
      median_observed_time_mo=median(get(tcol)), stage_III_IV=sum(stage_late==1), suboptimal=sum(debulk=="suboptimal"),
      age_median=median(age)), by=tma][, endpoint:=ep]
    if (fe=="polar_mean") {   # Figure 7a,b: KM of the outer polarization tertiles of this model cohort (middle dropped)
      ct<-quantile(sub$polar_mean, c(1/3, 2/3)); kk<-sub[polar_mean<=ct[1] | polar_mean>ct[2]]
      kk[, grp:=factor(ifelse(polar_mean>ct[2],"SecB","SecA"), levels=c("SecA","SecB"))]
      kf<-survfit(as.formula(paste0("Surv(",tcol,",",ecol,")~grp")), data=kk); s60<-summary(kf, times=60, extend=TRUE)
      lr<-survdiff(as.formula(paste0("Surv(",tcol,",",ecol,")~grp")), data=kk)
      km[[length(km)+1]]<-data.table(endpoint=ep, n_cohort=nrow(sub), n=nrow(kk), events=sum(kk[[ecol]]),
        cut_low=unname(ct[1]), cut_high=unname(ct[2]), n_SecA=sum(kk$grp=="SecA"), n_SecB=sum(kk$grp=="SecB"),
        events_SecA=sum(kk[[ecol]][kk$grp=="SecA"]), events_SecB=sum(kk[[ecol]][kk$grp=="SecB"]),
        surv5y_SecA=s60$surv[1], surv5y_SecA_lower=s60$lower[1], surv5y_SecA_upper=s60$upper[1],
        surv5y_SecB=s60$surv[2], surv5y_SecB_lower=s60$lower[2], surv5y_SecB_upper=s60$upper[2],
        logrank_p=pchisq(lr$chisq, 1, lower.tail=FALSE))
    }
  }
}
R<-rbindlist(res); R[, q_value:=p.adjust(p_value,"BH")]   # BH across the 6 tests, as Fig 7c,d
fwrite(R, file.path(OUT,"st11_feature_rows.csv"))
fwrite(rbindlist(full), file.path(OUT,"st11_full_coefficients.csv"))
fwrite(rbindlist(zph), file.path(OUT,"st11_cox_zph.csv"))
fwrite(rbindlist(coh), file.path(OUT,"st11_cohort_by_tma.csv"))
fwrite(rbindlist(km), file.path(OUT,"st11_km_tertiles.csv"))
cat("ST11_R_DONE\n")
'''

rfile = os.path.join(WORK, "st11_models.R")
with open(rfile, "w") as fh:
    fh.write(R_CODE)
subprocess.run(["Rscript", rfile], check=True, env={**os.environ, "ST11_WORK": WORK,
               "FINAL_TMA_CORES_R": os.path.join(REPO, "spatial", "00_setup", "final_tma_cores.R")})

R = pd.read_csv(f"{WORK}/st11_feature_rows.csv")
FULL = pd.read_csv(f"{WORK}/st11_full_coefficients.csv")
ZPH = pd.read_csv(f"{WORK}/st11_cox_zph.csv")
COH = pd.read_csv(f"{WORK}/st11_cohort_by_tma.csv")
KM = pd.read_csv(f"{WORK}/st11_km_tertiles.csv")

# ---- 2. Reproduction check against the promoted figure source -------------------------------------------
fig = pd.read_csv(FIG_CSV)
chk = fig.merge(R, on=["endpoint", "feature"], suffixes=("_fig", "_st"))
assert len(chk) == len(fig) == len(R), (len(chk), len(fig), len(R))
for c in ["HR", "HR_lower", "HR_upper", "p_value", "q_value", "n", "n_events"]:
    dmax = float(np.max(np.abs(chk[f"{c}_fig"] - chk[f"{c}_st"])))
    assert dmax < 1e-8, f"{c} differs from Fig7_ABCD_cox_2core_results.csv by {dmax}"
print(f"reproduction check vs {os.path.basename(FIG_CSV)}: all {len(chk)} rows identical")

# Figure 7c,d family: Benjamini-Hochberg across the 3 secretory measures x 2 endpoints (as R p.adjust "BH")
def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)[::-1]
    q = np.empty(m); q[o] = np.minimum.accumulate(p[o] * m / np.arange(m, 0, -1))
    return np.minimum(q, 1)
SEC_F = ["prop_SecA_of_sec", "prop_SecB_of_sec", "polar_mean"]
sec = R[R.feature.isin(SEC_F)].copy()
assert len(sec) == 6, len(sec)
sec["q_fig7c"] = bh(sec.p_value.values)
assert float(np.max(np.abs(sec.q_fig7c.values - sec.q_value.values))) < 1e-8   # the model q is already this family
R = R.merge(sec[["endpoint", "feature", "q_fig7c"]], on=["endpoint", "feature"], how="left")
fdr = pd.read_csv(FDR_CSV)
fdr["feature"] = fdr.measure.map({"SecA proportion": "prop_SecA_of_sec", "SecB proportion": "prop_SecB_of_sec",
                                  "Polarization": "polar_mean"})
c7 = fdr.merge(sec, on=["endpoint", "feature"], suffixes=("_fig", "_st"))
assert len(c7) == 6
for a, b in [("q_fdr", "q_fig7c"), ("HR_fig", "HR_st"), ("p_value_fig", "p_value_st")]:
    assert float(np.max(np.abs(c7[a] - c7[b]))) < 1e-8, f"{a} differs from Fig7_ABCD_forest_fdr.csv"
km_fig = pd.read_csv(KM_CSV)
ck = km_fig.merge(KM, on="endpoint", suffixes=("_fig", "_st"))
assert len(ck) == len(km_fig) == len(KM) == 2
for c in [c for c in KM.columns if c != "endpoint"]:
    assert float(np.max(np.abs(ck[f"{c}_fig"] - ck[f"{c}_st"]))) < 1e-8, f"KM {c} differs from Fig7_ABCD_km_logrank.csv"
print("reproduction check vs Fig7_ABCD_forest_fdr.csv (6 rows) and Fig7_ABCD_km_logrank.csv (2 rows): identical")

# ---- 3. Format ----------------------------------------------------------------------------------------
SECRETORY = {"prop_SecA_of_sec": "SecA proportion of secretory cells",
             "prop_SecB_of_sec": "SecB proportion of secretory cells",
             "polar_mean": "Mean polarization score (UCell)"}
def label(f): return SECRETORY.get(f, f.replace("dens_", "").replace("_", " ") + " density")
def fclass(f): return "Secretory measure" if f in SECRETORY else "Cell-type density"
def scaling(f):
    if f.startswith("dens_"):
        return "log1p(cells/mm2), z-scored; HR per +1 SD"
    if f == "polar_mean":
        return "raw mean cell-level polarization UCell score, z-scored; HR per +1 SD"
    return "raw proportion (0-1) of secretory cells, z-scored; HR per +1 SD"
def ci(lo, hi): return f"{lo:.2f}-{hi:.2f}"
EP = {"OS_5yr": "OS (5-year)", "PFS_5yr": "PFS (5-year)"}
fe_order = {f: i for i, f in enumerate(
    ["prop_SecA_of_sec", "prop_SecB_of_sec", "polar_mean"] + sorted(f for f in R.feature.unique() if f.startswith("dens_")))}
TERM = {"z": "feature (per +1 SD)", "stage_late": "FIGO stage III/IV vs I/II",
        "debulksuboptimal": "Residual disease: suboptimal vs optimal debulking", "age": "Age (per year)"}

full_w = FULL.copy()
full_w["term_key"] = full_w.term.map({"z": "feat", "stage_late": "stage", "debulksuboptimal": "resid", "age": "age"})
wide = full_w.pivot_table(index=["endpoint", "feature"], columns="term_key", values=["HR", "HR_lower", "HR_upper", "p_value"])
zph_w = ZPH.pivot_table(index=["endpoint", "feature"], columns="term", values="p_value")

rows = []
for r in R.itertuples():
    k = (r.endpoint, r.feature)
    w = wide.loc[k]
    row = {"endpoint": EP[r.endpoint], "feature": label(r.feature),
           "feature_variable": r.feature, "feature_scaling": scaling(r.feature),
           "SD_of_modelled_variable": round(r.feat_sd, 4),
           "HR": round(r.HR, 3), "CI95": ci(r.HR_lower, r.HR_upper),
           "HR_lower": round(r.HR_lower, 3), "HR_upper": round(r.HR_upper, 3),
           "p_value": float(f"{r.p_value:.3g}"), "FDR_BH": float(f"{r.q_value:.3g}"),
           "n": int(r.n), "events": int(r.n_events)}
    for key, nm in [("stage", "stage_III_IV"), ("resid", "residual_suboptimal"), ("age", "age_per_year")]:
        row[f"{nm}_HR"] = round(w[("HR", key)], 3)
        row[f"{nm}_CI95"] = ci(w[("HR_lower", key)], w[("HR_upper", key)])
        row[f"{nm}_p"] = float(f"{w[('p_value', key)]:.3g}")
    row["PH_test_feature_p"] = float(f"{zph_w.loc[k, 'z']:.3g}")
    row["PH_test_global_p"] = float(f"{zph_w.loc[k, 'GLOBAL']:.3g}")
    row["_o"] = (r.endpoint, fe_order[r.feature])
    rows.append(row)
main = pd.DataFrame(rows).sort_values("_o").drop(columns="_o").reset_index(drop=True)

fm = FULL.copy()
fm["term"] = fm.term.map(TERM)
fm.insert(1, "feature_label", fm.feature.map(label))
fm["endpoint"] = fm.endpoint.map(EP)
fm["_o"] = [(e, fe_order[f]) for e, f in zip(FULL.endpoint, FULL.feature)]
fm = fm.sort_values("_o", kind="stable").drop(columns="_o")
for c in ["coef", "se", "HR", "HR_lower", "HR_upper", "z_stat"]: fm[c] = fm[c].round(4)
fm["p_value"] = fm.p_value.map(lambda v: float(f"{v:.3g}"))
fm = fm.rename(columns={"feature": "feature_variable", "feature_label": "feature"})

zp = ZPH.copy()
zp["term"] = zp.term.map({**TERM, "debulk": TERM["debulksuboptimal"], "GLOBAL": "GLOBAL (all terms)"})   # cox.zph names the factor "debulk"
zp.insert(1, "feature_label", zp.feature.map(label))
zp["endpoint"] = zp.endpoint.map(EP)
zp["_o"] = [(e, fe_order[f]) for e, f in zip(ZPH.endpoint, ZPH.feature)]
zp = zp.sort_values("_o", kind="stable").drop(columns="_o")
zp["chisq"] = zp.chisq.round(3); zp["p_value"] = zp.p_value.map(lambda v: float(f"{v:.3g}"))
zp["PH_violation_p_lt_0.05"] = zp.p_value < 0.05
zp = zp.rename(columns={"feature": "feature_variable", "feature_label": "feature"})

coh = COH.copy()
coh["endpoint"] = coh.endpoint.map(EP)
coh["tma"] = coh.tma.map({"sfe_tma": "CHUM HGSC TMA", "sfe_COEUR_TMA_G": "COEUR TMA (Bloc G)",
                          "sfe_COEUR_TMA_H": "COEUR TMA (Bloc H)"})
tot = coh.groupby("endpoint", as_index=False)[["n", "events", "stage_III_IV", "suboptimal"]].sum()
tot["tma"] = "Total"
coh = pd.concat([coh, tot], ignore_index=True)[["endpoint", "tma", "n", "events", "median_observed_time_mo",
                                                "stage_III_IV", "suboptimal", "age_median"]]
coh = coh.sort_values(["endpoint", "tma"], key=lambda s: s.map(lambda v: "zz" if v == "Total" else v)).reset_index(drop=True)
coh["median_observed_time_mo"] = coh.median_observed_time_mo.round(1)

def pct(v, lo, hi): return f"{100 * v:.0f}% ({100 * lo:.0f}-{100 * hi:.0f})"
kmt = pd.DataFrame({
    "endpoint": KM.endpoint.map(EP), "n_model_cohort": KM.n_cohort,
    "lower_tertile_cut_polarization": KM.cut_low.round(4), "upper_tertile_cut_polarization": KM.cut_high.round(4),
    "n_compared": KM.n, "events": KM.events, "n_SecA_dominant": KM.n_SecA, "events_SecA_dominant": KM.events_SecA,
    "n_SecB_enriched": KM.n_SecB, "events_SecB_enriched": KM.events_SecB,
    "KM_5yr_SecA_dominant_95CI": [pct(*v) for v in KM[["surv5y_SecA", "surv5y_SecA_lower", "surv5y_SecA_upper"]].values],
    "KM_5yr_SecB_enriched_95CI": [pct(*v) for v in KM[["surv5y_SecB", "surv5y_SecB_lower", "surv5y_SecB_upper"]].values],
    "logrank_p": KM.logrank_p.map(lambda v: float(f"{v:.3g}"))})

n_ph = int((ZPH.term == "z").sum()); n_ph_viol = int(((ZPH.term == "z") & (ZPH.p_value < 0.05)).sum())
FP = final_patients()
N_FIN, N2_ALL = int(FP.final.sum()), int(FP.two_final_cores.sum())
N2 = int(pd.read_csv(f"{WORK}/st11_two_core_n.csv").n_two_core_clinical.iloc[0])
# Cox (complete-case) cohort per endpoint; the KM table counts only the two outer tertiles since 2026-10-20
CO = R[R.feature == "polar_mean"].set_index("endpoint")
KO, KP = CO.loc["OS_5yr"], CO.loc["PFS_5yr"]
g_viol = int(((ZPH.term == "GLOBAL") & (ZPH.p_value < 0.05)).sum())

readme = [
    "Cohort: the final TMA core set used by every TMA analysis in the paper (Table S8): tumour cores on the Xenium TMAs (CHUM HGSC TMA; COEUR TMA Blocs G and H) "
    "that pass core QC (run per array with the same method on the SingleR labels), from patients with primary HGSC on central review (unreviewed kept) "
    "who were treatment-naive at sampling (no neoadjuvant chemotherapy) and have >= 500 cells on those cores. "
    f"{N_FIN} patients; {N2_ALL} have >= 2 final cores ({N2} with clinical follow-up). Complete-case (time > 0, stage, residual disease and age available): "
    f"n = {int(KO.n)} for OS ({int(KO.n_events)} events) and n = {int(KP.n)} for PFS ({int(KP.n_events)} events).",
    "Endpoints: overall survival (OS) and progression-free survival (PFS), administratively censored at 5 years (60 months).",
    "Model: one Cox proportional-hazards model per feature and endpoint: Surv(time, event) ~ feature + FIGO stage (III/IV vs I/II) + residual disease "
    "(suboptimal vs optimal debulking) + age (years). 3 secretory measures x 2 endpoints = 6 models.",
    "Features: the three pre-specified secretory measures, pooled across the patient's final cores: SecA and SecB proportion of secretory (SecA + Intermediate + SecB) "
    "epithelial cells, and mean cell-level polarization UCell score across secretory cells. Cell-type densities are not modelled.",
    "Scaling: proportions and polarization are z-scored on their raw scale, using the mean/SD of the complete-case model cohort, so every feature HR is per +1 SD "
    "increase (SD_of_modelled_variable gives the SD). Covariate HRs are on their natural scale.",
    "FDR_BH: Benjamini-Hochberg adjustment across the 6 tests shown in Figure 7c,d (SecA proportion, SecB proportion and polarization; OS and PFS).",
    "KM_tertiles (Figure 7a,b): within each endpoint's model cohort, patients in the lowest polarization tertile (SecA-dominant, <= lower cut) are compared "
    "with those in the highest tertile (SecB-enriched, > upper cut); the middle tertile is left out, as for TCGA-OV (Figure 7g,h). "
    "Kaplan-Meier 5-year survival with 95% CI (log transform) and two-sided log-rank test (unadjusted).",
    "PH_test_*: Schoenfeld-residual test (survival::cox.zph, KM time transform) run on the same fitted models for this table (diagnostic; not shown in the figure). "
    f"Feature term p < 0.05 in {n_ph_viol} of {n_ph} models; global test p < 0.05 in {g_viol} of {n_ph} models. "
    "Offenders: " + "; ".join(f"{EP[e]} {label(f)} (p = {p:.3g})" for e, f, p in
                              ZPH[(ZPH.term == 'z') & (ZPH.p_value < 0.05)][['endpoint', 'feature', 'p_value']].itertuples(index=False)) + ".",
    "Univariable (unadjusted) Cox models: not part of the Figure 7 generator (all Cox estimates are adjusted) and are therefore not reported here; "
    "the only unadjusted comparison is the Figure 7a,b log-rank test (KM_tertiles).",
    "Sheets: Multivariable_Cox = feature HR per model with covariate HRs and PH tests (one row per model); Full_models = every coefficient of every model; "
    "PH_tests = full cox.zph output; Cohort = n/events/covariates per TMA for the complete-case cohort; KM_tertiles = Figure 7a,b.",
    "Source: model code from figures/figure7/01_fig7abcd_tma_cox_km_tables.R (Figure 7a-d plotted by figures/figure7/02_fig7abcd_tma_survival_row.R); the feature HR rows, "
    "the Figure 7c,d FDR and the KM rows are verified identical to figures/final_panels/figure7/Fig7_ABCD_cox_2core_results.csv, Fig7_ABCD_forest_fdr.csv and "
    "Fig7_ABCD_km_logrank.csv. Generated by tables/11_tableS11_xenium_tma_cox.py.",
    "Software: R survival::coxph (Efron ties), survfit and survdiff.",
]
write_table(11, "xenium_tma_multivariable_cox",
            title="Supplementary Table 11. Multivariable Cox regression of secretory polarization measures for 5-year OS and PFS in the Xenium HGSC TMA cohort (Figure 7a-d).",
            readme=readme,
            sheets=[("Multivariable_Cox", main), ("KM_tertiles", kmt), ("Full_models", fm), ("PH_tests", zp), ("Cohort", coh)])
