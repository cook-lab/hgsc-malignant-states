# Fig 7a-d model code (patients with two cores across TMA A + TMA B, multivariable adjusted Cox) on SingleR labels;
# Table S11 (tables/11_tableS11_xenium_tma_cox.py) embeds a copy of sections 1-4.
# Features: the three pre-specified secretory measures (SecA proportion, SecB proportion, polarization); the per-patient
# cell-type densities of section 1 are not modelled. Model = z(feature) + stage+debulk+age.
# Cohort: patients with >= 2 final TMA cores (spatial/02_qc/06b_final_tma_cores.py: per-array core QC on SingleR
# labels, primary HGSC, treatment-naive, >= 500 cells), without neoadjuvant chemotherapy or non-HGSC central histology.
# Writes cox_2core_results.csv (feature_type "adj"; 6 models = Fig 7c,d, q = BH across the 6) and
# tma_2core_patients.csv (Fig 7a,b KM) for figure7/02_fig7abcd_tma_survival_row.R.
suppressWarnings(suppressMessages({library(data.table); library(survival)}))
source("spatial/00_setup/final_tma_cores.R")   # run from the repo root
SG <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "singler_recompute")
SP <- file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "202609_final_hgsc_spatial/exploratory/survival_polarization_by_coregroup")
OUT<- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "figures/figure7")
SECA<-"SecA epithelium"; INT<-"Intermediate epithelium"; SECB<-"SecB epithelium"; SEC<-c(SECA,INT,SECB)

## ---- 1. Per-patient features on SingleR labels (densities pooled over the patient's cores) ----
feat_list<-list()
for (nm in c("tma","COEUR_TMA_G","COEUR_TMA_H")) {
  d <- fread(file.path(SG, paste0(nm,"_singler_cells.csv.gz")))
  d <- d[!is.na(cell_label_singler) & cell_label_singler!="" & sample_type=="tumour" &
         !is.na(core_id) & core_id!="Off core" & !is.na(patient_id)]
  d <- d[as.character(core_id) %in% final_core_ids(nm)]   # final TMA cores only; n_cores counts them
  # per-core area (mm^2), then pool per patient
  ca <- d[, .(area=sum(cell_area,na.rm=TRUE)/1e6), by=.(patient_id,core_id)]
  pa <- ca[, .(area=sum(area), n_cores=uniqueN(core_id)), by=patient_id]
  cnt<- d[, .N, by=.(patient_id,cell_label_singler)]
  cw <- dcast(cnt, patient_id~cell_label_singler, value.var="N", fill=0)
  f  <- merge(pa, cw, by="patient_id")
  ctcols <- setdiff(names(cw),"patient_id")
  for (ct in ctcols) set(f, j=paste0("dens_", gsub(" ","_",ct)), value=f[[ct]]/f$area)
  # secretory proportions + polar
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

## ---- 2. Clinical (TMA A: per_patient_features.csv; TMA B: coeur_clinical_perpatient_extracted.csv) ----
chum<-as.data.table(read.csv(file.path(Sys.getenv("DATA_ROOT", "/path/to/deposit"), "2026_final_xenium_analysis/output/10_clinical/per_patient_features.csv"), stringsAsFactors=FALSE))[, .(patient_id,os_time_5y,os_event_5y,pfs_time_5y,pfs_event_5y,stage,residual_binary,age)]
chum[, `:=`(patient_id=as.character(patient_id), tma="sfe_tma",
            os_time=os_time_5y, os_event=os_event_5y, pfs_time=pfs_time_5y, pfs_event=pfs_event_5y, age=as.numeric(age),
            stage_late=as.integer(suppressWarnings(as.numeric(stage))>=3),
            debulk=fifelse(residual_binary=="optimal","optimal",fifelse(residual_binary=="suboptimal","suboptimal",NA_character_)))]
coeur<-as.data.table(read.csv(file.path(SP,"coeur_clinical_perpatient_extracted.csv"), stringsAsFactors=FALSE))
coeur<-coeur[os_event %in% c(0,1) & !is.na(os_time)]
coeur[, patient_id:=as.character(patient)]   # raw id; the `tma` column disambiguates TMA A vs TMA B
coeur[, `:=`(os_event=as.integer(os_event==1 & as.numeric(os_time)<=60),
             pfs_event=as.integer(pfs_event==1 & suppressWarnings(as.numeric(pfs_time))<=60),
             os_time=pmin(as.numeric(os_time),60), pfs_time=pmin(suppressWarnings(as.numeric(pfs_time)),60),
             age=as.numeric(age), stage_late=as.integer(suppressWarnings(as.numeric(figo_code))>=3),
             debulk=fifelse(as.character(debulking)=="0","optimal",fifelse(as.character(debulking)=="1","suboptimal",NA_character_)))]
mkc<-function(x) x[, .(patient_id,tma,os_time,os_event,pfs_time,pfs_event,stage_late,debulk,age)]
cG<-copy(coeur)[,tma:="sfe_COEUR_TMA_G"]; cH<-copy(coeur)[,tma:="sfe_COEUR_TMA_H"]
clin<-rbind(mkc(chum),mkc(cG),mkc(cH))

## ---- 2b. Primary HGSC, treatment-naive at sampling ----
# Drops patients given neoadjuvant chemotherapy before the sampled tumour (TMA A treatment_status;
# TMA B Post_Chemo_Code) and TMA B cases whose central histology review is not HGSC
# (Histopathologyrevision 2-6); unreviewed histology is kept. Same sources as Table S8.
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
drop<-unique(merge(d2[,.(patient_id,tma)], ex, by=c("patient_id","tma")))
cat(sprintf("excluded (NACT or non-HGSC) among >=2-core patients: %d (%s)\n", uniqueN(drop[,.(patient_id,tma)]),
            paste(names(table(drop$why)), table(drop$why), collapse=", ")))
d2<-d2[!ex, on=.(patient_id, tma)]
# per-patient table for the Fig 7a-d survival row (figure7/02_fig7abcd_tma_survival_row.R); lab-internal (raw ids)
fwrite(d2[, .(patient_id, tma, n_cores, prop_SecA_of_sec, prop_SecB_of_sec, polar_mean,
              os_time, os_event, pfs_time, pfs_event, stage_late, debulk, age)], file.path(OUT, "tma_2core_patients.csv"))
cat(sprintf("2-core cohort: %d (CHUM %d, COEUR_G %d, COEUR_H %d)\n", nrow(d2),
  sum(d2$tma=="sfe_tma"),sum(d2$tma=="sfe_COEUR_TMA_G"),sum(d2$tma=="sfe_COEUR_TMA_H")))

## ---- 4. Multivariable Cox (adj: +stage+debulk+age), HR per +1 SD ----
FEATS<-c("prop_SecA_of_sec","prop_SecB_of_sec","polar_mean")   # pre-specified secretory measures only
zc<-function(x){s<-sd(x,na.rm=TRUE); if(is.na(s)||s==0) rep(NA,length(x)) else (x-mean(x,na.rm=TRUE))/s}
lbl<-function(f){ if(f=="prop_SecA_of_sec")return("SecA proportion"); if(f=="prop_SecB_of_sec")return("SecB proportion")
  if(f=="polar_mean")return("Polarization (mean UCell)"); gsub("_"," ",sub("dens_","",f)) |> (\(s) sub("^T cell$","T cell",s))() }
res<-list()
for (EPi in c("OS","PFS")){
  tcol<-if(EPi=="OS")"os_time" else "pfs_time"; ecol<-if(EPi=="OS")"os_event" else "pfs_event"
  for (fe in FEATS){
    sub<-copy(d2); covs<-c("stage_late","debulk","age")
    keep<-is.finite(sub[[fe]])&is.finite(sub[[tcol]])&!is.na(sub[[ecol]])&sub[[tcol]]>0
    for(cv in covs) keep<-keep & !is.na(sub[[cv]])
    sub<-sub[keep]; if(nrow(sub)<15||sum(sub[[ecol]])<8) next
    # log1p applies to density features only (none is modelled here); proportions and polarization enter raw
    val <- if(grepl("^dens_",fe)) log1p(sub[[fe]]) else sub[[fe]]
    sub[, z:=zc(val)]
    fit<-tryCatch(coxph(as.formula(paste0("Surv(",tcol,",",ecol,")~z+stage_late+debulk+age")),data=sub),error=function(e)NULL)
    if(is.null(fit)) next; s<-summary(fit)
    res[[length(res)+1]]<-data.table(endpoint=paste0(EPi,"_5yr"), feature=fe, feature_type="adj",
      n=s$n, n_events=s$nevent, HR=s$conf.int["z","exp(coef)"], HR_lower=s$conf.int["z","lower .95"],
      HR_upper=s$conf.int["z","upper .95"], p_value=s$coefficients["z","Pr(>|z|)"], feature_label=lbl(fe))
  }
}
R<-rbindlist(res); R[, q_value:=p.adjust(p_value,"BH")]   # BH across the 6 tests (3 measures x OS/PFS), as Fig 7c,d
setcolorder(R,c("endpoint","feature","feature_type","n","n_events","HR","HR_lower","HR_upper","p_value","q_value","feature_label"))
fwrite(R, file.path(OUT,"cox_2core_results.csv"))
cat(sprintf("OS n=%d ev=%d ; PFS n=%d ev=%d\n", R[endpoint=="OS_5yr"][1]$n, R[endpoint=="OS_5yr"][1]$n_events,
    R[endpoint=="PFS_5yr"][1]$n, R[endpoint=="PFS_5yr"][1]$n_events))
cat("headline (OS):\n"); print(R[endpoint=="OS_5yr" & feature %in% c("polar_mean","prop_SecB_of_sec","prop_SecA_of_sec"),
      .(feature, HR=round(HR,2), lo=round(HR_lower,2), hi=round(HR_upper,2), p=round(p_value,3))])
cat("FIG7AB_2CORE_DONE\n")
