#!/usr/bin/env Rscript
# PURPOSE: Fig 6 SecB-niche macrophage definition (SingleR labels) + pseudobulk niche-vs-control DEG, 8 whole tissues.
# Niche mac  = Macrophage with >=5 SecB epithelial cells within 50 um (whole section).
# Control mac= Macrophage with 0 SecB and >=5 SecA within 50 um.
# DEG: pseudobulk counts per section x group, limma-voom ~ section + group (paired), sections with >= 50 niche and
# >= 50 control macrophages (7 of 8). Counts = raw 477-gene panel (fig6_viewer/degdata count export).
# OUTPUTS (<fig6_singler>/): mac_niche_WT.csv.gz (group per macrophage; read by 05, 07 and fig6_roi_picker/roi_server.py),
#   mac_counts_WT.csv and mac_signature_deg_WT.csv (Table S10, tables/10_tableS10_secb_niche_macrophage_degs.py).
# The by_core branches of classify() / run_cohort() are not used by the whole-tissue call.
# MANUSCRIPT PANEL(S): Fig 6f, 6h, 6i (niche labels); Table S10. RUNTIME TIER: moderate.
suppressWarnings(suppressMessages({library(data.table); library(Matrix); library(dbscan); library(limma); library(edgeR)}))
SG  <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "singler_recompute")
DEG <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "fig6_viewer/degdata")
OUT <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "fig6_singler"); dir.create(OUT, showWarnings=FALSE)
R <- 50; MIN_SECB <- 5; MIN_SECA <- 5; MING <- 50; MINCORE <- 500
WT  <- c("WT01","WT02","WT03","WT04","WT05","WT06","WT07","WT08")

count_nb <- function(q, ref, eps){  # number of ref points within eps of each query point
  if (!nrow(ref) || !nrow(q)) return(integer(nrow(q)))
  lengths(frNN(x=as.matrix(ref), eps=eps, query=as.matrix(q))$id)
}
classify <- function(d, by_core){
  mac <- d[cell_label_singler=="Macrophage"]
  secb <- d[cell_label_singler=="SecB epithelium"]; seca <- d[cell_label_singler=="SecA epithelium"]
  if (!by_core) {
    mac[, nSecB := count_nb(mac[,.(x,y)], secb[,.(x,y)], R)]
    mac[, nSecA := count_nb(mac[,.(x,y)], seca[,.(x,y)], R)]
  } else {
    mac[, `:=`(nSecB=0L, nSecA=0L)]
    for (cc in unique(mac$core_id)) {
      i <- which(mac$core_id==cc)
      mac[i, nSecB := count_nb(mac[i,.(x,y)], secb[core_id==cc,.(x,y)], R)]
      mac[i, nSecA := count_nb(mac[i,.(x,y)], seca[core_id==cc,.(x,y)], R)]
    }
  }
  mac[, group := fifelse(nSecB>=MIN_SECB, "niche", fifelse(nSecB==0 & nSecA>=MIN_SECA, "ctrl", "other"))]
  mac
}
load_counts <- function(s){
  m <- readMM(gzfile(file.path(DEG, paste0(s,"_counts.mtx.gz"))))       # genes x cells
  rownames(m) <- readLines(file.path(DEG, paste0(s,"_genes.txt")))
  colnames(m) <- fread(file.path(DEG, paste0(s,"_cells.csv.gz")), select="barcode")$barcode
  as(m, "CsparseMatrix")
}

run_cohort <- function(samples, by_core, tag){
  macs <- list(); pb <- list()
  for (s in samples) {
    t0 <- Sys.time()
    d <- fread(file.path(SG, paste0(s,"_singler_cells.csv.gz")))
    d <- d[!is.na(cell_label_singler) & cell_label_singler!=""]
    if (by_core) {
      d <- d[sample_type=="tumour" & !is.na(core_id) & core_id!="Off core"]
      keep_core <- d[, .N, by=core_id][N>=MINCORE, core_id]; d <- d[core_id %in% keep_core]
    }
    mac <- classify(d, by_core); mac[, sample:=s]
    macs[[s]] <- mac[, .(sample, barcode, x, y, core_id=if(by_core) core_id else NA_character_, patient_id=if(by_core) patient_id else NA_character_, nSecB, nSecA, group)]
    cnt <- load_counts(s)
    for (g in c("niche","ctrl")) {
      bc <- intersect(mac[group==g, barcode], colnames(cnt))
      if (length(bc)) pb[[paste(s,g,sep="|")]] <- data.table(sample=s, group=g, n=length(bc), gene=rownames(cnt), count=as.numeric(Matrix::rowSums(cnt[, bc, drop=FALSE])))
    }
    message(sprintf("[%s] %-12s macs %6d  niche %5d  ctrl %5d  (%.1f min)", tag, s, nrow(mac),
                    sum(mac$group=="niche"), sum(mac$group=="ctrl"), as.numeric(difftime(Sys.time(),t0,units="mins"))))
    rm(cnt, d); gc(verbose=FALSE)
  }
  M <- rbindlist(macs); fwrite(M, file.path(OUT, paste0("mac_niche_",tag,".csv.gz")))
  P <- rbindlist(pb)
  nn <- unique(P[, .(sample, group, n)]); fwrite(nn, file.path(OUT, paste0("mac_counts_",tag,".csv")))
  ok <- nn[group=="niche" & n>=MING, sample]; ok <- intersect(ok, nn[group=="ctrl" & n>=MING, sample])
  message(sprintf("[%s] sections passing >=%d niche & ctrl macs: %s", tag, MING, paste(ok, collapse=", ")))
  W <- dcast(P[sample %in% ok], gene ~ sample + group, value.var="count", fill=0)
  mat <- as.matrix(W[, -1]); rownames(mat) <- W$gene
  cols <- colnames(mat); smp <- factor(sub("_(niche|ctrl)$","",cols)); grp <- factor(sub("^.*_","",cols), levels=c("ctrl","niche"))
  y <- DGEList(mat); keep <- filterByExpr(y, group=grp); y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE])
  des <- model.matrix(~ smp + grp); v <- voom(y, des); fit <- eBayes(lmFit(v, des))
  tt <- as.data.table(topTable(fit, coef="grpniche", number=Inf), keep.rownames="gene")
  setnames(tt, c("adj.P.Val"), c("FDR"))
  # per-section log2FC (CPM) + consistency
  cpm <- edgeR::cpm(y, log=FALSE)
  per <- sapply(levels(smp), function(s) log2((cpm[, paste0(s,"_niche")]+1)/(cpm[, paste0(s,"_ctrl")]+1)))
  per <- as.data.table(per, keep.rownames="gene")
  tt <- merge(tt, per, by="gene")
  sc <- setdiff(names(per),"gene")
  tt[, n_up := rowSums(.SD > 0.1), .SDcols=sc][, n_dn := rowSums(.SD < -0.1), .SDcols=sc][, n_sec := length(sc)]
  setorder(tt, P.Value)
  fwrite(tt, file.path(OUT, paste0("mac_signature_deg_",tag,".csv")))
  message(sprintf("[%s] genes tested %d, FDR<0.05: %d", tag, nrow(tt), sum(tt$FDR<0.05)))
  tt
}
wt  <- run_cohort(WT,  FALSE, "WT")
cat("\nFIG6_NICHE_DONE\n")
