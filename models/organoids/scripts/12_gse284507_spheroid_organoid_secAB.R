# Trimmed copy (2026-10-21) of 202603_organoids/scripts/12_gse284507_spheroid_organoid_secAB.R, the script that ran.
# Kept: 12_scores_all_sources.csv (read by 12h for Fig 3f-h) and the printed group medians. Removed: the spheroid vs
# organoid paired Wilcoxon contrast (not reported) and the saved R environment. The untrimmed script, which the
# header below describes, is in the lab project folder.

# 12 — SecB MARKER vs FUNCTION for GSE284507 (Exploiting cancer dormancy in EOC:
#   spheroid vs organoid; bulk RNA-seq, 7 iOvCa HGSC lines x {organoid,spheroid} x 2 reps = 28).
# Same matched-backbone method as 11u: atlas pseudobulk (SecA/Trans/SecB) + our PDO + DepMap +
# de la Puente cold co-culture + GSE284507, ALL in ONE matrix, one CPM, ONE UCell call at
# maxRank=FULL common gene space. UCell is rank-based per-sample so linear/log/CPM don't move it.
# Biology tested: EOC SPHEROIDS = the dormant, quiescent, matrix-detached peritoneal-metastasis
# state (paper: AMPK up, Akt down, cell-cycle down); ORGANOIDS = proliferative. Do spheroids
# shift toward SecB FUNCTION (HIF-induction - proliferation, i.e. quiescence) while both keep the
# SecB MARKER (keratins/SLPI/TACSTD2)?  Patient-matched within-line paired contrast (7 lines).
suppressPackageStartupMessages({library(Seurat);library(UCell);library(data.table);library(Matrix);library(dplyr);library(readxl);library(org.Hs.eg.db)})
set.seed(42)
proj<-"/Volumes/CookLab/Sarah/202603_organoids"; dd<-file.path(proj,"data/11_nextgen")
out<-file.path(proj,"output/12_GSE284507_spheroid_organoid"); dir.create(out,showWarnings=FALSE,recursive=TRUE)
atlasdir<-"/Volumes/CookLab/Sarah/202601_final_atlas/output/18_ucell_atlas"
shared<-readLines(file.path(atlasdir,"shared_gene_list.txt"))

# canonical modules — IDENTICAL to 11q/11u
MOD<-list(
  SecB_marker  = c("KRT17","KRT19","KRT7","LCN2","PRSS22","SLPI","TACSTD2"),
  Hypoxia_HIF  = c("NDRG1","SLC16A3","VEGFA","CA9","BNIP3","P4HA1","PDK1","SLC2A1"),
  Proliferation= c("MKI67","TOP2A","PCNA","CCNB1","CDK1"),
  SecA         = c("MECOM","FBXO21","LGR5","LPAR3","PBX1","SOX17","RCN2"))

# --- atlas pseudobulk (raw UMI) ---
apb<-fread(file.path(dd,"atlas_pseudobulk_counts.csv")); A<-as.matrix(apb[,-1]); rownames(A)<-apb$gene
ameta<-fread(file.path(dd,"atlas_pseudobulk_meta.csv"))
# --- our PDO pseudobulk (raw UMI) ---
o<-readRDS(file.path(proj,"seurat_untreated_baseline.rds")); o<-JoinLayers(o)
cts<-GetAssayData(o,assay="RNA",layer="counts"); mdl<-o$pdo_model
O<-sapply(sort(unique(mdl)),function(m) Matrix::rowSums(cts[,mdl==m,drop=FALSE])); rownames(O)<-rownames(cts)
# --- DepMap next-gen organoids (un-log to TPM) ---
de<-fread(file.path(dd,"next_gen_expression.csv")); setnames(de,1,"ModelID")
D<-t(as.matrix(de[,-1])); rownames(D)<-sub(" \\(.*","",rownames(D)); D<-D[!duplicated(rownames(D)),]; D<-2^D-1; colnames(D)<-de$ModelID
# --- de la Puente cold co-culture epithelial -> pseudobulk per condition ---
dp<-readRDS(file.path(proj,"output/10_HGSC_public_SecAB_scoring/10a_seurat_epithelial.rds")); dp<-JoinLayers(dp)
dc<-GetAssayData(dp,assay="RNA",layer="counts"); cond<-as.character(dp$condition); conds<-sort(unique(cond))
P<-sapply(conds,function(cc) Matrix::rowSums(dc[,cond==cc,drop=FALSE])); rownames(P)<-rownames(dc); colnames(P)<-paste0("delaPuente_",conds)

# --- GSE284507 (this study): normalized bulk counts, ENTREZID-keyed -> collapse to symbol ---
g<-as.data.frame(read_excel(file.path(proj,"data/GSE284507/GSE284507_normalized_counts.xlsx")))
colnames(g)<-sub("_combined: Counts","",colnames(g))            # e.g. "182-ORG-1"
ez<-as.character(g$ENTREZID); G<-as.matrix(g[,-1])
sym<-mapIds(org.Hs.eg.db,keys=ez,column="SYMBOL",keytype="ENTREZID",multiVals="first")
keep<-!is.na(sym); G<-G[keep,,drop=FALSE]; sym<-sym[keep]
G<-rowsum(G,group=sym)                                          # collapse duplicate symbols (sum)
cat(sprintf("GSE284507: %d entrez -> %d symbols x %d samples\n",length(ez),nrow(G),ncol(G)))
gsm<-data.frame(sample=colnames(G),
                line=sub("-.*","",colnames(G)),
                state=ifelse(grepl("-ORG",colnames(G)),"organoid","spheroid"),
                stringsAsFactors=FALSE)

# --- ONE matched matrix, common gene space, CPM, ONE UCell call at maxRank=full ---
common<-Reduce(intersect,list(shared,rownames(A),rownames(O),rownames(D),rownames(P),rownames(G)))
comb<-cbind(A[common,],O[common,],D[common,],P[common,],G[common,])
grp<-c(rep("atlas",ncol(A)),rep("our",ncol(O)),rep("depmap",ncol(D)),rep("delapuente",ncol(P)),rep("gse284507",ncol(G)))
comb<-sweep(comb,2,colSums(comb),"/")*1e6
for(m in names(MOD)){mi<-intersect(MOD[[m]],rownames(comb)); cat(sprintf("  %-14s %d/%d genes present\n",m,length(mi),length(MOD[[m]]))); MOD[[m]]<-mi}
sc<-ScoreSignatures_UCell(comb,features=MOD,maxRank=length(common),w_neg=1,ties.method="average")
res<-as.data.frame(sc); colnames(res)<-sub("_UCell","",colnames(res)); res$sample<-colnames(comb); res$grp<-grp
res$sub<-NA
res$sub[grp=="atlas"]<-ameta$celltype_nmf[match(res$sample[grp=="atlas"],ameta$group)]
meta<-fread(file.path(dd,"depmap_model_metadata.csv"))
res$sub[grp=="depmap"]<-meta$OncotreeSubtype[match(res$sample[grp=="depmap"],meta$ModelID)]
res$sub[grp=="delapuente"]<-sub("delaPuente_","",res$sample[grp=="delapuente"])
res$sub[grp=="gse284507"]<-gsm$state[match(res$sample[grp=="gse284507"],gsm$sample)]
res$line<-NA; res$line[grp=="gse284507"]<-gsm$line[match(res$sample[grp=="gse284507"],gsm$sample)]
res$hgsc<-grepl("High-Grade Serous",res$sub)
# Function axis (z-scored across ALL samples, as 11s): HIF induction minus proliferation
z<-function(x)(x-mean(x))/(sd(x)+1e-9)
res$Function<-z(res$Hypoxia_HIF)-z(res$Proliferation)
write.csv(res,file.path(out,"12_scores_all_sources.csv"),row.names=FALSE)
cat(sprintf("\ncommon genes=%d | total samples=%d\n",length(common),ncol(comb)))

# ---------- reference medians (atlas backbone) ----------
amed<-function(s,c) median(res[[c]][res$grp=="atlas"&res$sub==s])
cat("\n=== ATLAS BACKBONE (reference) ===\n")
for(s in c("SecA epithelium","Transitioning epithelium","SecB epithelium"))
  cat(sprintf("  %-24s SecB_marker=%.3f Hypoxia_HIF=%.3f Prolif=%.3f SecA=%.3f Function=%+.2f\n",
      s,amed(s,"SecB_marker"),amed(s,"Hypoxia_HIF"),amed(s,"Proliferation"),amed(s,"SecA"),amed(s,"Function")))

# ---------- GSE284507 spheroid vs organoid ----------
gr<-res[res$grp=="gse284507",]
grpmed<-function(st) {m<-gr$sub==st; cat(sprintf("  %-9s (n=%d) SecB_marker=%.3f Hypoxia_HIF=%.3f Prolif=%.3f SecA=%.3f Function=%+.2f\n",
      st,sum(m),median(gr$SecB_marker[m]),median(gr$Hypoxia_HIF[m]),median(gr$Proliferation[m]),median(gr$SecA[m]),median(gr$Function[m]))) }
cat("\n=== GSE284507 medians by state ===\n"); grpmed("organoid"); grpmed("spheroid")
