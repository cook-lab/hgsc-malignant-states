#!/usr/bin/env Rscript
# PURPOSE: Fig 6 per-section Xenium effects, BOTH partners, same sections (8 whole tissues, SingleR labels):
#   mac_niche_vs_ctrl : niche macrophages vs control macrophages (mac side; 6f bottom row, 6h macrophage rows)
#   secb_contact      : SecB <=30 um of a niche mac vs SecB with no macrophage within 30 um (6h contact row)
#   secb_vs_seca      : SecB vs SecA epithelium (6h SecB rows)
# Genes: log2FC of mean CP10k. Pathway sets: difference in mean set score (mean z of log1p CP10k over set genes).
# Summary across sections: n_up/n, Wilcoxon signed-rank (groups need >= 20 cells per section).
# INPUTS: singler_recompute cell tables; fig6_viewer/degdata count export; mac_niche_WT.csv.gz (04);
#   <fig6_singler>/inputs/pathway_gene_sets_v2.csv.
# OUTPUTS (<fig6_singler>/): reciprocal_persample.csv (per-section record), reciprocal_summary.csv (read by 11 and 13).
# MANUSCRIPT PANEL(S): Fig 6f, 6h; the contact statements of the Results. RUNTIME TIER: moderate.
suppressWarnings(suppressMessages({library(data.table); library(Matrix); library(dbscan)}))
SG<-file.path(Sys.getenv("OUTPUT_ROOT", "output"), "singler_recompute"); DEG<-file.path(Sys.getenv("OUTPUT_ROOT", "output"), "fig6_viewer/degdata")
OUT<-file.path(Sys.getenv("OUTPUT_ROOT", "output"), "fig6_singler")
WT<-c("WT01","WT02","WT03","WT04","WT05","WT06","WT07","WT08")
PS<-fread(file.path(OUT,"inputs/pathway_gene_sets_v2.csv"))
SETS<-list(nfkb="nfkb", survival="survival_antiapoptotic", lactate="lactate_metabolism", hypoxia="hypoxia",
           complement_reg="complement_regulation", checkpoint="checkpoint_inhibitory")
SETS<-lapply(SETS, function(p) unique(PS[grepl(paste0("^",p,"$"),pathway), gene]))
GENES<-c("TNF","NFKBIA","NFKB1","NFKB2","RELA","BIRC3","BCL2L1","BCL2","VEGFA","CD55","CD274","LGALS9",
         "SLC16A3","LDHA","SLC2A1","ICAM1","CD44","MMP7","ITGB6","VCAM1","C7","CTSL","MMP9","MCL1","MRC1","C1QC","TREM2","LGALS3","IL32","TGFB1","INHBA","CD83","MXD1","TGM2","NDRG1","PGK1","ENO1","ALDOA","PDK1","ITGAV","A2M","HAVCR2","LAG3")
R_CT<-30; niche<-fread(file.path(OUT,"mac_niche_WT.csv.gz")); res<-list()
for (s in WT) {
  d<-as.data.table(read.csv(gzfile(file.path(SG,paste0(s,"_singler_cells.csv.gz")))))[,.(barcode,x,y,cell_label_singler)]
  d<-d[!is.na(cell_label_singler)&cell_label_singler!=""]; d[niche[sample==s],on="barcode",grp:=i.group]
  m<-readMM(gzfile(file.path(DEG,paste0(s,"_counts.mtx.gz")))); g<-readLines(file.path(DEG,paste0(s,"_genes.txt")))
  cb<-fread(file.path(DEG,paste0(s,"_cells.csv.gz")),select="barcode")$barcode; lib<-pmax(Matrix::colSums(m),1)
  need<-intersect(unique(c(GENES,unlist(SETS))),g); gi<-match(need,g)
  CP<-as.matrix(t(m[gi,,drop=FALSE])); colnames(CP)<-need; CP<-sweep(CP,1,lib,"/")*1e4; rownames(CP)<-cb; rm(m); gc(verbose=FALSE)
  L<-log1p(CP); Z<-scale(L); Z[is.na(Z)]<-0
  SS<-sapply(SETS,function(gs){gs<-intersect(gs,colnames(Z)); if(length(gs)) rowMeans(Z[,gs,drop=FALSE]) else rep(NA_real_,nrow(Z))})
  d<-d[barcode %in% cb]; xy<-as.matrix(d[,.(x,y)])
  mac_all<-which(d$cell_label_singler=="Macrophage"); ni<-which(d$grp=="niche")
  secb<-which(d$cell_label_singler=="SecB epithelium"); seca<-which(d$cell_label_singler=="SecA epithelium")
  near<-lengths(frNN(x=xy[ni,,drop=FALSE],eps=R_CT,query=xy[secb,,drop=FALSE])$id)>0
  anym<-lengths(frNN(x=xy[mac_all,,drop=FALSE],eps=R_CT,query=xy[secb,,drop=FALSE])$id)>0
  CMP<-list(mac_niche_vs_ctrl=list(d$barcode[d$grp %in% "niche"], d$barcode[d$grp %in% "ctrl"]),
            secb_contact=list(d$barcode[secb][near], d$barcode[secb][!anym]),
            secb_vs_seca=list(d$barcode[secb], d$barcode[seca]))
  for (cn in names(CMP)) { a<-CMP[[cn]][[1]]; b<-CMP[[cn]][[2]]
    if (length(a)<20||length(b)<20) next
    for (gg in intersect(GENES,colnames(CP)))
      res[[length(res)+1]]<-data.table(sample=s,comparison=cn,feature=gg,type="gene",n_a=length(a),n_b=length(b),
        effect=log2(mean(CP[a,gg])+0.01)-log2(mean(CP[b,gg])+0.01))
    for (st in names(SETS)) if (!all(is.na(SS[,st])))
      res[[length(res)+1]]<-data.table(sample=s,comparison=cn,feature=paste0("set:",st),type="set",n_a=length(a),n_b=length(b),
        effect=mean(SS[a,st])-mean(SS[b,st]))
  }
  message(s," done")
}
X<-rbindlist(res); fwrite(X,file.path(OUT,"reciprocal_persample.csv"))
S<-X[,.(n_sec=.N,n_up=sum(effect>0),mean_effect=mean(effect),
        p=if(.N>=3) suppressWarnings(wilcox.test(effect,mu=0)$p.value) else NA_real_),by=.(comparison,feature,type)]
fwrite(S,file.path(OUT,"reciprocal_summary.csv"))
W<-dcast(S,feature~comparison,value.var=c("mean_effect","n_up","n_sec","p"))
print(W[,.(feature, mac=sprintf("%+.2f %d/%d p=%.3f",mean_effect_mac_niche_vs_ctrl,n_up_mac_niche_vs_ctrl,n_sec_mac_niche_vs_ctrl,p_mac_niche_vs_ctrl),
             secb_contact=sprintf("%+.2f %d/%d p=%.3f",mean_effect_secb_contact,n_up_secb_contact,n_sec_secb_contact,p_secb_contact),
             secb_vs_seca=sprintf("%+.2f %d/%d p=%.3f",mean_effect_secb_vs_seca,n_up_secb_vs_seca,n_sec_secb_vs_seca,p_secb_vs_seca))])
cat("RECIP_DONE\n")
