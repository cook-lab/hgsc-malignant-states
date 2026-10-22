# One-time export per sample for the interactive ROI-DEG server.
# Writes to fig6_viewer/degdata/<sample>_{counts.mtx.gz, cells.csv.gz, genes.txt}
# Counts = RAW SFE counts (477 genes incl MMP7/LGALS3). Coords = cellSeg centroid (matches viewer .json extent frame).
# Labels = base cell_label with SPLIT epithelial overlay (SecA/Int/SecB). usage: deg_precompute.R <sfe_name>
# run with cwd = hgsc-malignant-states repo (sources its spatial/00_setup + config)
.here<-"figures/figure4"; source(file.path("spatial","00_setup","00_setup.R"))
suppressPackageStartupMessages({library(SpatialFeatureExperiment); library(Matrix); library(data.table); library(sf)})
args<-commandArgs(trailingOnly=TRUE); nm<-args[1]
sid<-if(length(args)>=2 && nzchar(args[2])) args[2] else NA_character_   # optional colData$sample_id filter (split TMAs)
OUT<-cfg_path("output_root","fig6_viewer","degdata"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
s<-sub("^sfe_","",nm)
ob<-if(length(args)>=3 && nzchar(args[3])) args[3] else if(!is.na(sid)) paste0(s,"_",sid) else s   # output basename
sc<-fread(cfg_path("output_root","split_recompute","purified_polarization_scores.csv"),
          select=c("sample","barcode","cell_label_split"))
message("[",Sys.time(),"] ",nm," sid=",sid," -> ",ob)
sfe<-load_sfe(nm)
lab<-as.character(sfe$cell_label); lab[lab=="Transitioning epithelium"]<-"Intermediate epithelium"
sub<-sc[sample==nm]; m<-match(colnames(sfe),sub$barcode); lab[!is.na(m)]<-sub$cell_label_split[m[!is.na(m)]]
keep<- !is.na(lab) & lab!=""
if(!is.na(sid)) keep<-keep & (as.character(sfe$sample_id)==sid)
sfe<-sfe[,keep]; lab<-lab[keep]
xy<-st_coordinates(st_centroid(st_geometry(colGeometry(sfe,"cellSeg"))))
cnt<-as(counts(sfe),"CsparseMatrix")   # genes x cells, raw integer
genes<-rownames(cnt)
# write
mtx<-file.path(OUT,paste0(ob,"_counts.mtx"))
writeMM(cnt, mtx); system2("gzip",c("-f",shQuote(mtx)))
cells<-data.table(barcode=colnames(cnt), x=xy[,1], y=xy[,2], cell_label=lab)
fwrite(cells, file.path(OUT,paste0(ob,"_cells.csv.gz")))
writeLines(genes, file.path(OUT,paste0(ob,"_genes.txt")))
message("  ",ob,": ",ncol(cnt)," cells x ",nrow(cnt)," genes")
cat("DEGPREP_DONE ",ob,"\n")
