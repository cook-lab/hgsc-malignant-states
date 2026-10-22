#!/usr/bin/env Rscript
# ============================================================================
# Fig 5j -- MACROPHAGE MORPHOLOGY vs the local epithelial polarization niche
# ----------------------------------------------------------------------------
#   * cell identity = SingleR cell_label_singler (spatial/00_setup/singler_labels.R)
#   * niche = mean polarization_UCell of SECRETORY epithelial neighbours
#     (SecA / Intermediate / SecB) within 50 um of the macrophage (frNN on
#     spatialCoords); macrophages with < 5 epithelial neighbours dropped
#   * per-unit tertiles of niche polarization; unit = tissue (WT) or patient (TMA)
#   * TMA units: the final TMA core set (spatial/00_setup/final_tma_cores.R: per-array
#     core QC on SingleR labels, primary HGSC, treatment-naive, >= 500 cells), patients keyed by TMA
#     ("<tma>|<patient_id>": TMA A and TMA B reuse some patient numbers), and the >=500-cell rule
#     counted per TMA patient on the retained cores
#   * morphometrics from cellSeg / nucSeg polygons, matched to barcodes BY NAME
#     (sfe_tma carries polygons for only 309,797 / 671,453 cells)
#   * 11 features (cell and nucleus area, perimeter, circularity, eccentricity; N:C ratio,
#     cell solidity, nucleus-cell centroid offset); implausible values filtered (cell area
#     10-5,000 um^2, nucleus 5-1,000 um^2, N:C <= 1). Fig 5j shows cell_area.
# Stats: per-unit median in low vs high tertile, paired Wilcoxon, Bonferroni x11.
# Eligibility: >= 10 macrophages in both the low and high tertile.
# OUTPUTS (output_root/split_recompute/; folder and file names are historical):
#   fig5e_mac_morphology_percell_presplit.csv.gz, fig5e_mac_morphology_stats_presplit.csv,
#   fig5e_mac_morphology_unitmedians_presplit.csv (drawn and tested by 11_fig5_compose_panels.py)
# MANUSCRIPT PANEL(S): Fig 5j. RUNTIME TIER: heavy (polygon geometry of 11 SFEs).
# ============================================================================
source("spatial/00_setup/00_setup.R")
source("spatial/00_setup/singler_labels.R")   # singler_labels(): SingleR labels
source("spatial/00_setup/final_tma_cores.R")
suppressWarnings(suppressMessages({
  library(SpatialFeatureExperiment); library(SummarizedExperiment)
  library(sf); library(dbscan); library(data.table)
}))
set.seed(CFG$seed)
OUT <- cfg_path("output_root","split_recompute")
RADIUS <- 50; MIN_EPI_NB <- 5; MIN_PER_TERTILE <- 10; MIN_PATIENT_CELLS <- 500
EPI <- c("SecA epithelium","Intermediate epithelium","SecB epithelium")
WT  <- paste0("sfe_", CFG$cohort$whole_tissue)
TMA <- c("sfe_tma","sfe_COEUR_TMA_G","sfe_COEUR_TMA_H")
FEATURES <- c("cell_area","nuc_area","nc_ratio","cell_perimeter","nuc_perimeter",
              "cell_circularity","nuc_circularity","cell_solidity",
              "cell_eccentricity","nuc_eccentricity","nc_centroid_offset")

ecc <- function(g) vapply(seq_along(g), function(i) {
  cc <- sf::st_coordinates(g[i]); if (nrow(cc) < 3) return(NA_real_)
  m <- cov(cc[,1:2]); if (any(!is.finite(m))) return(NA_real_)
  ev <- pmax(eigen(m, symmetric=TRUE, only.values=TRUE)$values, 0)
  if (max(ev) <= 0) NA_real_ else 1 - sqrt(min(ev)/max(ev)) }, numeric(1))

res <- list(); ptot <- list()
for (nm in c(WT, TMA)) {
  t0 <- Sys.time(); message("[", format(t0), "] ", nm)
  sfe <- load_sfe(nm)
  L <- singler_labels(sfe, nm); lab <- L$lab; pol <- L$polarization
  bc <- colnames(sfe); xy <- spatialCoords(sfe)
  unit <- if (nm %in% WT) rep(nm, length(bc)) else {
    s <- data.table::fread(file.path(cfg_path("output_root"),"singler_recompute",
           paste0(sub("^sfe_","",nm),"_singler_cells.csv.gz")),
           select=c("barcode","patient_id","core_id","sample_type"))
    s[, core := as.character(core_id)]
    s <- s[!is.na(patient_id) & !is.na(core) & core %in% final_core_ids(nm)]   # final TMA cores
    s[, unit := paste0(sub("^sfe_","",nm), "|", as.character(patient_id))]
    ptot[[nm]] <- s[, .N, by = unit]                      # all cell types on retained cores
    s$unit[match(bc, s$barcode)] }
  cs <- colGeometry(sfe,"cellSeg"); ns <- colGeometry(sfe,"nucSeg")
  stopifnot(identical(rownames(cs), rownames(ns)))
  rm(sfe); gc(verbose=FALSE)

  is_epi <- lab %in% EPI & is.finite(pol)
  mac <- which(lab == "Macrophage" & bc %in% rownames(cs))
  nn <- dbscan::frNN(x = xy[is_epi,,drop=FALSE], eps = RADIUS, query = xy[mac,,drop=FALSE], sort = FALSE)
  pe <- pol[is_epi]
  n_epi <- lengths(nn$id)
  niche <- vapply(nn$id, function(i) if (length(i)) mean(pe[i]) else NA_real_, numeric(1))
  ok <- n_epi >= MIN_EPI_NB & !is.na(unit[mac])
  mac <- mac[ok]; n_epi <- n_epi[ok]; niche <- niche[ok]
  message(sprintf("   macrophages with polygons + >=%d epi neighbours: %d", MIN_EPI_NB, length(mac)))

  gi <- match(bc[mac], rownames(cs))
  mc <- sf::st_make_valid(sf::st_geometry(cs)[gi]); mn <- sf::st_make_valid(sf::st_geometry(ns)[gi])
  v <- sf::st_is_valid(mc) & sf::st_is_valid(mn) &
       sf::st_geometry_type(mc) %in% c("POLYGON","MULTIPOLYGON") &
       sf::st_geometry_type(mn) %in% c("POLYGON","MULTIPOLYGON")
  mc <- mc[v]; mn <- mn[v]; mac <- mac[v]; n_epi <- n_epi[v]; niche <- niche[v]
  ca <- as.numeric(sf::st_area(mc)); na <- as.numeric(sf::st_area(mn))
  cp <- as.numeric(sf::st_length(sf::st_boundary(mc))); np <- as.numeric(sf::st_length(sf::st_boundary(mn)))
  dt <- data.table(cohort = ifelse(nm %in% WT,"WT","TMA"), sample = nm, unit = unit[mac],
                   barcode = bc[mac], n_epi_nb = n_epi, niche_pol = niche,
                   cell_area = ca, nuc_area = na, nc_ratio = na/pmax(ca,1e-6),
                   cell_perimeter = cp, nuc_perimeter = np,
                   cell_circularity = ifelse(cp>0, 4*pi*ca/cp^2, NA_real_),
                   nuc_circularity  = ifelse(np>0, 4*pi*na/np^2, NA_real_),
                   cell_solidity = ca/pmax(as.numeric(sf::st_area(sf::st_convex_hull(mc))),1e-6),
                   cell_eccentricity = ecc(mc), nuc_eccentricity = ecc(mn),
                   nc_centroid_offset = as.numeric(sf::st_distance(sf::st_centroid(mc), sf::st_centroid(mn), by_element=TRUE)))
  n0 <- nrow(dt)
  dt <- dt[is.finite(cell_area) & is.finite(nuc_area) & cell_area >= 10 & cell_area <= 5000 &
           nuc_area >= 5 & nuc_area <= 1000 & nc_ratio <= 1 &
           is.finite(cell_perimeter) & is.finite(nuc_perimeter)]
  message(sprintf("   filter kept %d / %d ; %.1f min", nrow(dt), n0,
                  as.numeric(difftime(Sys.time(), t0, units="mins"))))
  res[[nm]] <- dt; rm(cs, ns, mc, mn); gc(verbose=FALSE)
}
all <- rbindlist(res)
# cohort rule: >= 500 cells per patient, counted over ALL cells on the patient's
# retained tumour cores, per TMA (as 03_fig5b_nuclear_area.R)
ptot <- rbindlist(ptot)
keep_pt <- ptot[N >= MIN_PATIENT_CELLS, unit]
message(sprintf("TMA patients with >=%d cells: %d / %d", MIN_PATIENT_CELLS, length(keep_pt), nrow(ptot)))
all <- all[cohort == "WT" | unit %in% keep_pt]
fwrite(all, file.path(OUT, "fig5e_mac_morphology_percell_presplit.csv.gz"))
# rank-based tertiles within unit (quantile breaks are non-unique for TMA
# patients with 1-2 macrophages; those fail the >=10-per-tertile rule anyway)
all[, tertile := c("low","mid","high")[pmin(3L, ceiling(3 * frank(niche_pol, ties.method="first") / .N))],
    by = .(cohort, unit)]

stat_rows <- list(); med_rows <- list()
for (co in c("WT","TMA")) {
  d <- all[cohort == co & tertile %in% c("low","high")]
  nt <- dcast(d[, .N, by=.(unit,tertile)], unit ~ tertile, value.var="N", fill=0)
  elig <- nt[low >= MIN_PER_TERTILE & high >= MIN_PER_TERTILE, unit]
  message(sprintf("%s eligible units: %d / %d", co, length(elig), nrow(nt)))
  med <- d[unit %in% elig, lapply(.SD, median, na.rm=TRUE), by=.(unit,tertile), .SDcols=FEATURES]
  med_rows[[co]] <- cbind(cohort=co, med)
  for (f in FEATURES) {
    w <- dcast(med, unit ~ tertile, value.var=f)[is.finite(low) & is.finite(high)]
    dl <- w$high - w$low; dir <- ifelse(median(dl) > 0, "up", "down")
    p <- wilcox.test(w$high, w$low, paired=TRUE, exact = nrow(w) <= 50)$p.value
    stat_rows[[paste(co,f)]] <- data.table(cohort=co, feature=f, n_units=nrow(w), direction=dir,
      n_concordant = if (dir=="up") sum(dl>0) else sum(dl<0),
      median_low = median(w$low), median_high = median(w$high),
      median_delta = median(dl), p = p, p_bonf = min(1, p*length(FEATURES)))
  }
}
st <- rbindlist(stat_rows)
fwrite(st, file.path(OUT, "fig5e_mac_morphology_stats_presplit.csv"))
fwrite(rbindlist(med_rows), file.path(OUT, "fig5e_mac_morphology_unitmedians_presplit.csv"))
print(st[, .(cohort, feature, n_units, direction, n_concordant, median_low=signif(median_low,4),
             median_high=signif(median_high,4), p=signif(p,3), p_bonf=signif(p_bonf,3))])
message("DONE")
