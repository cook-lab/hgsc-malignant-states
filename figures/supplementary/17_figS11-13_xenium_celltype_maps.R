#!/usr/bin/env Rscript
# =============================================================================
# Figures S11 / S12 / S13 — Xenium cell-type spatial maps on the final SingleR labels
# -----------------------------------------------------------------------------
# PURPOSE
#   S11  : the 8 whole HGSC tissue sections (core-8), one canvas, all at the
#          SAME physical scale, labelled by de-identified tissue code (WT01-08).
#   S12a : the final HGSC cores on TMA A (CHUM), each patient's duplicate cores side by side.
#   S12b : the 15 fallopian-tube epithelium (FTE) cores (on the TMA A slide).
#   S13  : the final HGSC cores on TMA B (COEUR blocks G + H).
#   + a shared 18-cell-type legend (saved as SuppFig11_/12_/13_celltype_legend).
#   Internal labels carry earlier figure numbers: console headers and the `panel` column of
#   SF10_SF11_SF12_counts.csv say SF10 (= S11), SF11a / SF11b (= S12a / S12b) and SF12 (= S13).
#
#   Colour = cell_label_singler (18 SingleR types; secretory split into
#   SecA / Intermediate / SecB on polarization_UCell with the frozen 06f
#   thresholds), house palette `ref_palette` from spatial/00_setup/00_setup.R
#   (the palette of the Figure 4 and Figure 5 cell-type maps).
#
#   EVERY cell is drawn (no subsampling); the point layer is rasterized
#   (ggrastr) so the SVG stays small with vector text. Draw order is shuffled
#   (seed from config) so no cell type is systematically hidden.
#
# COHORT
#   Whole tissues : core-8 (all > 80k cells).
#   TMA cores     : the final TMA core set (spatial/02_qc/06b_final_tma_cores.py, read through
#                   spatial/00_setup/final_tma_cores.R; author decision 2026-10-19), the cores behind
#                   every TMA panel: tumour cores that pass per-array core QC on the SingleR labels, from
#                   patients with primary HGSC, treatment-naive at sampling, with >= 500 cells on those
#                   cores -> 378 patients / 503 cores (TMA A 95 / 177, TMA B 283 / 326). Cores outside
#                   the set (core QC, neoadjuvant chemotherapy, non-HGSC histology) are not drawn; they
#                   are listed with the reason in Table S8.
#   FTE cores     : the 15 FTE cores on TMA A (sample_type == "fallopian").
#
# CORE ORDER ("grid order")
#   TMA A: core k and k+101 are the same patient (duplicate core sets 1-101 /
#          102-202). One unit per patient = two adjacent tiles (set-1 core left,
#          set-2 core right, blank if absent), units in CHUM map order of set 1,
#          9 pairs per row, code label centred over the pair with a thin bracket.
#   TMA B: physical row-major order per block (G, then H): rows are chains of
#          mutual nearest right-hand neighbours (see grid_order()), merged top ->
#          bottom, cores ordered left -> right within a row.
#          A 2-core patient (both cores always on the same block) is drawn as an
#          adjacent pair at the position of its first core in that order.
#   FTE  : by de-identified code (FTE01..FTE15).
#   Tiles are re-centred on each core's median x/y; all tiles share one scale.
#
# DE-IDENTIFICATION
#   Labels use the study codes from the PRIVATE key
#   ($OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv). The key
#   and source IDs are read here but never written to any output (file names,
#   titles, text). Tile labels show only the code.
#
# INPUTS
#   output_root/singler_recompute/<sample>_singler_cells.csv.gz
#   $OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv
#   $OUTPUT_ROOT/final_panels/supplementary/tables/SuppTable8_xenium_clinical_cohort_Whole_tissues.csv
#       (public, de-identified; treatment status per WT code)
# OUTPUTS  (output_root/figures/supplementary/SF10_SF11_singler/)
#   SuppFig11_xenium_whole_tissue_celltype.{png,svg}
#   SuppFig12a_xenium_TMA_A_paired_cores_celltype.{png,svg}
#   SuppFig12b_xenium_FTE_cores_celltype.{png,svg}
#   SuppFig13_xenium_TMA_B_cores_celltype.{png,svg}
#   SuppFig1{1,2,3}_celltype_legend.{png,svg}   (identical images)
#   SF10_SF11_SF12_counts.csv   (de-identified n's per panel)
# MANUSCRIPT PANEL(S): Figures S11, S12a, S12b, S13.
#
# RUN (repo root):
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OUTPUT_ROOT=<outputs>
#   Rscript figures/supplementary/17_figS11-13_xenium_celltype_maps.R [all|wt|tmaA|fte|tmaB|legend]
# RUNTIME TIER: moderate (~5 min; 4.3M points rendered).
# =============================================================================

suppressMessages(source("spatial/00_setup/00_setup.R"))   # CFG, cfg_path, ref_palette
source("spatial/00_setup/final_tma_cores.R")              # final_core_ids()
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrastr) })

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args)) args[1] else "all"

SG      <- cfg_path("output_root", "singler_recompute")
WT_TAB  <- file.path(Sys.getenv("OUTPUT_ROOT", "output"), "final_panels/supplementary/tables/SuppTable8_xenium_clinical_cohort_Whole_tissues.csv")
OUT     <- cfg_path("output_root", "figures", "supplementary", "SF10_SF11_singler")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DPI     <- 450
MIN_CELLS <- 500
LEVELS  <- names(ref_palette)
stopifnot(length(LEVELS) == 18)

source(file.path("config", "study_codes.R"))   # study codes via the private key (never written out)
counts <- list()

# ---------------------------------------------------------------- helpers
point_layer <- function(size) {
  geom_point_rast(aes(x, y, colour = lab), shape = 16, size = size, stroke = 0,
                  raster.dpi = DPI, dev = "ragg")
}
colour_scale <- scale_colour_manual(values = ref_palette, limits = LEVELS, drop = FALSE,
                                    guide = "none")

save_both <- function(p, stem, w, h) {
  png <- file.path(OUT, paste0(stem, ".png")); svg <- file.path(OUT, paste0(stem, ".svg"))
  ggsave(png, p, width = w, height = h, dpi = DPI, bg = "white", device = ragg::agg_png)
  ggsave(svg, p, width = w, height = h, bg = "white", device = svglite::svglite)
  cat(sprintf("  saved %s  (%.1f x %.1f in; png %.1f MB, svg %.1f MB)\n", stem, w, h,
              file.size(png) / 1e6, file.size(svg) / 1e6))
}

read_sg <- function(nm, cols) {
  d <- fread(file.path(SG, paste0(nm, "_singler_cells.csv.gz")), select = cols,
             colClasses = list(character = intersect(cols, c("patient_id", "core_id"))))
  d[, lab := factor(cell_label_singler, levels = LEVELS)]
  stopifnot(!anyNA(d$lab))                         # every cell has one of the 18 labels
  d
}

# scale bar in data units (µm); label below the bar
scalebar <- function(x1, y0, um, lwd, txt_pt, lab_gap) {
  lab <- if (um >= 1000) paste0(um / 1000, " mm") else paste0(um, " µm")
  list(annotate("segment", x = x1 - um, xend = x1, y = y0, yend = y0,
                linewidth = lwd, lineend = "butt", colour = "black"),
       annotate("text", x = x1 - um / 2, y = y0 - lab_gap, label = lab, vjust = 1,
                size = txt_pt / .pt, colour = "black"))
}

# ============================================================== Figure S11: whole tissues (internal SF10)
if (MODE %in% c("all", "wt")) {
  cat("\n=== SF10: whole tissues ===\n")
  wkey <- data.table(source_id = CFG$cohort$whole_tissue)[, study_code := study_code_of("Whole tissue", source_id)][order(study_code)]
  trt  <- fread(WT_TAB)[, .(code = `Tissue code`, trt = `Treatment status at sampling`)]
  trt[, trt_s := c("Treatment-naive" = "treatment-naive", "Post-chemotherapy" = "post-chemotherapy")[trt]]
  GAP_X <- 700; GAP_Y <- 1500; LAB_H <- 700               # µm
  tis <- list(); boxes <- list()
  for (i in seq_len(nrow(wkey))) {
    d <- read_sg(wkey$source_id[i], c("x", "y", "cell_label_singler"))
    d[, `:=`(x = x - min(x), y = y - min(y))]
    boxes[[i]] <- data.table(code = wkey$study_code[i], w = max(d$x), h = max(d$y), n = nrow(d))
    d[, code := wkey$study_code[i]]
    tis[[i]] <- d
  }
  boxes <- rbindlist(boxes)
  boxes[, row := ifelse(seq_len(.N) <= 4, 1L, 2L)]
  row_h <- boxes[, .(rh = max(h)), by = row]
  # row 1 on top (y up); each row: cells bottom-aligned on a common baseline
  base_y <- c(`2` = 0, `1` = row_h[row == 2, rh] + LAB_H + GAP_Y)
  boxes[, x0 := c(0, cumsum(w + GAP_X))[seq_len(.N)], by = row]
  boxes[, y0 := base_y[as.character(row)] + (row_h$rh[match(row, row_h$row)] - h) / 2]
  canvas <- rbindlist(lapply(seq_along(tis), function(i) {
    b <- boxes[i]; tis[[i]][, .(x = x + b$x0, y = y + b$y0, lab)]
  }))
  set.seed(CFG$seed); canvas <- canvas[sample(.N)]
  W <- max(boxes[, x0 + w]); H <- max(boxes[, y0 + h]) + LAB_H
  ymin <- -1100
  w_in <- 10; h_in <- w_in * (H - ymin) / W
  lbl <- boxes[, .(x = x0 + w / 2, y = base_y[as.character(row)] + row_h$rh[match(row, row_h$row)] + LAB_H * 0.35,
                   label = paste0(code, "  (", trt$trt_s[match(code, trt$code)], ")"))]
  p <- ggplot(canvas) + point_layer(0.16) + colour_scale +
    annotate("text", x = lbl$x, y = lbl$y, label = lbl$label, size = 7 / .pt, vjust = 0, fontface = "bold") +
    coord_fixed(xlim = c(-200, W + 200), ylim = c(ymin, H), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(4, 4, 4, 4))
  for (i in seq_len(nrow(boxes))) {
    b <- boxes[i]
    p <- p + scalebar(b$x0 + b$w, b$y0 - 250, 1000, 0.5, 5.5, 120)
  }
  save_both(p, "SuppFig11_xenium_whole_tissue_celltype", w_in, h_in)
  counts$wt <- boxes[, .(panel = "SF10", unit = code, cells = n)]
  print(boxes[, .(code, n, w_mm = round(w / 1000, 1), h_mm = round(h / 1000, 1))])
  rm(tis, canvas); invisible(gc(verbose = FALSE))
}

# ============================================================== TMA / FTE galleries
# rows of cores in physical order for one block. The COEUR arrays are hand-built
# and locally warped (neighbour angles spread over +/-10 deg), so a single global
# de-rotation merges rows. Instead: link every core to its nearest right-hand
# neighbour (400-1100 µm, |dy| < 0.45 dx; mutual-nearest so each core has <= 1
# link each way) -> row chains; chains are then merged top -> bottom into rows when
# their mean y is within 350 µm and their x ranges do not overlap. Order = row, x.
grid_order <- function(cen) {
  x <- cen$cx; y <- cen$cy; n <- length(x)
  right <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    dx <- x - x[i]; dy <- y - y[i]; dd <- sqrt(dx^2 + dy^2)
    ok <- which(dd > 400 & dd < 1100 & dx > 0 & abs(dy) < 0.45 * dx)
    if (length(ok)) right[i] <- ok[which.min(dd[ok])]
  }
  # keep a link only if i is also the nearest left-hand candidate of right[i]
  for (j in unique(na.omit(right))) {
    src <- which(right == j)
    if (length(src) > 1) { keep <- src[which.min(abs(x[src] - x[j]) + abs(y[src] - y[j]))]
                           right[setdiff(src, keep)] <- NA }
  }
  chain <- rep(NA_integer_, n); k <- 0L
  for (st in setdiff(seq_len(n), right)) {           # chain starts = cores nobody links to
    k <- k + 1L; i <- st
    while (!is.na(i) && is.na(chain[i])) { chain[i] <- k; i <- right[i] }
  }
  stopifnot(!anyNA(chain))
  ch <- data.table(chain = chain, x, y)[, .(my = mean(y), x0 = min(x), x1 = max(x)), by = chain][order(-my)]
  rows <- list()
  for (r in seq_len(nrow(ch))) {
    placed <- FALSE
    for (g in seq_along(rows)) {
      G <- rows[[g]]
      if (abs(mean(G$my) - ch$my[r]) < 350 && all(ch$x1[r] < G$x0 - 300 | ch$x0[r] > G$x1 + 300)) {
        rows[[g]] <- rbind(G, ch[r]); placed <- TRUE; break }
    }
    if (!placed) rows[[length(rows) + 1]] <- ch[r]
  }
  rows <- rows[order(-sapply(rows, function(G) mean(G$my)))]
  row_of_chain <- rep(seq_along(rows), sapply(rows, nrow)); names(row_of_chain) <- unlist(lapply(rows, `[[`, "chain"))
  r <- row_of_chain[as.character(chain)]
  cat(sprintf("    %d rows (cores/row: %s)\n", length(rows), paste(tabulate(r), collapse = " ")))
  order(r, x)
}

# Flow-layout gallery. `units` = data.table(unit, slot, core_key, label) in display
# order: a unit is one patient/core occupying 1 or 2 tiles (slot 0/1); core_key NA =
# blank tile. Units flow left -> right and wrap so a unit never breaks across rows.
# The label is centred over the unit; 2-tile units get a thin bracket under the label.
# `sections` = list(list(title=, units=)).
gallery <- function(cells, sections, row_w, pitch, unit_gap, lab_pt, title_pt, w_in, stem,
                    point_size, bar_um = 500) {
  cells <- copy(cells)
  cen <- cells[, .(mx = median(x), my = median(y)), by = core_key]
  cells[cen, on = "core_key", `:=`(x = x - i.mx, y = y - i.my)]
  half <- cells[, .(e = max(abs(c(x, y)))), by = core_key]
  cat(sprintf("    tile pitch %.0f µm; max core half-extent %.0f µm (%d cores exceed pitch/2)\n",
              pitch, max(half$e), sum(half$e > pitch / 2)))
  LAB <- pitch * 0.24; TITLE <- if (length(sections) > 1) pitch * 0.6 else 0
  tiles <- list(); ulab <- list(); titles <- list(); y_cursor <- 0; Wmax <- 0
  for (s in sections) {
    if (TITLE > 0) { titles[[length(titles) + 1]] <- data.table(y = y_cursor - TITLE * 0.35, label = s$title)
                     y_cursor <- y_cursor - TITLE }
    x <- 0; row <- 0L
    for (u in unique(s$units$unit)) {
      uu <- s$units[unit == u]; w <- nrow(uu) * pitch
      if (x > 0 && x + w > row_w + 1e-6) { row <- row + 1L; x <- 0 }
      ty <- y_cursor - LAB - row * (pitch + LAB) - pitch / 2
      tiles[[length(tiles) + 1]] <- uu[, .(core_key, tx = x + slot * pitch + pitch / 2, ty = ty)]
      ulab[[length(ulab) + 1]] <- data.table(x0 = x, x1 = x + w, y = ty + pitch / 2, label = uu$label[1],
                                             pair = nrow(uu) > 1)
      Wmax <- max(Wmax, x + w); x <- x + w + unit_gap
    }
    y_cursor <- y_cursor - (row + 1) * (pitch + LAB) - pitch * 0.15
  }
  tiles <- rbindlist(tiles); ulab <- rbindlist(ulab); titles <- rbindlist(titles)
  tiles <- tiles[!is.na(core_key)]
  stopifnot(!anyDuplicated(tiles$core_key), setequal(unique(cells$core_key), tiles$core_key))
  cells[tiles, on = "core_key", `:=`(x = x + i.tx, y = y + i.ty)]
  set.seed(CFG$seed); cells <- cells[sample(.N)]
  W <- Wmax; ybot <- y_cursor - pitch * 0.45
  p <- ggplot(cells) + point_layer(point_size) + colour_scale +
    annotate("text", x = (ulab$x0 + ulab$x1) / 2, y = ulab$y + LAB * 0.22, label = ulab$label,
             size = lab_pt / .pt, vjust = 0, colour = "grey15") +
    scalebar(W - pitch * 0.1, ybot + pitch * 0.32, bar_um, 0.45, 5.5, pitch * 0.06) +
    coord_fixed(xlim = c(0, W), ylim = c(ybot, 0), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(4, 4, 4, 4))
  pr <- ulab[pair == TRUE]
  if (nrow(pr)) p <- p + annotate("segment", x = pr$x0 + pitch * 0.12, xend = pr$x1 - pitch * 0.12,
                                  y = pr$y + LAB * 0.12, yend = pr$y + LAB * 0.12,
                                  linewidth = 0.2, colour = "grey55")
  if (nrow(titles)) p <- p + annotate("text", x = 0, y = titles$y, label = titles$label, hjust = 0,
                                      size = title_pt / .pt, fontface = "bold")
  h_in <- w_in * (0 - ybot) / W
  save_both(p, stem, w_in, h_in)
}
singletons <- function(core_key, label) data.table(unit = seq_along(core_key), slot = 0L, core_key, label)

TMA_MODES <- c("all", "tmaA", "tmaB")
if (MODE %in% c(TMA_MODES, "fte")) {
  cols <- c("x", "y", "cell_label_singler", "patient_id", "core_id", "sample_type")
  tA <- read_sg("tma", cols)[, blk := "A"]
  tG <- read_sg("COEUR_TMA_G", cols)[, blk := "G"]
  tH <- read_sg("COEUR_TMA_H", cols)[, blk := "H"]
  all_tma <- rbindlist(list(tA, tG, tH))
  all_tma[, TMA := ifelse(blk == "A", "A", "B")]
  rm(tA, tG, tH); invisible(gc(verbose = FALSE))
}

if (MODE %in% TMA_MODES) {
  cat("\n=== HGSC TMA cohort ===\n")
  tum <- all_tma[sample_type == "tumour" & !is.na(patient_id) & patient_id != "" &
                   !is.na(core_id) & core_id != "Off core"]
  off <- all_tma[(is.na(sample_type) | sample_type != "fallopian")][, .N, by = blk][
    tum[, .(on = .N), by = blk], on = "blk"][, paste0(blk, "=", N - on)]
  cat("  TMA cells not on an assigned tumour core (not drawn):", paste(off, collapse = ", "), "\n")
  # final TMA core set only (core QC pass, primary HGSC, treatment-naive, patient >= 500 cells on those cores)
  fin <- rbindlist(lapply(c(A = "tma", G = "COEUR_TMA_G", H = "COEUR_TMA_H"),
                          function(a) data.table(core_id = final_core_ids(a))), idcol = "blk")
  tum <- tum[fin, on = .(blk, core_id), nomatch = 0]
  stopifnot(uniqueN(tum[, .(blk, core_id)]) == nrow(fin))
  keep <- unique(tum[, .(TMA, patient_id)])
  tum[, core_key := paste(blk, core_id, sep = ":")]
  # de-identified code per patient
  tum[, code := study_code_of(paste("TMA", TMA[1]), patient_id), by = TMA]   # codes pass through
  stopifnot(!anyNA(tum$code))
  cores <- tum[, .(n = .N, code = code[1], cx = median(x), cy = median(y)), by = .(TMA, blk, core_key, core_id)]
  cat(sprintf("  cohort: %d patients (A %d, B %d), %d cores (A %d, B %d), %s cells; %d cores < %d cells\n",
              nrow(keep), keep[TMA == "A", .N], keep[TMA == "B", .N], nrow(cores),
              cores[TMA == "A", .N], cores[TMA == "B", .N], format(nrow(tum), big.mark = ","),
              cores[n < MIN_CELLS, .N], MIN_CELLS))
  cores[, ncore := .N, by = code]
}

if (MODE %in% c("all", "tmaA")) {
  # --- Figure S12a (internal SF11a) TMA A: duplicate cores side by side. Core k and k+101 are the same
  #     patient; units in CHUM map order of core set 1 (k), set-1 core left, set-2 core
  #     right; a single-core patient keeps its side with a blank tile on the other.
  cat("\n=== SF11a: TMA A paired cores ===\n")
  cA <- cores[TMA == "A"][, num := as.integer(core_id)]
  stopifnot(!anyNA(cA$num), all(cA$num %in% 1:202), all(cA$ncore <= 2))
  cA[, `:=`(pos = ifelse(num <= 101, num, num - 101L), slot = ifelse(num <= 101, 0L, 1L))]
  stopifnot(cA[, uniqueN(code), by = pos]$V1 == 1)           # both cores of a position = one patient
  uA <- CJ(pos = sort(unique(cA$pos)), slot = 0:1)[cA, on = .(pos, slot), `:=`(core_key = i.core_key)]
  uA[cA, on = "pos", label := i.code]
  uA <- uA[, .(unit = pos, slot, core_key, label)]
  cat(sprintf("  %d patients: %d paired, %d single-core\n", uniqueN(uA$unit),
              uA[, sum(!is.na(core_key)), by = unit][V1 == 2, .N], uA[, sum(!is.na(core_key)), by = unit][V1 == 1, .N]))
  counts$tmaA <- cA[order(code, num)][, .(panel = "SF11a", unit = paste0(code, "_core", seq_len(.N)), cells = n), by = code][, code := NULL]
  gallery(tum[TMA == "A", .(core_key, x, y, lab)], list(list(title = "", units = uA)),
          row_w = 9 * 2 * 900 + 8 * 300, pitch = 900, unit_gap = 300, lab_pt = 4.5, title_pt = 6.5, w_in = 9,
          stem = "SuppFig12a_xenium_TMA_A_paired_cores_celltype", point_size = 0.22)
}

if (MODE %in% c("all", "tmaB")) {
  # --- Figure S13 (internal SF12) TMA B: physical row-major order per block; a 2-core patient's cores are
  #     placed side by side at the position of its first core in that order (no patient
  #     spans blocks G and H).
  cat("\n=== SF12: TMA B cores ===\n")
  sec_B <- list()
  for (b in c("G", "H")) {
    cat("  block", b, "\n")
    cb <- cores[blk == b]
    cb <- cb[grid_order(cb)]
    cb[, ord := seq_len(.N)]
    cb[, `:=`(unit = min(ord), slot = rank(ord) - 1L), by = code]
    ub <- cb[order(unit, slot), .(unit, slot, core_key, label = code)]
    sec_B[[b]] <- list(title = paste0("Block ", b), units = ub)
  }
  stopifnot(cores[TMA == "B", uniqueN(blk), by = code]$V1 == 1)
  npair <- cores[TMA == "B" & ncore == 2, uniqueN(code)]
  cat(sprintf("  %d two-core patients placed side by side; %d single-core patients\n", npair,
              cores[TMA == "B" & ncore == 1, .N]))
  counts$tmaB <- cores[TMA == "B"][order(code, core_key)][, .(panel = "SF12", unit = paste0(code, "_core", seq_len(.N)), cells = n), by = code][, code := NULL]
  gallery(tum[TMA == "B", .(core_key, x, y, lab)], sec_B, row_w = 22 * 900 + 10 * 200, pitch = 900,
          unit_gap = 200, lab_pt = 4, title_pt = 6.5, w_in = 9,
          stem = "SuppFig13_xenium_TMA_B_cores_celltype", point_size = 0.19)
}

if (MODE %in% c("all", "fte")) {
  cat("\n=== SF11b: FTE cores ===\n")
  fte <- all_tma[sample_type == "fallopian"]
  fte[, code := study_code_of("FTE core", core_id)]   # FTE cores keyed by core id; codes pass through
  stopifnot(!anyNA(fte$code), uniqueN(fte$code) == 15)
  fte[, core_key := code]
  fc <- fte[, .N, by = code][order(code)]
  cat(sprintf("  %d FTE cores, %s cells (min %d, max %d per core)\n", nrow(fc),
              format(sum(fc$N), big.mark = ","), min(fc$N), max(fc$N)))
  counts$fte <- fc[, .(panel = "SF11b", unit = code, cells = N)]
  gallery(fte[, .(core_key, x, y, lab)], list(list(title = "", units = singletons(fc$code, fc$code))),
          row_w = 5 * 900 + 4 * 100, pitch = 900, unit_gap = 100, lab_pt = 7, title_pt = 7, w_in = 6.5,
          stem = "SuppFig12b_xenium_FTE_cores_celltype", point_size = 0.42, bar_um = 200)
}

# ============================================================== legend
if (MODE %in% c("all", "legend")) {
  cat("\n=== legend ===\n")
  lg <- data.table(lab = factor(LEVELS, levels = LEVELS), i = seq_along(LEVELS))
  lg[, `:=`(col = (i - 1) %/% 9, row = (i - 1) %% 9)]
  p <- ggplot(lg) +
    geom_point(aes(x = col * 2.3, y = -row, colour = lab), size = 2.6, shape = 16) +
    geom_text(aes(x = col * 2.3 + 0.18, y = -row, label = lab), hjust = 0, size = 7 / .pt) +
    annotate("text", x = -0.1, y = 1.1, label = "Cell type (SingleR)", hjust = 0, size = 7.5 / .pt,
             fontface = "bold") +
    scale_colour_manual(values = ref_palette, guide = "none") +
    coord_cartesian(xlim = c(-0.2, 4.9), ylim = c(-8.6, 1.6), expand = FALSE) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
  for (stem in c("SuppFig11_celltype_legend", "SuppFig12_celltype_legend", "SuppFig13_celltype_legend"))
    save_both(p, stem, 3.4, 1.75)
}

if (length(counts)) {
  cf <- file.path(OUT, "SF10_SF11_SF12_counts.csv")
  new <- rbindlist(counts)
  if (file.exists(cf) && MODE != "all") new <- rbind(fread(cf)[!panel %in% unique(new$panel)], new)
  fwrite(new, cf)
  print(new[, .(units = .N, cells = sum(cells), min = min(cells), median = median(cells)), by = panel])
}
cat("\nDone.\n")
