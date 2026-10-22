# Study codes for samples and patients (de-identification).
#
# Everything published uses study codes: WT01-WT08 (whole tissues), TMA-P### (TMA patients) and FTE##
# (fallopian-tube cores). Data held in the lab may still carry the source identifiers; the mapping lives in
# a private key file that is never stored in this repository (path: $STUDY_CODE_KEY, by default
# $OUTPUT_ROOT/private_keys/SuppTable8_code_key_PRIVATE.csv).
#
# study_code_of(arm, x) maps identifiers to codes and passes codes through, so scripts run unchanged on
# the lab's data and on the deposited, already-coded data (where the key is not needed). Use it rather
# than joining on study_code_table(arms): that is the key itself (arm, source_id, study_code, plus each
# code mapped to itself) and is empty without the key file. FTE cores are keyed by core id.
#
#   source("config/study_codes.R")
#   study_code_of("Whole tissue", CFG$cohort$whole_tissue)
#   tum[, code := study_code_of(paste("TMA", TMA[1]), patient_id), by = TMA]
suppressPackageStartupMessages(library(data.table))

STUDY_CODE_KEY <- Sys.getenv("STUDY_CODE_KEY",
  file.path(Sys.getenv("OUTPUT_ROOT", "output"), "private_keys", "SuppTable8_code_key_PRIVATE.csv"))
STUDY_CODE_RX <- c("Whole tissue" = "^WT[0-9]{2}$", "TMA A" = "^TMA-P[0-9]{3}$",
                   "TMA B" = "^TMA-P[0-9]{3}$", "FTE core" = "^FTE[0-9]{2}$")

study_code_table <- function(arms) {
  k <- if (file.exists(STUDY_CODE_KEY)) fread(STUDY_CODE_KEY, colClasses = "character")
       else data.table(arm = character(), study_code = character(), source_id = character())
  k <- k[arm %in% arms]
  k[arm == "FTE core", source_id := sub(".*core ", "", source_id)]
  ident <- unique(k[, .(arm, source_id = study_code, study_code)])
  unique(rbind(k[, .(arm, source_id, study_code)], ident))
}

study_code_of <- function(arm, x) {
  x <- sub("\\.0$", "", as.character(x))
  k <- study_code_table(arm)
  out <- k$study_code[match(x, k$source_id)]
  coded <- is.na(out) & grepl(STUDY_CODE_RX[[arm]], x)
  out[coded] <- x[coded]
  if (anyNA(out)) stop("no study code for ", arm, ": ", paste(head(x[is.na(out)], 3), collapse = ", "),
                       " (is the private key at $STUDY_CODE_KEY?)")
  out
}
