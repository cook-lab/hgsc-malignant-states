# ============================================================================
# 01_parse_metadata.R
# Parse the GSE140082 series_matrix.txt.gz file into a tidy sample-level
# metadata table.
#
# Outputs:
#   data/processed/sample_metadata.tsv  (n=380 rows)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
})

ROOT <- "/Users/dpcook/My Drive/Lab/Manuscripts/hgsc_epitypes/supplemental_analysis/ICON7"
SM_FILE <- file.path(ROOT, "data/raw/GSE140082_series_matrix.txt.gz")
OUT_FILE <- file.path(ROOT, "data/processed/sample_metadata.tsv")

dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

# Series matrix files have a "!Sample_*" header block (key TAB val1 TAB val2 ...).
# Read only the sample header section (before the !series_matrix_table_begin line).
con <- gzfile(SM_FILE, "rt")
raw <- readLines(con)
close(con)

i_begin <- grep("^!series_matrix_table_begin", raw)
if (length(i_begin) == 0) stop("series_matrix_table_begin marker not found")
header_lines <- raw[seq_len(i_begin - 1)]
sample_lines <- grep("^!Sample_", header_lines, value = TRUE)

# Split each line on tab
parse_line <- function(ln) {
  parts <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  key <- sub("^!", "", parts[1])
  vals <- gsub('^"|"$', '', parts[-1])
  list(key = key, vals = vals)
}
parsed <- lapply(sample_lines, parse_line)

# Pull out scalar (one-per-sample) fields keyed by Sample_*
# Multiple Sample_characteristics_ch1 rows: each row holds one "tag: value" pair
# per sample. We need to distribute them by the tag prefix.
n_samples <- length(parsed[[1]]$vals)

base_fields <- list()
char_rows <- list()
for (p in parsed) {
  if (p$key == "Sample_characteristics_ch1") {
    char_rows[[length(char_rows) + 1]] <- p$vals
  } else if (length(p$vals) == n_samples) {
    # Some keys repeat (e.g. Sample_supplementary_file_1/_2). Keep the first.
    if (is.null(base_fields[[p$key]])) base_fields[[p$key]] <- p$vals
  }
}

# Build a char-tag matrix: each row of char_rows holds entries like
# "treatment: bevacizumab" or "age: 53". Parse each cell into tag and value.
parse_tag <- function(x) {
  # Find the FIRST ": " split
  m <- regexpr(": ", x, fixed = TRUE)
  if (any(m < 0)) {
    # rows that don't have ":" are unusual — keep tag NA, value=raw
    out_tag <- ifelse(m < 0, NA_character_, substr(x, 1, m - 1))
    out_val <- ifelse(m < 0, x, substring(x, m + 2))
  } else {
    out_tag <- substr(x, 1, m - 1)
    out_val <- substring(x, m + 2)
  }
  list(tag = out_tag, val = out_val)
}

char_long <- data.table()
for (ri in seq_along(char_rows)) {
  parsed_row <- parse_tag(char_rows[[ri]])
  char_long <- rbind(char_long,
    data.table(sample_index = seq_len(n_samples),
               tag = parsed_row$tag,
               val = parsed_row$val))
}

# A given tag may appear only on some samples (different row positions). We
# collapse by sample_index + tag; if duplicate tag for one sample, keep first.
char_long <- char_long[!is.na(tag)]
char_wide <- dcast(char_long, sample_index ~ tag,
                   value.var = "val",
                   fun.aggregate = function(v) v[1])

# Add base fields
base_df <- data.table(sample_index = seq_len(n_samples))
for (k in names(base_fields)) {
  base_df[[k]] <- base_fields[[k]]
}

meta <- merge(base_df, char_wide, by = "sample_index", all.x = TRUE)

# Common renaming
setnames(meta,
  old = c("Sample_title", "Sample_geo_accession"),
  new = c("sample_title", "gsm"),
  skip_absent = TRUE)

# Standardize numeric columns
numeric_cols <- intersect(
  c("age", "final_pfstm", "final_pfsid", "final_ostm", "final_osid",
    "manuscript_analysis359"),
  names(meta))
for (cc in numeric_cols) {
  meta[[cc]] <- suppressWarnings(as.numeric(meta[[cc]]))
}

# Cleanup whitespace
chr_cols <- setdiff(names(meta), numeric_cols)
for (cc in chr_cols) {
  if (is.character(meta[[cc]])) {
    meta[[cc]] <- trimws(meta[[cc]])
  }
}

message(sprintf("Parsed metadata for %d samples", nrow(meta)))
message("Column summary:")
for (cc in names(meta)) {
  if (is.numeric(meta[[cc]])) {
    message(sprintf("  %-30s numeric n_nonNA=%d range=[%s, %s]",
                    cc, sum(!is.na(meta[[cc]])),
                    format(suppressWarnings(min(meta[[cc]], na.rm = TRUE))),
                    format(suppressWarnings(max(meta[[cc]], na.rm = TRUE)))))
  } else {
    u <- unique(meta[[cc]])
    show_u <- if (length(u) <= 8) paste(u, collapse = " | ") else paste(length(u), "unique")
    message(sprintf("  %-30s chr     %s", cc, show_u))
  }
}

fwrite(meta, OUT_FILE, sep = "\t")
message(sprintf("Wrote: %s", OUT_FILE))
