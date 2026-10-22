# Trimmed copy (2026-10-21) of 2026_final_xenium_analysis/scripts/sandbox/clean_clinical_data.R, the script that ran.
# Kept: all the code that writes data/clinical_data_clean.csv (read by explore_clinical_associations.Rmd).
# Removed: the closing console summary (column types and missingness).
# The untrimmed script is in the lab project folder.
# =============================================================================
# Clean clinical data for survival and clinical analyses
# Input:  data/clinical_data_updated.xlsx
# Output: data/clinical_data_clean.csv
# =============================================================================

library(readxl)
library(dplyr)
library(readr)

# --- Load raw data -----------------------------------------------------------
raw <- read_excel("data/clinical_data_updated.xlsx")
cat("Raw data:", nrow(raw), "rows x", ncol(raw), "cols\n")

# --- Drop rows with no clinical data ----------------------------------------
# Rows where everything except patient_id is NA
has_data <- !is.na(raw$stage) | !is.na(raw$survival_months)
df <- raw[has_data, ]
cat("After dropping empty rows:", nrow(df), "rows\n")
cat("Dropped patient_ids with no data:",
    paste(raw$patient_id[!has_data], collapse = ", "), "\n\n")

# --- Clean column names ------------------------------------------------------
df <- df %>%
  rename(treatment_to_recurrence_months = `from treatment to recurrence`)

# --- Fix treatment_status encoding -------------------------------------------
df$treatment_status <- gsub("na\u00efve", "naive", df$treatment_status)
df$treatment_status <- trimws(df$treatment_status)
cat("treatment_status values:", paste(unique(df$treatment_status), collapse = ", "), "\n")

# --- Clean ca125 -------------------------------------------------------------
# Issues: European decimals (commas), ">1000", "-99" (missing code)
df$ca125_raw <- df$ca125

df$ca125 <- gsub(",", ".", df$ca125)         # European decimal -> period
df$ca125 <- ifelse(df$ca125 == "-99", NA, df$ca125)  # -99 = missing
df$ca125_above_threshold <- grepl(">", df$ca125)     # flag >1000 entries
df$ca125 <- gsub(">", "", df$ca125)           # remove > for numeric
df$ca125 <- as.numeric(df$ca125)

cat("\nCA-125 summary:\n")
print(summary(df$ca125))
cat("Entries with >threshold flag:", sum(df$ca125_above_threshold, na.rm = TRUE), "\n")

# --- Standardize residual disease --------------------------------------------
# Original categories (French surgical notes):
#   "Non (CHIR)" = no residual / complete resection (R0)
#   "Microscopique (CHIR)" = microscopic only
#   "<= 1 cm (CHIR)" or "<1 cm (CHIR)" = optimal (<=1cm)
#   "> 1 cm, <= 2 cm (CHIR)" = suboptimal (1-2cm)
#   "> 2 cm (CHIR)" = suboptimal (>2cm)
#   "Miliaire (CHIR)" = miliary (unresectable/suboptimal)
#   "Oui (CHIR)" or "Oui (Biop)" = residual present, size unspecified

df$residual_raw <- df$residual

standardize_residual <- function(x) {
  if (is.na(x)) return(NA_character_)
  x <- trimws(x)

  # No residual disease
  if (grepl("^Non", x)) return("R0")

  # Miliary = gross suboptimal
  if (grepl("Miliaire", x) & !grepl("<=|< 1|> 1|> 2", x)) return("suboptimal_miliary")

  # Size-based categories (check most specific first)
  if (grepl("> 2 cm", x)) return("suboptimal_gt2cm")
  if (grepl("> 1 cm", x)) return("suboptimal_1to2cm")
  if (grepl("<= 1 cm|< ?1 cm|<1 cm", x)) return("optimal_le1cm")
  if (grepl("Microscopique", x)) return("optimal_microscopic")

  # "Oui" without size = residual present, unspecified
  if (grepl("Oui", x)) return("residual_unspecified")

  return(NA_character_)
}

df$residual_category <- sapply(df$residual, standardize_residual)

# Binary: optimal (R0, <=1cm, microscopic) vs suboptimal (>1cm, miliary)
df$residual_binary <- case_when(
  df$residual_category %in% c("R0", "optimal_le1cm", "optimal_microscopic") ~ "optimal",
  df$residual_category %in% c("suboptimal_gt2cm", "suboptimal_1to2cm",
                               "suboptimal_miliary") ~ "suboptimal",
  df$residual_category == "residual_unspecified" ~ NA_character_,
  TRUE ~ NA_character_
)

cat("\nResidual disease mapping:\n")
print(table(df$residual_category, useNA = "ifany"))
cat("\nBinary residual:\n")
print(table(df$residual_binary, useNA = "ifany"))

# --- Code chemo response numerically -----------------------------------------
df$chemo_response_6mo <- case_when(
  df$chemo_status_6months == "chemoresistant" ~ 0L,
  df$chemo_status_6months == "chemosensitive" ~ 1L,
  TRUE ~ NA_integer_
)

df$chemo_response_12mo <- case_when(
  df$chemo_status_1year == "chemoresistant" ~ 0L,
  df$chemo_status_1year == "chemosensitive" ~ 1L,
  TRUE ~ NA_integer_
)

cat("\nChemo response 6mo:\n")
print(table(df$chemo_response_6mo, useNA = "ifany"))
cat("Chemo response 12mo:\n")
print(table(df$chemo_response_12mo, useNA = "ifany"))

# --- Ensure numeric types ----------------------------------------------------
df$stage <- as.integer(df$stage)
df$progression_months <- as.numeric(df$progression_months)
df$progression_outcome <- as.integer(df$progression_outcome)
df$survival_months <- as.numeric(df$survival_months)
df$survival_outcome <- as.integer(df$survival_outcome)
df$age <- as.numeric(df$age)
df$treatment_to_recurrence_months <- as.numeric(df$treatment_to_recurrence_months)

# --- Stage as factor (for analyses) ------------------------------------------
df$stage_figo <- factor(df$stage,
                        levels = c(1, 2, 3, 4),
                        labels = c("I", "II", "III", "IV"))

# --- Select and order final columns -----------------------------------------
clinical_clean <- df %>%
  select(
    patient_id,
    age,
    stage,
    stage_figo,
    ca125,
    ca125_above_threshold,
    residual_category,
    residual_binary,
    treatment_status,
    # PFS
    progression_months,
    progression_outcome,
    # OS
    survival_months,
    survival_outcome,
    # Chemo response
    treatment_to_recurrence_months,
    chemo_status_6months,
    chemo_status_1year,
    chemo_response_6mo,
    chemo_response_12mo,
    # Raw values for reference
    ca125_raw = ca125_raw,
    residual_raw = residual_raw
  )

# --- Write output ------------------------------------------------------------
write_csv(clinical_clean, "data/clinical_data_clean.csv")
cat("\nWritten: data/clinical_data_clean.csv\n")
