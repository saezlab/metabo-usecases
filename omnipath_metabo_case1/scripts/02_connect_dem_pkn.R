library(readxl)
library(dplyr)
library(tidyr)
library(knitr)

DATA_FILE <- "/Users/priscilla/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project/CCLE_case_study/data/ExtendedDataTable_DifferentialAnalysis_Shorthouse_LungMutations.xlsx"
PKN_DIR   <- "/Users/priscilla/Library/CloudStorage/OneDrive-UniversitätHeidelberg/00_project/CCLE_case_study/data"


# ── Helper ────────────────────────────────────────────────────────────────────

trunc_col <- function(x, n = 40) ifelse(nchar(x) > n, paste0(substr(x, 1, n), "…"), x)


# Expand semicolon-separated ChEBI column to a flat vector of IDs
extract_chebi <- function(dem_df) {
  dem_df |>
    filter(!is.na(chebi)) |>
    pull(chebi) |>
    strsplit(";") |>
    unlist() |>
    trimws() |>
    unique()
}

# Filter PKN edges where either endpoint is a DEM ChEBI
screen_pkn <- function(pkn, chebi_ids) {
  pkn |> filter(source %in% chebi_ids | target %in% chebi_ids)
}

# ── Load PKN ──────────────────────────────────────────────────────────────────

pkn_allosteric   <- read.csv(file.path(PKN_DIR, "pkn_allosteric.csv"))
pkn_enzyme_metab <- read.csv(file.path(PKN_DIR, "pkn_enzyme_metabolite.csv"))

# Select top 10 up-regulated (highest t) and top 10 down-regulated (lowest t)
top10_dem <- function(sheet) {
  raw <- read_excel(DATA_FILE, sheet = sheet) |>
    select(metabolite_name, logFC, t, P.Value, chebi)
  bind_rows(
    slice_max(raw, t, n = 10),
    slice_min(raw, t, n = 10)
  ) |> arrange(desc(t))
}

# ── EGFR ──────────────────────────────────────────────────────────────────────

egfr_dem <- top10_dem("EGFR_filt_limma")

egfr_chebi <- extract_chebi(egfr_dem)

egfr_allosteric <- screen_pkn(pkn_allosteric, egfr_chebi)
egfr_enzyme <- screen_pkn(pkn_enzyme_metab, egfr_chebi)

# ── KRAS ──────────────────────────────────────────────────────────────────────

kras_dem <- top10_dem("KRAS_filt_limma")

kras_chebi <- extract_chebi(kras_dem)

kras_allosteric <- screen_pkn(pkn_allosteric, kras_chebi)

kras_enzyme <- screen_pkn(pkn_enzyme_metab, kras_chebi)

