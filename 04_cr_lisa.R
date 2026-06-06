# =============================================================================
# VALIDATE Joint Analysis — 04_cr_lisa.R
# =============================================================================
#
# Author: Alessandro Bosi
#
# Project:
#   Joint analysis of the VALIDATE-SWEDEHEART trial combining trial and
#   observational data to estimate causal treatment effects
#
# Purpose:
#   Loads raw socioeconomic data from Swedish national registers (SCB),
#   applies func_formatting() to derive SES variables for the target
#   population, and merges all SES variables to target_population by lopnr.
#
# Input:
#   target_population   — dataset created in 03_cr_eligibility.R
#   Swedish register data (civil status, LISA, birth country, migration,
#   education, occupation) — internal KI servers, not publicly available
#
# Output:
#   target_population   — extended with SES variables
#
# Dependencies:
#   01_functions.R (func_formatting)
#
# Note:
#   Raw data paths point to internal KI servers and are not publicly available.
#   The analysis uses Swedish national health register data (LISA, SCB)
#   which are subject to data sharing restrictions.
#
# Called by: 00_master.R
# =============================================================================

#-------------------------------------------------------------------------------
# 1. LOAD RAW SCB DATA
#-------------------------------------------------------------------------------

cat("=== Loading SCB raw data ===\n")

# Civil status data (1990-2019)
years      <- 1990:2019
file_names <- paste0("W:/C6_Berglund/data/rawdata/SCB VALIDATE/validate_lev_civil_", years, ".sas7bdat")

data_list_civil <- lapply(seq_along(file_names), function(i) {
  df       <- read_sas(file_names[i])
  df$year  <- years[i]
  return(df)
})

cr_civil <- bind_rows(data_list_civil) %>%
  rename_all(tolower) %>%
  select(lopnr, civil, year)

rm(data_list_civil)
cat("cr_civil loaded\n")

# Demographic data
cr_demographic <- read_sas("W:/C6_Berglund/data/rawdata/SCB VALIDATE/validate_lev_demografi.sas7bdat") %>%
  rename_all(tolower) %>%
  select(lopnr, varldsdelnamn) %>%
  mutate(varldsdelnamn = ifelse(varldsdelnamn == "", NA, varldsdelnamn)) %>%
  filter(!is.na(varldsdelnamn))

cat("cr_demographic loaded\n")

# Birth country data
cr_birthcountry <- read_sas("W:/C6_Berglund/data/rawdata/SCB VALIDATE/FodLand 2024-11-29/validate_lev_fodland_eu27_2020.sas7bdat") %>%
  rename_all(tolower) %>%
  select(lopnr, fodgreg5)

cat("cr_birthcountry loaded\n")

# Region and urbanization data
region_size <- read_sas("W:/C6_Berglund/data/rawdata/SCB VALIDATE/degurba.sas7bdat") %>%
  rename_all(tolower)

# LISA data (1990-2019)
years      <- 1990:2019
file_names <- paste0("W:/C6_Berglund/data/rawdata/SCB VALIDATE/validate_lev_lisa_", years, ".sas7bdat")

data_list_lisa <- lapply(seq_along(file_names), function(i) {
  df      <- read_sas(file_names[i])
  df$year <- years[i]
  return(df)
})

cr_lisa <- bind_rows(data_list_lisa) %>%
  rename_all(tolower) %>%
  relocate(year, .after = lopnr) %>%
  mutate_if(is.character, ~ na_if(., "")) %>%
  left_join(select(region_size, kommun, degurba_code), by = "kommun")

rm(data_list_lisa, region_size)
cat("cr_lisa loaded\n")

# Migration data
cr_migration <- read_sas("W:/C6_Berglund/data/rawdata/SCB TASTE/Migrationer_240926/validate_lev_migr_240926.sas7bdat") %>%
  rename_all(tolower) %>%
  mutate_if(is.character, ~ na_if(., "")) %>%
  select(-varldsdelnamn) %>%
  mutate(migration_date = ymd(datum))

cat("cr_migration loaded\n")

# Education and occupation data (1990-2022)
years      <- 1990:2022
file_names <- paste0("W:/C6_Berglund/data/rawdata/SCB VALIDATE/validate_Utbildning/validate_lev_lisa_", years, ".sas7bdat")

data_list_lisa_utbildning <- lapply(seq_along(file_names), function(i) {
  df      <- read_sas(file_names[i])
  df$year <- years[i]
  return(df)
})

cr_lisa_utbildning <- bind_rows(data_list_lisa_utbildning) %>%
  mutate_if(is.character, ~ na_if(., "")) %>%
  mutate(across(everything(), ~ ifelse(grepl("^\\*", .), NA, .)))

# Education dataset
cr_lisa_education <- cr_lisa_utbildning %>%
  mutate(sun2000niva_old_complete = coalesce(Sun2000niva_old, Sun2000niva_Old,
                                             Sun2000Niva_old, Sun2020Niva_Old)) %>%
  select(LopNr, year, sun2000niva_old_complete, Sun2000Inr) %>%
  rename_all(tolower) %>%
  filter(!is.na(sun2000niva_old_complete) | !is.na(sun2000inr))

# Occupation dataset
cr_lisa_ssyk <- cr_lisa_utbildning %>%
  mutate(SsykAr        = as.numeric(SsykAr),
         yseg_complete = coalesce(YSEG, YSeG)) %>%
  filter(is.na(SsykAr) | SsykAr <= year) %>%
  select(LopNr, year, Ssyk3, yseg_complete) %>%
  filter(!is.na(Ssyk3) | !is.na(yseg_complete)) %>%
  rename_all(tolower)

rm(data_list_lisa_utbildning, cr_lisa_utbildning)
cat("cr_lisa_education and cr_lisa_ssyk loaded\n")

#-------------------------------------------------------------------------------
# 2. APPLY FUNC_FORMATTING AND MERGE
#-------------------------------------------------------------------------------

cat("\n=== Applying func_formatting ===\n")
cat("Input: target_population (n =", nrow(target_population), ")\n")

cr_validate_base <- list(
  cr_civil          = cr_civil,
  cr_demographic    = cr_demographic,
  cr_birthcountry   = cr_birthcountry,
  cr_lisa           = cr_lisa,
  cr_migration      = cr_migration,
  cr_lisa_education = cr_lisa_education,
  cr_lisa_ssyk      = cr_lisa_ssyk
)

ses_data <- func_formatting(input_data   = cr_validate_base,
                            cr_inclusion = target_population)

cat("SES variables formatted\n\n")

#-------------------------------------------------------------------------------
# 3. MERGE SES VARIABLES TO TARGET POPULATION
#-------------------------------------------------------------------------------

cat("=== Merging SES variables ===\n")

target_population <- target_population %>%
  left_join(ses_data$cr_civil_formatted,          by = "lopnr") %>%
  left_join(ses_data$cr_birthcountry_formatted,   by = "lopnr") %>%
  left_join(ses_data$cr_lisa_formatted,           by = "lopnr") %>%
  left_join(ses_data$cr_migration_formatted,      by = "lopnr") %>%
  left_join(ses_data$cr_lisa_education_formatted, by = "lopnr") %>%
  left_join(ses_data$cr_lisa_ssyk_formatted,      by = "lopnr")

# Clean up intermediate objects
rm(cr_civil, cr_demographic, cr_birthcountry, cr_lisa, cr_migration,
   cr_lisa_education, cr_lisa_ssyk, cr_validate_base, ses_data,
   eligible_population, trial_participants)