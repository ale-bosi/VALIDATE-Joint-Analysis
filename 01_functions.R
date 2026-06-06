# =============================================================================
# VALIDATE Joint Analysis — 01_functions.R
# =============================================================================
#
# Author: Alessandro Bosi
#
# Project:
#   Joint analysis of the VALIDATE-SWEDEHEART trial combining trial and
#   observational data to estimate causal treatment effects
#
# Purpose:
#   Defines all functions used across the analysis pipeline. This script
#   must be sourced before any other script.
#
#   Functions are grouped by domain:
#     1. SES data formatting     — func_formatting()
#     2. Standardization         — standard.trial, standard.obs,
#                                  standard.adjust.S, standard.transport,
#                                  standardisation, and single-model variants
#     3. IPW                     — IPW.trial.function, IPW.obs.function,
#                                  IPW.transport, IPW.adjust.S,
#                                  IPW.pooled.function
#     4. AIPW                    — AIPW.trial.function, AIPW.obs.function,
#                                  AIPW.transport, AIPW.adjust.S,
#                                  AIPW.pooled.function
#
# Note:
#   All estimator functions accept a data argument and an index argument
#   for use with boot::boot(). The confounders object must be defined in
#   the calling environment before sourcing the analysis script.
#
# Called by: 00_master.R
# =============================================================================


# =============================================================================
# 1. SES DATA FORMATTING
# =============================================================================

#' Format socioeconomic status (SES) variables from Swedish national registers
#' for a given study population. Processes civil status, birth country, LISA
#' (income, employment, family structure, region), migration, education, and
#' occupation data. All variables are derived using the 3 years prior to the
#' index date (PCI date).
#'
#' Note: SCB data is collected at the end of each calendar year. To ensure
#' temporality, all observations occurring in the same year as or after the
#' index intervention date are excluded.
#'
#' @param input_data   Named list of register dataframes (cr_civil,
#'                     cr_birthcountry, cr_lisa, cr_migration,
#'                     cr_lisa_education, cr_lisa_ssyk).
#' @param cr_inclusion Dataframe of the study population with columns:
#'                     lopnr, year_interdat.
#' @return             Named list of formatted dataframes, one per data source,
#'                     ready for left_join to the study population.
func_formatting <- function(input_data, cr_inclusion) {
  
  cr_data <- input_data
  
  # --- Civil status ---
  cr_civil_formatted <- select(cr_inclusion, lopnr, year_interdat) %>%
    left_join(cr_data$cr_civil, by = "lopnr") %>%
    mutate(civil_status = case_when(
      civil == "OG" ~ 1,  # Single
      civil == "G"  ~ 2,  # Married
      civil == "RP" ~ 2,  # Registered partner -> combine with married
      civil == "SP" ~ 3,  # Separated partner  -> combine with divorced
      civil == "S"  ~ 3,  # Divorced
      civil == "Ä"  ~ 4,  # Widow              -> combine with surviving partner
      civil == "EP" ~ 4   # Surviving partner
    )) %>%
    filter(year < year_interdat, year_interdat - year <= 3) %>%
    group_by(lopnr, year) %>%
    mutate(status_discrepancy = n_distinct(civil_status) > 1) %>%
    ungroup() %>%
    select(-status_discrepancy) %>%
    group_by(lopnr) %>%
    mutate(baseline_civil_status = civil_status[which.min(abs(year - year_interdat))]) %>%
    arrange(lopnr, year) %>%
    mutate(
      separation_3y = ifelse(lag(civil_status) == 2 & civil_status == 3, 1, NA),
      widow_3y      = ifelse(lag(civil_status) == 2 & civil_status == 4, 1, NA)
    ) %>%
    mutate(
      separation_3y = case_when(
        any(separation_3y == 1, na.rm = TRUE)  ~ 1,
        is.na(baseline_civil_status)            ~ NA_real_,
        TRUE                                    ~ 0
      ),
      widow_3y = case_when(
        any(widow_3y == 1, na.rm = TRUE)        ~ 1,
        is.na(baseline_civil_status)            ~ NA_real_,
        TRUE                                    ~ 0
      ),
      separation_widow_3y = case_when(
        separation_3y == 1 | widow_3y == 1      ~ 1,
        is.na(baseline_civil_status)            ~ NA_real_,
        TRUE                                    ~ 0
      )
    ) %>%
    ungroup()
  
  # Secondary check for separation and spousal death using diff()
  cr_civil_formatted <- cr_civil_formatted[order(cr_civil_formatted$lopnr,
                                                 cr_civil_formatted$year), ]
  cr_civil_formatted$change_detected_seperation <- with(
    cr_civil_formatted, ave(civil_status, lopnr,
                            FUN = function(x) any(diff(x) == 1 &
                                                    x[-length(x)] == 2 &
                                                    x[-1] == 3)))
  cr_civil_formatted$change_detected_widow <- with(
    cr_civil_formatted, ave(civil_status, lopnr,
                            FUN = function(x) any(diff(x) == 2 &
                                                    x[-length(x)] == 2 &
                                                    x[-1] == 4)))
  
  cr_civil_formatted <- select(cr_civil_formatted, lopnr, baseline_civil_status,
                               separation_3y, widow_3y, separation_widow_3y) %>%
    group_by(lopnr) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  # --- Birth country ---
  cr_birthcountry_formatted <- select(cr_inclusion, lopnr) %>%
    left_join(cr_data$cr_birthcountry, by = "lopnr") %>%
    mutate(birthcountry_cat = case_when(
      fodgreg5 == "Sverige"                     ~ 1,
      fodgreg5 == "Norden utom Sverige"         ~ 2,
      fodgreg5 == "EU27 utom Norden"            ~ 3,
      fodgreg5 == "Europa utom EU27 och Norden" ~ 3,
      fodgreg5 == "Nordamerika"                 ~ 4,
      fodgreg5 == "Sydamerika"                  ~ 4,
      fodgreg5 == "Afrika"                      ~ 5,
      fodgreg5 == "Asien"                       ~ 5,
      fodgreg5 == "Sovjetunionen"               ~ 5,
      fodgreg5 == "Oceanien"                    ~ 5,
      fodgreg5 == "Okänt"                       ~ NA,
      TRUE                                      ~ NA
    ))
  
  # --- LISA (family, region, income, employment) ---
  cr_lisa_formatted <- select(cr_inclusion, lopnr, year_interdat) %>%
    left_join(cr_data$cr_lisa, by = "lopnr") %>%
    filter(year < year_interdat, year_interdat - year <= 3) %>%
    group_by(lopnr, year) %>%
    mutate(
      status_discrepancy2 = n_distinct(alosdag) > 1,
      status_discrepancy3 = n_distinct(famstf) > 1,
      status_discrepancy4 = n_distinct(famtypf) > 1,
      status_discrepancy5 = n_distinct(antflyttkommun) > 1,
      status_discrepancy6 = n_distinct(kommun) > 1
    ) %>%
    ungroup() %>%
    select(-starts_with("status_discrepancy")) %>%
    mutate(
      income_disposable        = coalesce(dispink04, dispink),
      income_disposable_family = coalesce(dispinkfam04, dispinkfam),
      income_disposable_cu     = coalesce(dispinkke04, (dispinkfam / konsviktf))
    ) %>%
    select(-c(starts_with("disp"), konsviktf, konsviktf04)) %>%
    mutate(employment_status = as.numeric(coalesce(syssstat11, syssstatj,
                                                   syssstat, syssstatg))) %>%
    select(-starts_with("sysss")) %>%
    group_by(lopnr) %>%
    mutate(
      income_disposable_3ymean        = mean(income_disposable,        na.rm = TRUE),
      income_disposable_family_3ymean = mean(income_disposable_family, na.rm = TRUE),
      income_disposable_cu_3ymean     = mean(income_disposable_cu,     na.rm = TRUE),
      days_unemployed_3years_before   = sum(alosdag, na.rm = TRUE),
      unemployed_3ybefore             = ifelse(days_unemployed_3years_before >= 1, 1, 0),
      family_position = case_when(
        grepl("^1", famstf) ~ 1, grepl("^2", famstf) ~ 2,
        grepl("^3", famstf) ~ 3, grepl("^4", famstf) ~ 4,
        grepl("^0", famstf) ~ 5, TRUE ~ NA
      ),
      family_overall = case_when(
        family_position == 4 ~ 1, family_position == 3 ~ 1,
        family_position == 2 ~ 2, family_position == 1 ~ 3,
        family_position == 5 ~ 4
      ),
      family_type = case_when(
        family_position == 1 & famtypf %in% c(11, 21)          ~ 1,
        family_position == 1 & famtypf %in% c(12, 13, 22, 23)  ~ 2,
        family_position == 2 & famtypf %in% c(31, 32, 41, 42)  ~ 3,
        family_position == 3 & famtypf >= 12 & famtypf <= 42   ~ 4,
        family_position == 4 & famtypf == 50                    ~ 5,
        family_position == 5 & famtypf == 0                     ~ 6,
        TRUE                                                     ~ NA
      ),
      baseline_family_position  = family_position[which.min(abs(year - year_interdat))],
      baseline_family_overall   = family_overall[which.min(abs(year - year_interdat))],
      baseline_family_type      = family_type[which.min(abs(year - year_interdat))],
      region_number = substr(kommun, 1, 2),
      region = case_when(
        region_number == "01" ~ 1,  region_number == "03" ~ 2,
        region_number == "04" ~ 3,  region_number == "05" ~ 4,
        region_number == "06" ~ 5,  region_number == "07" ~ 6,
        region_number == "08" ~ 7,  region_number == "09" ~ 8,
        region_number == "10" ~ 9,  region_number == "12" ~ 10,
        region_number == "13" ~ 11, region_number == "14" ~ 12,
        region_number == "17" ~ 13, region_number == "18" ~ 14,
        region_number == "19" ~ 15, region_number == "20" ~ 16,
        region_number == "21" ~ 17, region_number == "22" ~ 18,
        region_number == "23" ~ 19, region_number == "24" ~ 20,
        region_number == "25" ~ 21, TRUE ~ NA
      ),
      baseline_region = region[which.min(abs(year - year_interdat))],
      swedish_region = case_when(
        substr(kommun, 1, 2) %in% c("01", "09")                              ~ 1,
        substr(kommun, 1, 2) %in% c("03","04","17","18","19","20","21")       ~ 2,
        substr(kommun, 1, 2) %in% c("05","06","08")                           ~ 3,
        substr(kommun, 1, 2) %in% c("07","10","11","12") |
          kommun %in% c("1315","1380","1381")                                 ~ 4,
        substr(kommun, 1, 2) %in% c("14","15","16") |
          kommun %in% c("1382","1383","1384") |
          substr(kommun, 1, 2) %in% c("22","23","24","25")                    ~ 5,
        TRUE                                                                   ~ NA
      ),
      baseline_swedish_region   = swedish_region[which.min(abs(year - year_interdat))],
      baseline_degurba_code     = degurba_code[which.min(abs(year - year_interdat))],
      moves_year_before         = antflyttkommun[which.min(abs(year - year_interdat))],
      moved_year_before         = ifelse(moves_year_before >= 1, 1, moves_year_before),
      baseline_employmentstatus = employment_status[which.min(abs(year - year_interdat))]
    ) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(
      income_disposable_3ymean_q5        = ntile(income_disposable_3ymean, 5),
      income_disposable_family_3ymean_q5 = ntile(income_disposable_family_3ymean, 5),
      income_disposable_cu_3ymean_q5     = ntile(income_disposable_cu_3ymean, 5)
    ) %>%
    select(lopnr, baseline_family_position, baseline_family_type, baseline_region,
           baseline_family_overall, baseline_swedish_region, baseline_degurba_code,
           moved_year_before, unemployed_3ybefore, income_disposable_3ymean,
           income_disposable_3ymean_q5, income_disposable_family_3ymean,
           income_disposable_family_3ymean_q5, income_disposable_cu_3ymean,
           income_disposable_cu_3ymean_q5, baseline_employmentstatus)
  
  # --- Migration ---
  cr_migration_formatted <- select(cr_inclusion, lopnr, year_interdat) %>%
    left_join(cr_data$cr_migration, by = "lopnr") %>%
    filter(!is.na(migration_date)) %>%
    mutate(migration_year = year(migration_date)) %>%
    filter(migration_year < year_interdat, year_interdat - migration_year <= 3) %>%
    filter(posttyp == "Inv") %>%
    group_by(lopnr) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(migration_3ybefore = 1) %>%
    select(lopnr, migration_3ybefore)
  
  # --- Education ---
  cr_lisa_education_formatted <- select(cr_inclusion, lopnr, year_interdat) %>%
    left_join(cr_data$cr_lisa_education, by = "lopnr") %>%
    filter(year < year_interdat, year_interdat - year <= 3) %>%
    mutate(
      education_length = case_when(
        sun2000niva_old_complete >= 1 & sun2000niva_old_complete <= 2 ~ 1,
        sun2000niva_old_complete >= 3 & sun2000niva_old_complete <= 4 ~ 2,
        sun2000niva_old_complete >= 5 & sun2000niva_old_complete <= 7 ~ 3,
        TRUE ~ NA
      ),
      highest_education_level = case_when(
        grepl("^0", sun2000inr) ~ 0, grepl("^1", sun2000inr) ~ 1,
        grepl("^2", sun2000inr) ~ 2, grepl("^3", sun2000inr) ~ 3,
        grepl("^4", sun2000inr) ~ 4, grepl("^5", sun2000inr) ~ 5,
        grepl("^6", sun2000inr) ~ 6, grepl("^7", sun2000inr) ~ 7,
        grepl("^8", sun2000inr) ~ 8, TRUE ~ NA
      ),
      highest_education_cat = case_when(
        highest_education_level == 0                    ~ 1,
        highest_education_level %in% c(1, 2, 3)        ~ 2,
        highest_education_level %in% c(4, 5, 6, 7)     ~ 3,
        highest_education_level == 8                    ~ 4,
        TRUE                                            ~ NA
      )
    ) %>%
    group_by(lopnr) %>%
    mutate(
      baseline_education_length = education_length[which.min(abs(year - year_interdat))],
      baseline_education_level  = highest_education_level[which.min(abs(year - year_interdat))],
      baseline_education_cat    = highest_education_cat[which.min(abs(year - year_interdat))]
    ) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(lopnr, baseline_education_length, baseline_education_level,
           baseline_education_cat)
  
  # --- Occupation ---
  cr_lisa_ssyk_formatted <- select(cr_inclusion, lopnr, year_interdat) %>%
    left_join(cr_data$cr_lisa_ssyk, by = "lopnr") %>%
    filter(year < year_interdat, year_interdat - year <= 3) %>%
    group_by(lopnr) %>%
    mutate(
      baseline_yseg  = yseg_complete[which.min(abs(year - year_interdat))],
      baseline_ssyk3 = ssyk3[which.min(abs(year - year_interdat))]
    ) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    mutate(
      baseline_yseg_cat = case_when(
        baseline_yseg %in% c("01", "02") ~ 1,
        baseline_yseg %in% c("04", "05") ~ 2,
        baseline_yseg %in% c("06", "08") ~ 3,
        baseline_yseg %in% c("03", "07") ~ 4,
        baseline_yseg == "09"            ~ 5,
        baseline_yseg == "10"            ~ 6,
        TRUE                             ~ NA
      ),
      baseline_ssyk3_value = substr(baseline_ssyk3, 1, 1)
    ) %>%
    select(lopnr, baseline_yseg, baseline_yseg_cat, baseline_ssyk3_value)
  
  list(
    cr_civil_formatted          = cr_civil_formatted,
    cr_birthcountry_formatted   = cr_birthcountry_formatted,
    cr_lisa_formatted           = cr_lisa_formatted,
    cr_migration_formatted      = cr_migration_formatted,
    cr_lisa_education_formatted = cr_lisa_education_formatted,
    cr_lisa_ssyk_formatted      = cr_lisa_ssyk_formatted
  )
}


# =============================================================================
# 2. STANDARDIZATION (G-FORMULA) FUNCTIONS
# =============================================================================

#' G-formula using trial data only — separate outcome models per arm.
#' Fits outcome models on trial participants, predicts for trial participants.
#' @param data    Dataframe passed by boot::boot().
#' @param index   Bootstrap indices.
#' @return        Named vector: Y1, Y0, effect (all in %).
standard.trial <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(trial, treat.received == "Intervention")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(trial, treat.received == "Control")[, c("death", confounders)],
                         family = binomial())
  
  Y1 <- mean(predict(outcome.treated, newdata = trial, type = "response")) * 100
  Y0 <- mean(predict(outcome.control, newdata = trial, type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula using trial data only — single outcome model with treatment as covariate.
#' @inheritParams standard.trial
standard.trial.single.model <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  
  outcome.model <- glm(death == "Yes" ~ .,
                       data   = data[, c("death", "treat.received", confounders)],
                       family = binomial())
  
  Y1 <- mean(predict(outcome.model,
                     newdata = mutate(trial, treat.received = "Intervention"),
                     type = "response")) * 100
  Y0 <- mean(predict(outcome.model,
                     newdata = mutate(trial, treat.received = "Control"),
                     type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula using observational data only — separate outcome models per arm.
#' Fits and predicts on non-randomised individuals only.
#' @inheritParams standard.trial
standard.obs <- function(data, index) {
  data <- data[index, ]
  obs  <- data |> filter(enrol == "Not randomised")
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(obs, treat.received == "Intervention")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(obs, treat.received == "Control")[, c("death", confounders)],
                         family = binomial())
  
  Y1 <- mean(predict(outcome.treated, newdata = obs, type = "response")) * 100
  Y0 <- mean(predict(outcome.control, newdata = obs, type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula using observational data only — single outcome model.
#' @inheritParams standard.trial
standard.obs.single.model <- function(data, index) {
  data <- data[index, ]
  obs  <- data |> filter(enrol == "Not randomised")
  
  outcome.model <- glm(death == "Yes" ~ .,
                       data   = obs[, c("death", "treat.received", confounders)],
                       family = binomial())
  
  Y1 <- mean(predict(outcome.model,
                     newdata = mutate(obs, treat.received = "Intervention"),
                     type = "response")) * 100
  Y0 <- mean(predict(outcome.model,
                     newdata = mutate(obs, treat.received = "Control"),
                     type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula — stratified joint analysis.
#' Fits separate outcome models for trial participants and non-participants,
#' combines predictions across the full target population.
#' @inheritParams standard.trial
standard.adjust.S <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  obs   <- data |> filter(enrol == "Not randomised")
  
  outcome.trial.treated <- glm(death == "Yes" ~ .,
                               data   = filter(trial, treat.received == "Intervention")[, c("death", confounders)],
                               family = binomial())
  outcome.trial.control <- glm(death == "Yes" ~ .,
                               data   = filter(trial, treat.received == "Control")[, c("death", confounders)],
                               family = binomial())
  outcome.obs.treated   <- glm(death == "Yes" ~ .,
                               data   = filter(obs, treat.received == "Intervention")[, c("death", confounders)],
                               family = binomial())
  outcome.obs.control   <- glm(death == "Yes" ~ .,
                               data   = filter(obs, treat.received == "Control")[, c("death", confounders)],
                               family = binomial())
  
  Y1 <- mean(c(predict(outcome.trial.treated, newdata = trial, type = "response"),
               predict(outcome.obs.treated,   newdata = obs,   type = "response"))) * 100
  Y0 <- mean(c(predict(outcome.trial.control, newdata = trial, type = "response"),
               predict(outcome.obs.control,   newdata = obs,   type = "response"))) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula — generalizability analysis.
#' Fits outcome models on trial participants only, predicts for all individuals
#' in the target population.
#' @inheritParams standard.trial
standard.transport <- function(data, index) {
  data <- data[index, ]
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Intervention",
                                         enrol == "Randomised")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Control",
                                         enrol == "Randomised")[, c("death", confounders)],
                         family = binomial())
  
  Y1 <- mean(predict(outcome.treated, newdata = data, type = "response")) * 100
  Y0 <- mean(predict(outcome.control, newdata = data, type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula — generalizability analysis, single outcome model.
#' Fits a single outcome model on trial participants including treatment as
#' covariate, predicts for all individuals in the target population.
#' @inheritParams standard.trial
standard.single.model.transport <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  
  outcome.model <- glm(death == "Yes" ~ .,
                       data   = trial[, c("death", "treat.received", confounders)],
                       family = binomial())
  
  Y1 <- mean(predict(outcome.model,
                     newdata = mutate(data, treat.received = "Intervention"),
                     type = "response")) * 100
  Y0 <- mean(predict(outcome.model,
                     newdata = mutate(data, treat.received = "Control"),
                     type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula — pooled joint analysis, separate outcome models per arm.
#' Fits outcome models on all individuals together, predicts for everyone.
#' @inheritParams standard.trial
standardisation <- function(data, index) {
  data <- data[index, ]
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Intervention")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Control")[, c("death", confounders)],
                         family = binomial())
  
  Y1 <- mean(predict(outcome.treated, newdata = data, type = "response")) * 100
  Y0 <- mean(predict(outcome.control, newdata = data, type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' G-formula — pooled joint analysis, single outcome model.
#' @inheritParams standard.trial
standard.single.model <- function(data, index) {
  data <- data[index, ]
  
  outcome.model <- glm(death == "Yes" ~ .,
                       data   = data[, c("death", "treat.received", confounders)],
                       family = binomial())
  
  Y1 <- mean(predict(outcome.model,
                     newdata = mutate(data, treat.received = "Intervention"),
                     type = "response")) * 100
  Y0 <- mean(predict(outcome.model,
                     newdata = mutate(data, treat.received = "Control"),
                     type = "response")) * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


# =============================================================================
# 3. INVERSE PROBABILITY WEIGHTING (IPW) FUNCTIONS
# =============================================================================

#' IPW using trial data only.
#' Estimates treatment weights among trial participants, fits a weighted
#' outcome model, predicts risks under each treatment.
#' @inheritParams standard.trial
IPW.trial.function <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  
  model.weight <- glm(treat.received == "Intervention" ~ .,
                      data   = trial[, c("treat.received", confounders)],
                      family = binomial())
  weights <- ifelse(trial$treat.received == "Intervention",
                    1 / model.weight$fitted.values,
                    1 / (1 - model.weight$fitted.values))
  weights <- weights / mean(weights)
  
  IPW.model <- glm(death == "Yes" ~ treat.received == "Intervention",
                   data    = trial,
                   family  = binomial(),
                   weights = weights)
  
  preds <- predict(IPW.model, newdata = trial, type = "response")
  Y1 <- mean(preds[trial$treat.received == "Intervention"]) * 100
  Y0 <- mean(preds[trial$treat.received == "Control"])      * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' IPW using observational data only.
#' @inheritParams standard.trial
IPW.obs.function <- function(data, index) {
  data <- data[index, ]
  obs  <- data |> filter(enrol == "Not randomised")
  
  model.weight <- glm(treat.received == "Intervention" ~ .,
                      data   = obs[, c("treat.received", confounders)],
                      family = binomial())
  weights <- ifelse(obs$treat.received == "Intervention",
                    1 / model.weight$fitted.values,
                    1 / (1 - model.weight$fitted.values))
  weights <- weights / mean(weights)
  
  IPW.model <- glm(death == "Yes" ~ treat.received == "Intervention",
                   data    = obs,
                   family  = binomial(),
                   weights = weights)
  
  preds <- predict(IPW.model, newdata = obs, type = "response")
  Y1 <- mean(preds[obs$treat.received == "Intervention"]) * 100
  Y0 <- mean(preds[obs$treat.received == "Control"])      * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' IPW — generalizability analysis.
#' Combines trial participation weights and treatment weights for trial
#' participants; non-participants receive weight 0.
#' @inheritParams standard.trial
IPW.transport <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  
  model.participation  <- glm(enrol == "Randomised" ~ .,
                              data   = data[, c("enrol", confounders)],
                              family = binomial())
  predict.participation <- predict(model.participation, newdata = data, type = "response")
  weight.participation  <- ifelse(data$enrol == "Randomised",
                                  1 / predict.participation, 0)
  
  model.treatment  <- glm(treat.received == "Intervention" ~ .,
                          data   = trial[, c("treat.received", confounders)],
                          family = binomial())
  predict.treatment <- predict(model.treatment, newdata = data, type = "response")
  weight.treatment  <- ifelse(data$treat.received == "Intervention",
                              1 / predict.treatment, 1 / (1 - predict.treatment))
  weight.treatment[data$enrol == "Not randomised"] <- 0
  
  weights <- weight.participation * weight.treatment
  weights <- weights / mean(weights)
  
  IPW.model <- glm(death == "Yes" ~ treat.received == "Intervention",
                   data    = data,
                   family  = binomial(),
                   weights = weights)
  
  preds <- predict(IPW.model, newdata = data, type = "response")
  Y1 <- mean(preds[data$treat.received == "Intervention"]) * 100
  Y0 <- mean(preds[data$treat.received == "Control"])      * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' IPW — stratified joint analysis.
#' Estimates treatment weights separately for trial participants and
#' non-participants, combines into a single weight vector.
#' @inheritParams standard.trial
IPW.adjust.S <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  obs   <- data |> filter(enrol == "Not randomised")
  
  model.participants <- glm(treat.received == "Intervention" ~ .,
                            data   = trial[, c("treat.received", confounders)],
                            family = binomial())
  w.trial <- ifelse(trial$treat.received == "Intervention",
                    1 / model.participants$fitted.values,
                    1 / (1 - model.participants$fitted.values))
  
  model.non.participants <- glm(treat.received == "Intervention" ~ .,
                                data   = obs[, c("treat.received", confounders)],
                                family = binomial())
  w.obs <- ifelse(obs$treat.received == "Intervention",
                  1 / model.non.participants$fitted.values,
                  1 / (1 - model.non.participants$fitted.values))
  
  weights <- rep(NA, nrow(data))
  weights[data$enrol == "Randomised"]     <- w.trial
  weights[data$enrol == "Not randomised"] <- w.obs
  weights <- weights / mean(weights)
  
  IPW.model <- glm(death == "Yes" ~ treat.received == "Intervention",
                   data    = data,
                   family  = binomial(),
                   weights = weights)
  
  preds <- predict(IPW.model, newdata = data, type = "response")
  Y1 <- mean(preds[data$treat.received == "Intervention"]) * 100
  Y0 <- mean(preds[data$treat.received == "Control"])      * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


#' IPW — pooled joint analysis.
#' Marginal treatment probability derived by combining treatment probabilities
#' in trial and non-trial strata weighted by participation probability:
#' P(A=a|X) = P(A=a|X,S=1)*P(S=1|X) + P(A=a|X,S=0)*P(S=0|X)
#' @inheritParams standard.trial
IPW.pooled.function <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  obs   <- data |> filter(enrol == "Not randomised")
  
  model.participation <- glm(enrol == "Randomised" ~ .,
                             data   = data[, c("enrol", confounders)],
                             family = binomial())
  prob.S1 <- model.participation$fitted.values
  prob.S0 <- 1 - prob.S1
  
  model.participants <- glm(treat.received == "Intervention" ~ .,
                            data   = trial[, c("treat.received", confounders)],
                            family = binomial())
  p.trt.trial <- predict(model.participants, newdata = data, type = "response")
  prob.A.trial <- ifelse(data$treat.received == "Intervention",
                         p.trt.trial, 1 - p.trt.trial)
  
  model.non.participants <- glm(treat.received == "Intervention" ~ .,
                                data   = obs[, c("treat.received", confounders)],
                                family = binomial())
  p.trt.obs <- predict(model.non.participants, newdata = data, type = "response")
  prob.A.obs <- ifelse(data$treat.received == "Intervention",
                       p.trt.obs, 1 - p.trt.obs)
  
  prob.treatment <- prob.A.trial * prob.S1 + prob.A.obs * prob.S0
  weights        <- 1 / prob.treatment
  weights        <- weights / mean(weights)
  
  IPW.model <- glm(death == "Yes" ~ treat.received == "Intervention",
                   data    = data,
                   family  = binomial(),
                   weights = weights)
  
  preds <- predict(IPW.model, newdata = data, type = "response")
  Y1 <- mean(preds[data$treat.received == "Intervention"]) * 100
  Y0 <- mean(preds[data$treat.received == "Control"])      * 100
  
  return(c(Y1 = Y1, Y0 = Y0, effect = Y1 - Y0))
}


# =============================================================================
# 4. AUGMENTED INVERSE PROBABILITY WEIGHTING (AIPW) FUNCTIONS
# =============================================================================

#' AIPW using trial data only.
#' Doubly robust estimator combining IPW weights and standardization
#' predictions. Consistent if either the outcome model or the treatment
#' probability model is correctly specified.
#' @inheritParams standard.trial
#' @return Named vector: AIPW1, AIPW0, AIPW (all in %).
AIPW.trial.function <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  
  model.weight <- glm(treat.received == "Intervention" ~ .,
                      data   = trial[, c("treat.received", confounders)],
                      family = binomial())
  weights <- ifelse(trial$treat.received == "Intervention",
                    1 / model.weight$fitted.values,
                    1 / (1 - model.weight$fitted.values))
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(trial, treat.received == "Intervention")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(trial, treat.received == "Control")[, c("death", confounders)],
                         family = binomial())
  
  predicted.Y1 <- predict(outcome.treated, newdata = trial, type = "response")
  predicted.Y0 <- predict(outcome.control, newdata = trial, type = "response")
  obs.treat    <- ifelse(trial$treat.received == "Intervention", 1, 0)
  obs.Y        <- ifelse(trial$death == "Yes", 1, 0)
  
  AIPW1 <- (sum(obs.treat * weights))^-1 *
    sum(obs.treat * weights * (obs.Y - predicted.Y1)) + mean(predicted.Y1)
  AIPW0 <- (sum((1 - obs.treat) * weights))^-1 *
    sum((1 - obs.treat) * weights * (obs.Y - predicted.Y0)) + mean(predicted.Y0)
  
  return(c(AIPW1 = AIPW1 * 100, AIPW0 = AIPW0 * 100, AIPW = (AIPW1 - AIPW0) * 100))
}


#' AIPW using observational data only.
#' @inheritParams AIPW.trial.function
AIPW.obs.function <- function(data, index) {
  data <- data[index, ]
  obs  <- data |> filter(enrol == "Not randomised")
  
  model.weight <- glm(treat.received == "Intervention" ~ .,
                      data   = obs[, c("treat.received", confounders)],
                      family = binomial())
  weights <- ifelse(obs$treat.received == "Intervention",
                    1 / model.weight$fitted.values,
                    1 / (1 - model.weight$fitted.values))
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(obs, treat.received == "Intervention")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(obs, treat.received == "Control")[, c("death", confounders)],
                         family = binomial())
  
  predicted.Y1 <- predict(outcome.treated, newdata = obs, type = "response")
  predicted.Y0 <- predict(outcome.control, newdata = obs, type = "response")
  obs.treat    <- ifelse(obs$treat.received == "Intervention", 1, 0)
  obs.Y        <- ifelse(obs$death == "Yes", 1, 0)
  
  AIPW1 <- (sum(obs.treat * weights))^-1 *
    sum(obs.treat * weights * (obs.Y - predicted.Y1)) + mean(predicted.Y1)
  AIPW0 <- (sum((1 - obs.treat) * weights))^-1 *
    sum((1 - obs.treat) * weights * (obs.Y - predicted.Y0)) + mean(predicted.Y0)
  
  return(c(AIPW1 = AIPW1 * 100, AIPW0 = AIPW0 * 100, AIPW = (AIPW1 - AIPW0) * 100))
}


#' AIPW — generalizability analysis.
#' Uses trial participation weights and treatment weights from trial participants.
#' Outcome models fitted on trial participants, predictions applied to everyone.
#' @inheritParams AIPW.trial.function
AIPW.transport <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  
  model.participation  <- glm(enrol == "Randomised" ~ .,
                              data   = data[, c("enrol", confounders)],
                              family = binomial())
  predict.participation <- predict(model.participation, newdata = data, type = "response")
  weight.participation  <- ifelse(data$enrol == "Randomised",
                                  1 / predict.participation, 0)
  
  model.treatment  <- glm(treat.received == "Intervention" ~ .,
                          data   = data[data$enrol == "Randomised",
                                        c("treat.received", confounders)],
                          family = binomial())
  predict.treatment <- predict(model.treatment, newdata = data, type = "response")
  weight.treatment  <- ifelse(data$treat.received == "Intervention",
                              1 / predict.treatment, 1 / (1 - predict.treatment))
  weight.treatment[data$enrol == "Not randomised"] <- 0
  weights <- weight.participation * weight.treatment
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Intervention",
                                         enrol == "Randomised")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Control",
                                         enrol == "Randomised")[, c("death", confounders)],
                         family = binomial())
  
  predicted.Y1 <- predict(outcome.treated, newdata = data, type = "response")
  predicted.Y0 <- predict(outcome.control, newdata = data, type = "response")
  enrol.ind    <- ifelse(data$enrol == "Randomised", 1, 0)
  obs.treat    <- ifelse(data$treat.received == "Intervention", 1, 0)
  obs.Y        <- ifelse(data$death == "Yes", 1, 0)
  
  AIPW1 <- (sum(enrol.ind * obs.treat * weights))^-1 *
    sum(enrol.ind * obs.treat * weights * (obs.Y - predicted.Y1)) + mean(predicted.Y1)
  AIPW0 <- (sum(enrol.ind * (1 - obs.treat) * weights))^-1 *
    sum(enrol.ind * (1 - obs.treat) * weights * (obs.Y - predicted.Y0)) + mean(predicted.Y0)
  
  return(c(AIPW1 = AIPW1 * 100, AIPW0 = AIPW0 * 100, AIPW = (AIPW1 - AIPW0) * 100))
}


#' AIPW — stratified joint analysis.
#' Separate treatment weight models and outcome models for trial participants
#' and non-participants. Predictions combined across the full target population.
#' @inheritParams AIPW.trial.function
AIPW.adjust.S <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  obs   <- data |> filter(enrol == "Not randomised")
  
  model.participants <- glm(treat.received == "Intervention" ~ .,
                            data   = trial[, c("treat.received", confounders)],
                            family = binomial())
  w.trial <- ifelse(trial$treat.received == "Intervention",
                    1 / model.participants$fitted.values,
                    1 / (1 - model.participants$fitted.values))
  
  model.non.participants <- glm(treat.received == "Intervention" ~ .,
                                data   = obs[, c("treat.received", confounders)],
                                family = binomial())
  w.obs <- ifelse(obs$treat.received == "Intervention",
                  1 / model.non.participants$fitted.values,
                  1 / (1 - model.non.participants$fitted.values))
  
  weights <- rep(NA, nrow(data))
  weights[data$enrol == "Randomised"]     <- w.trial
  weights[data$enrol == "Not randomised"] <- w.obs
  
  outcome.trial.treated <- glm(death == "Yes" ~ .,
                               data   = filter(trial, treat.received == "Intervention")[, c("death", confounders)],
                               family = binomial())
  outcome.trial.control <- glm(death == "Yes" ~ .,
                               data   = filter(trial, treat.received == "Control")[, c("death", confounders)],
                               family = binomial())
  outcome.obs.treated   <- glm(death == "Yes" ~ .,
                               data   = filter(obs, treat.received == "Intervention")[, c("death", confounders)],
                               family = binomial())
  outcome.obs.control   <- glm(death == "Yes" ~ .,
                               data   = filter(obs, treat.received == "Control")[, c("death", confounders)],
                               family = binomial())
  
  predicted.Y1 <- rep(NA, nrow(data))
  predicted.Y0 <- rep(NA, nrow(data))
  predicted.Y1[data$enrol == "Randomised"]     <- predict(outcome.trial.treated, newdata = trial, type = "response")
  predicted.Y1[data$enrol == "Not randomised"] <- predict(outcome.obs.treated,   newdata = obs,   type = "response")
  predicted.Y0[data$enrol == "Randomised"]     <- predict(outcome.trial.control, newdata = trial, type = "response")
  predicted.Y0[data$enrol == "Not randomised"] <- predict(outcome.obs.control,   newdata = obs,   type = "response")
  
  obs.treat <- ifelse(data$treat.received == "Intervention", 1, 0)
  obs.Y     <- ifelse(data$death == "Yes", 1, 0)
  
  AIPW1 <- (sum(obs.treat * weights))^-1 *
    sum(obs.treat * weights * (obs.Y - predicted.Y1)) + mean(predicted.Y1)
  AIPW0 <- (sum((1 - obs.treat) * weights))^-1 *
    sum((1 - obs.treat) * weights * (obs.Y - predicted.Y0)) + mean(predicted.Y0)
  
  return(c(AIPW1 = AIPW1 * 100, AIPW0 = AIPW0 * 100, AIPW = (AIPW1 - AIPW0) * 100))
}


#' AIPW — pooled joint analysis.
#' Marginal treatment probability derived by combining treatment probabilities
#' in trial and non-trial strata. Outcome models fitted on all individuals.
#' @inheritParams AIPW.trial.function
AIPW.pooled.function <- function(data, index) {
  data  <- data[index, ]
  trial <- data |> filter(enrol == "Randomised")
  obs   <- data |> filter(enrol == "Not randomised")
  
  model.participation <- glm(enrol == "Randomised" ~ .,
                             data   = data[, c("enrol", confounders)],
                             family = binomial())
  prob.S1 <- model.participation$fitted.values
  prob.S0 <- 1 - prob.S1
  
  model.participants <- glm(treat.received == "Intervention" ~ .,
                            data   = trial[, c("treat.received", confounders)],
                            family = binomial())
  p.trt.trial <- predict(model.participants, newdata = data, type = "response")
  prob.A.trial <- ifelse(data$treat.received == "Intervention",
                         p.trt.trial, 1 - p.trt.trial)
  
  model.non.participants <- glm(treat.received == "Intervention" ~ .,
                                data   = obs[, c("treat.received", confounders)],
                                family = binomial())
  p.trt.obs <- predict(model.non.participants, newdata = data, type = "response")
  prob.A.obs <- ifelse(data$treat.received == "Intervention",
                       p.trt.obs, 1 - p.trt.obs)
  
  prob.treatment <- prob.A.trial * prob.S1 + prob.A.obs * prob.S0
  weights        <- 1 / prob.treatment
  
  outcome.treated <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Intervention")[, c("death", confounders)],
                         family = binomial())
  outcome.control <- glm(death == "Yes" ~ .,
                         data   = filter(data, treat.received == "Control")[, c("death", confounders)],
                         family = binomial())
  
  predicted.Y1 <- predict(outcome.treated, newdata = data, type = "response")
  predicted.Y0 <- predict(outcome.control, newdata = data, type = "response")
  obs.treat    <- ifelse(data$treat.received == "Intervention", 1, 0)
  obs.Y        <- ifelse(data$death == "Yes", 1, 0)
  
  AIPW1 <- (sum(obs.treat * weights))^-1 *
    sum(obs.treat * weights * (obs.Y - predicted.Y1)) + mean(predicted.Y1)
  AIPW0 <- (sum((1 - obs.treat) * weights))^-1 *
    sum((1 - obs.treat) * weights * (obs.Y - predicted.Y0)) + mean(predicted.Y0)
  
  return(c(AIPW1 = AIPW1 * 100, AIPW0 = AIPW0 * 100, AIPW = (AIPW1 - AIPW0) * 100))
}