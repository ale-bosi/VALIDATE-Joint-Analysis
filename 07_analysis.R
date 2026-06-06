# =============================================================================
# VALIDATE Joint Analysis — 07_analysis.R
# =============================================================================
#
# Author: Alessandro Bosi
#
# Project:
#   Joint analysis of the VALIDATE-SWEDEHEART trial combining trial and
#   observational data to estimate causal treatment effects
#
# Purpose:
#   Runs all 15 analyses (3 estimators x 5 identification strategies) using
#   1000 nonparametric bootstrap replicates. Produces a results table and
#   a combined forest plot.
#
#   Estimators:  Standardization (G-formula), IPW, AIPW
#   Strategies:  Trial only, Observational only, Generalizability,
#                Stratified joint, Pooled joint
#
#   Outcome: All-cause mortality at 180 days.
#
# Input:
#   target_population   — prepared dataset from 05_cr_preparation.R
#   covariates          — character vector from 05_cr_preparation.R
#
# Output:
#   results             — 3x5 results table (estimate and 95% CI)
#   combined_forest_plot.png
#   analysis_results.RData
#
# Dependencies:
#   01_functions.R
#
# Called by: 00_master.R
# =============================================================================

# =============================================================================
# STEP 1: Prepare analysis dataset
# =============================================================================

validate.cc <- target_population %>%
  mutate(
    enrol          = enrol,
    treat.received = treat_received,
    death          = ifelse(primary_ev180 == 1, "Yes", "No")
  ) %>%
  select(all_of(covariates), enrol, treat.received, death)

# =============================================================================
# STEP 2: Check for missing data
# =============================================================================

missing_summary <- validate.cc %>%
  summarise(across(all_of(covariates), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0)

if (nrow(missing_summary) > 0) {
  print(missing_summary)
  stop("Missing data detected in covariates — check 05_cr_preparation.R.")
}

# =============================================================================
# STEP 3: Set confounders
# =============================================================================

confounders <- covariates

# =============================================================================
# STEP 4: Run all 15 analyses (R = 1000 bootstrap replicates)
# =============================================================================

# Standardization
set.seed(2025); std.trial    <- boot(validate.cc, standard.trial,     R = 1000)
set.seed(2026); std.obs.only <- boot(validate.cc, standard.obs,       R = 1000)
set.seed(2027); std.transp   <- boot(validate.cc, standard.transport, R = 1000)
set.seed(2028); std.obs      <- boot(validate.cc, standard.adjust.S,  R = 1000)
set.seed(2029); std.joint    <- boot(validate.cc, standardisation,    R = 1000)

# IPW
set.seed(2030); ipw.trial    <- boot(validate.cc, IPW.trial.function,  R = 1000)
set.seed(2031); ipw.obs.only <- boot(validate.cc, IPW.obs.function,    R = 1000)
set.seed(2032); ipw.transp   <- boot(validate.cc, IPW.transport,       R = 1000)
set.seed(2033); ipw.obs      <- boot(validate.cc, IPW.adjust.S,        R = 1000)
set.seed(2034); ipw.joint    <- boot(validate.cc, IPW.pooled.function, R = 1000)

# AIPW
set.seed(2035); aipw.trial    <- boot(validate.cc, AIPW.trial.function,  R = 1000)
set.seed(2036); aipw.obs.only <- boot(validate.cc, AIPW.obs.function,    R = 1000)
set.seed(2037); aipw.transp   <- boot(validate.cc, AIPW.transport,       R = 1000)
set.seed(2038); aipw.obs      <- boot(validate.cc, AIPW.adjust.S,        R = 1000)
set.seed(2039); aipw.joint    <- boot(validate.cc, AIPW.pooled.function, R = 1000)

# =============================================================================
# STEP 5: Extract results
# =============================================================================

extract_res <- function(boot_obj, index = 3) {
  c(estimate = boot_obj$t0[index],
    ci_lower  = as.numeric(quantile(boot_obj$t[, index], 0.025, na.rm = TRUE)),
    ci_upper  = as.numeric(quantile(boot_obj$t[, index], 0.975, na.rm = TRUE)))
}

std_trial_res      <- extract_res(std.trial)
std_obs_res        <- extract_res(std.obs.only)
std_transp_res     <- extract_res(std.transp)
std_stratified_res <- extract_res(std.obs)
std_pooled_res     <- extract_res(std.joint)

ipw_trial_res      <- extract_res(ipw.trial)
ipw_obs_res        <- extract_res(ipw.obs.only)
ipw_transp_res     <- extract_res(ipw.transp)
ipw_stratified_res <- extract_res(ipw.obs)
ipw_pooled_res     <- extract_res(ipw.joint)

aipw_trial_res      <- extract_res(aipw.trial)
aipw_obs_res        <- extract_res(aipw.obs.only)
aipw_transp_res     <- extract_res(aipw.transp)
aipw_stratified_res <- extract_res(aipw.obs)
aipw_pooled_res     <- extract_res(aipw.joint)

# =============================================================================
# STEP 6: Create results table
# =============================================================================

fmt <- function(res) {
  paste0(format(round(res[1], 3), nsmall = 3), " (",
         format(round(res[2], 3), nsmall = 3), ", ",
         format(round(res[3], 3), nsmall = 3), ")")
}

results <- data.frame(
  `Trial only`         = c(fmt(std_trial_res),      fmt(ipw_trial_res),      fmt(aipw_trial_res)),
  `Observational only` = c(fmt(std_obs_res),        fmt(ipw_obs_res),        fmt(aipw_obs_res)),
  `Generalizability`   = c(fmt(std_transp_res),     fmt(ipw_transp_res),     fmt(aipw_transp_res)),
  `Stratified joint`   = c(fmt(std_stratified_res), fmt(ipw_stratified_res), fmt(aipw_stratified_res)),
  `Pooled joint`       = c(fmt(std_pooled_res),     fmt(ipw_pooled_res),     fmt(aipw_pooled_res)),
  row.names            = c("Standardization", "IPW", "AIPW"),
  check.names          = FALSE
)

print(results)

# =============================================================================
# STEP 7: Forest plot
# =============================================================================

yi <- c(
  std_trial_res[1],      ipw_trial_res[1],      aipw_trial_res[1],
  std_obs_res[1],        ipw_obs_res[1],        aipw_obs_res[1],
  std_transp_res[1],     ipw_transp_res[1],     aipw_transp_res[1],
  std_stratified_res[1], ipw_stratified_res[1], aipw_stratified_res[1],
  std_pooled_res[1],     ipw_pooled_res[1],     aipw_pooled_res[1]
)

ci.lb <- c(
  std_trial_res[2],      ipw_trial_res[2],      aipw_trial_res[2],
  std_obs_res[2],        ipw_obs_res[2],        aipw_obs_res[2],
  std_transp_res[2],     ipw_transp_res[2],     aipw_transp_res[2],
  std_stratified_res[2], ipw_stratified_res[2], aipw_stratified_res[2],
  std_pooled_res[2],     ipw_pooled_res[2],     aipw_pooled_res[2]
)

ci.ub <- c(
  std_trial_res[3],      ipw_trial_res[3],      aipw_trial_res[3],
  std_obs_res[3],        ipw_obs_res[3],        aipw_obs_res[3],
  std_transp_res[3],     ipw_transp_res[3],     aipw_transp_res[3],
  std_stratified_res[3], ipw_stratified_res[3], aipw_stratified_res[3],
  std_pooled_res[3],     ipw_pooled_res[3],     aipw_pooled_res[3]
)

rows <- c(23:21, 18:16, 13:11, 8:6, 3:1)
slab <- rep(c("Standardization", "IPW", "AIPW"), 5)

png("W:/C6_Berglund/abosi/VALIDATE_Joint_analysis/Output/combined_forest_plot.png",
    width = 800, height = 1000)

metafor::forest(yi, ci.lb = ci.lb, ci.ub = ci.ub,
                rows   = rows,
                header = c("", "Risk difference [95% CI]"),
                slab   = slab,
                alim   = c(-5, 5), xlim = c(-10, 8),
                xlab   = "Estimated treatment effect (%)",
                pch    = 19, psize = 1.2, cex = 2)

abline(h = max(rows) + 1, col = "white", lwd = 3)
text(-10, 24, "Trial only",               pos = 4, font = 2, cex = 2)
text(-10, 19, "Observational only",        pos = 4, font = 2, cex = 2)
text(-10, 14, "Generalizability analysis", pos = 4, font = 2, cex = 2)
text(-10, 9,  "Stratified joint analysis", pos = 4, font = 2, cex = 2)
text(-10, 4,  "Pooled joint analysis",     pos = 4, font = 2, cex = 2)

dev.off()

# =============================================================================
# STEP 8: Save workspace
# =============================================================================

save(std.trial, std.obs.only, std.transp, std.obs, std.joint,
     ipw.trial, ipw.obs.only, ipw.transp, ipw.obs, ipw.joint,
     aipw.trial, aipw.obs.only, aipw.transp, aipw.obs, aipw.joint,
     results,
     file = "W:/C6_Berglund/abosi/VALIDATE_Joint_analysis/Output/analysis_results.RData")