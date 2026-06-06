# =============================================================================
# VALIDATE Joint Analysis — 03_cr_eligibility.R
# =============================================================================
#
# Author: Alessandro Bosi
#
# Project:
#   Joint analysis of the VALIDATE-SWEDEHEART trial combining trial and
#   observational data to estimate causal treatment effects
#
# Purpose:
#   Defines the target population by identifying eligible non-participants
#   using the combo variable, which captures the reason for non-participation.
#   Eligible non-participants are those who met the trial eligibility criteria
#   but declined participation, could not be asked, or were not approached.
#   The target population is formed by combining trial participants and
#   eligible non-participants.
#
# Input:
#   cr_trial            — cleaned trial dataset from 02_cr_trial.R
#
# Output:
#   trial_participants  — randomised individuals only
#   eligible_population — eligible non-participants only
#   target_population   — combined dataset (trial + eligible non-participants)
#
# Called by: 00_master.R
# =============================================================================

# Trial participants (randomised individuals)
trial_participants <- cr_trial %>%
  filter(enrol == "Randomised")

# Eligible non-participants: individuals who met eligibility criteria but
# declined participation, could not be asked, or were not approached
eligible_population <- cr_trial %>%
  filter(enrol == "Not randomised") %>%
  filter(combo %in% c("Uppfyllde kriterierna: Avböjde deltagande",
                      "Uppfyllde kriterierna: Kan ej tillfrågas",
                      "Uppfyllde kriterierna: Har ej tillfrågats"))

# Target population: trial participants + eligible non-participants
target_population <- bind_rows(trial_participants, eligible_population)

# Clean up environment
rm(list = setdiff(ls(), c("trial_participants", "eligible_population",
                          "target_population")))