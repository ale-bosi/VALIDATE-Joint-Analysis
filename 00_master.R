# =============================================================================
# VALIDATE Joint Analysis — 00_master.R
# =============================================================================
#
# Author: Alessandro Bosi
#
# Project:
#   Joint analysis of the VALIDATE-SWEDEHEART trial combining trial and
#   observational data to estimate causal treatment effects
#
# Purpose:
#   Master script for the full analysis pipeline. Sources all sub-scripts
#   in order, from functions and raw data loading to analysis and results.
#
# Note:
#   Raw data paths point to internal KI servers and are not publicly available.
#   The analysis uses Swedish national health register data which are subject
#   to data sharing restrictions.
# =============================================================================

# PACKAGES ---------------------------------------------------------------------
library(dplyr)
library(haven)
library(tidyr)
library(tidyverse)
library(lubridate)
library(survival)
library(boot)
library(mada)
library(table1)
library(smd)
library(metafor)

# SET PATHS --------------------------------------------------------------------
setwd("W:/C6_Berglund/abosi/VALIDATE_Joint_analysis/Script/")

# FUNCTIONS --------------------------------------------------------------------
source("01_functions.R")

# DATA LOADING -----------------------------------------------------------------
cr_validate_trial_raw <- read_sas("W:/C6_Berglund/data/rawdata/validate-trial/perpatient.sas7bdat")
cr_validate_pid       <- read_sas("W:/C6_Berglund/data/rawdata/SCB Validate/validate_pid.sas7bdat")

# DATA CLEANING ----------------------------------------------------------------
source("02_cr_trial.R")
source("03_cr_eligibility.R")
source("04_cr_lisa.R")
source("05_cr_preparation.R")

# DESCRIPTIVE ------------------------------------------------------------------
source("06_table1.R")

# ANALYSIS ---------------------------------------------------------------------
source("07_analysis.R")