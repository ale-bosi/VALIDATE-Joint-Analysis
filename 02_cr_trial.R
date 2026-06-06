# =============================================================================
# VALIDATE Joint Analysis — 02_cr_trial.R
# =============================================================================
#
# Author: Alessandro Bosi
#
# Project:
#   Joint analysis of the VALIDATE-SWEDEHEART trial combining trial and
#   observational data to estimate causal treatment effects
#
# Purpose:
#   Loads and cleans the raw VALIDATE trial dataset. Selects and renames
#   baseline characteristics, creates treatment received and enrollment
#   status variables, recodes clinical variables, and excludes individuals
#   with missing treatment or outcome data.
#
# Input:
#   cr_validate_trial_raw   — raw VALIDATE trial data (per-patient SAS file)
#   cr_validate_pid         — PID-lopnr linkage file
#
# Output:
#   cr_trial                — cleaned trial dataset with all eligible
#                             individuals (randomised and non-randomised)
#
# Note:
#   Raw data paths point to internal KI servers and are not publicly available.
#   The analysis uses Swedish national health register data which are subject
#   to data sharing restrictions.
#
# Called by: 00_master.R
# =============================================================================

cr_trial <- cr_validate_trial_raw %>%
  
  # Step 1: Merge with lopnr linkage file
  left_join(cr_validate_pid, by = "PID") %>%
  rename_all(tolower) %>%
  
  # Step 2: Exclude individuals with withdrawn consent (missing lopnr)
  filter(!is.na(lopnr)) %>%
  
  # Step 3: Select and rename baseline characteristics
  select(
    # Identifiers
    lopnr,
    pid,
    
    # Trial status
    rnd3,                                             # Randomization assignment (0=Heparin, 1=Bivalirudin)
    bivund,                                           # Bivalirudin received during PCI
    hepund,                                           # Heparin received during PCI
    randstrat,                                        # Randomization strata (1=STEMI, 2=NSTEMI)
    interdat,                                         # Intervention date
    combo,                                            # Eligibility combination variable
    
    # Demographics
    gender             = d_gender,                   # Gender
    age,                                              # Age
    bmi,                                              # BMI
    weight_low         = weight_lo,                  # Weight <60kg
    
    # Risk factors
    smoking_status,                                   # Smoking status
    diabetes,                                         # Diabetes
    hypertension       = hyperton,                   # Hypertension
    hyperlipidemia     = hyperlip_treat_reg_xx,      # Hyperlipidemia
    
    # Previous cardiovascular disease
    previous_mi        = tidinf,                     # Previous MI
    previous_pci       = tidpci,                     # Previous PCI
    previous_cabg      = tidcabg,                    # Previous CABG
    previous_stroke,                                  # Previous stroke
    
    # Clinical presentation
    killipklass,                                      # Killip class
    killipklass2,                                     # Killip class (binary)
    cpr_before_hospital,                              # CPR before hospital
    heart_rate,                                       # Heart rate
    systolic_blood_pressure,                          # Systolic BP
    diastolic_blood_pressure,                         # Diastolic BP
    
    # Pre-PCI medications
    prepci_aspirin     = asafore,                    # Pre-PCI Aspirin
    prepci_clopidogrel = clofore,                    # Pre-PCI Clopidogrel
    prepci_ticagrelor  = ticfore,                    # Pre-PCI Ticagrelor
    prepci_prasugrel   = prafore,                    # Pre-PCI Prasugrel
    prepci_heparin     = hepfore,                    # Pre-PCI Heparin
    prepci_warfarin    = warfore,                    # Pre-PCI Warfarin
    prepci_fibrinolytic = trofore,                   # Pre-PCI Fibrinolytic/Thrombolysis
    prepci_bivalirudin = bivfore,                    # Pre-PCI Bivalirudin
    ace_inhibitors     = ace_inhibitors_reg,         # ACE inhibitors regular use
    
    # Angiography findings
    fynd               = fynd2,                      # Angiography findings
    stenosis_summary   = stenossum,                  # Stenosis summary
    stenosis_class     = stenosklass,                # Stenosis class
    
    # Additional clinical
    renal_failure      = renalfailure_x,             # Renal failure
    serum_creatinine   = skreatinin,                 # Serum creatinine (µmol/L)
    access_site        = punkt_x,                    # Vascular access site
    left_ventricular_function = left_ventricular_function_x,  # LV function
    
    # Outcome variables
    dead,                                             # All-cause death
    survtime                                          # Time to death or censoring
  ) %>%
  
  # Step 4: Create enrollment status and treatment received variables
  mutate(
    enrol = ifelse(is.na(rnd3), "Not randomised", "Randomised"),
    
    # Treatment received applies to all patients (enrolled and non-enrolled)
    treat_received = case_when(
      bivund == 1               ~ "Intervention",   # Received Bivalirudin
      bivund == 0 & hepund == 1 ~ "Control",        # Received Heparin only
      TRUE                      ~ NA_character_
    ),
    
    # Year of intervention for SES data linkage
    year_interdat = year(interdat)
  ) %>%
  
  # Step 5: Convert character and binary variables to factors
  mutate(across(c(gender, randstrat, smoking_status, diabetes, hypertension,
                  hyperlipidemia, previous_mi, previous_pci, previous_cabg,
                  previous_stroke, killipklass, killipklass2, cpr_before_hospital,
                  weight_low, prepci_aspirin, prepci_clopidogrel, prepci_ticagrelor,
                  prepci_prasugrel, prepci_heparin, prepci_warfarin,
                  prepci_fibrinolytic, prepci_bivalirudin, ace_inhibitors,
                  renal_failure, enrol, treat_received, access_site,
                  left_ventricular_function),
                as.factor)) %>%
  
  # Step 6: Set missing values to 0 for binary clinical variables
  mutate(across(c(diabetes, hypertension, hyperlipidemia, previous_mi, previous_pci,
                  previous_cabg, previous_stroke, cpr_before_hospital, weight_low,
                  prepci_aspirin, prepci_clopidogrel, prepci_ticagrelor, prepci_prasugrel,
                  prepci_heparin, prepci_warfarin, prepci_fibrinolytic,
                  prepci_bivalirudin, ace_inhibitors, renal_failure),
                ~ if_else(is.na(.), factor(0), .))) %>%
  
  # Step 7: Recode smoking status (set unknown/missing to NA)
  mutate(smoking_status = case_when(
    smoking_status == 9 ~ NA,
    smoking_status == 0 ~ factor(0),   # Never smoker
    smoking_status == 1 ~ factor(1),   # Former smoker
    smoking_status == 2 ~ factor(2),   # Current smoker
    TRUE ~ NA
  )) %>%
  
  # Step 8: Create simplified Killip class (binary: I vs II-IV)
  mutate(killip_binary = case_when(
    killipklass2 == 1  ~ factor(0),    # Killip I
    killipklass2 == 10 ~ factor(1),    # Killip II-IV
    TRUE ~ NA
  )) %>%
  
  # Step 9: Create stenosis class variable (A, B, C) from stenosis summary
  mutate(stenosclass = case_when(
    grepl("C", stenosis_summary) ~ "C",
    grepl("B", stenosis_summary) ~ "B",
    grepl("A", stenosis_summary) ~ "A",
    TRUE ~ NA_character_
  )) %>%
  
  # Step 10: Create simplified angiography finding (combine normal/inconclusive/missing)
  mutate(fynd_mod = case_when(
    fynd %in% c(0, 1) ~ 1,
    is.na(fynd)       ~ 1,
    TRUE              ~ fynd
  )) %>%
  
  # Step 11: Rename outcome variables
  mutate(
    primary_ev180 = dead,
    primary_time  = survtime
  ) %>%
  
  # Step 12: Reorder key variables
  relocate(lopnr, pid, enrol, treat_received, year_interdat, interdat, rnd3,
           randstrat, primary_ev180, primary_time) %>%
  
  # Step 13: Exclude individuals with unknown treatment received
  filter(!is.na(treat_received)) %>%
  
  # Step 14: Exclude randomised individuals without death outcome data
  filter(!(enrol == "Randomised" & is.na(dead)))