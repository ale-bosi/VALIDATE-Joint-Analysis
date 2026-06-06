# =============================================================================
# VALIDATE Joint Analysis — 05_cr_preparation.R
# =============================================================================
#
# Author: Alessandro Bosi
#
# Project:
#   Joint analysis of the VALIDATE-SWEDEHEART trial combining trial and
#   observational data to estimate causal treatment effects
#
# Purpose:
#   Prepares the target population for analysis. Steps include: recoding
#   clinical variables and pre-PCI medications, relabeling SES variables,
#   checking and collapsing sparse categories, defining the covariate list,
#   and selecting the final variables for analysis.
#
# Input:
#   target_population   — dataset from 04_cr_lisa.R
#
# Output:
#   target_population   — fully prepared dataset
#   covariates          — character vector of covariate names
#
# Called by: 00_master.R
# =============================================================================

cat("=== Preparing target_population for analysis ===\n")
cat("Initial observations:", nrow(target_population), "\n\n")

# =============================================================================
# STEP 1: Recode clinical variables
# =============================================================================

cat("Step 1: Recoding clinical variables...\n")

target_population <- target_population %>%
  mutate(
    # MI type
    stemi = case_when(
      randstrat == "1" ~ "STEMI",
      randstrat == "2" ~ "NSTEMI",
      is.na(randstrat)  ~ "Unknown",
      TRUE              ~ "Unknown"
    ),
    stemi = factor(stemi, levels = c("STEMI", "NSTEMI", "Unknown")),
    
    # Sex
    female = case_when(
      gender == "1" ~ "Male",
      gender == "2" ~ "Female",
      TRUE          ~ "Male"
    ),
    female = factor(female, levels = c("Male", "Female")),
    
    # Age >= 65
    age_65 = case_when(
      is.na(age)  ~ "Unknown",
      age >= 65   ~ "Yes",
      age < 65    ~ "No",
      TRUE        ~ "Unknown"
    ),
    age_65 = factor(age_65, levels = c("No", "Yes", "Unknown")),
    
    # BMI categories
    bmi_cat = case_when(
      is.na(bmi)              ~ "Unknown",
      bmi < 18.5              ~ "Underweight",
      bmi >= 18.5 & bmi < 25 ~ "Normal",
      bmi >= 25 & bmi < 30   ~ "Overweight",
      bmi >= 30               ~ "Obesity",
      TRUE                    ~ "Unknown"
    ),
    bmi_cat = factor(bmi_cat, levels = c("Underweight", "Normal", "Overweight",
                                         "Obesity", "Unknown")),
    
    # Smoking status
    smoking_status = case_when(
      smoking_status == "0"                     ~ "Never smoker",
      smoking_status == "1"                     ~ "Previous smoker",
      smoking_status == "2"                     ~ "Current smoker",
      smoking_status == "9" | is.na(smoking_status) ~ "Missing",
      TRUE                                      ~ "Missing"
    ),
    smoking_status = factor(smoking_status, levels = c("Never smoker", "Previous smoker",
                                                       "Current smoker", "Missing")),
    
    # Comorbidities
    diabetes = case_when(
      diabetes == "0" ~ "No", diabetes == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    diabetes = factor(diabetes, levels = c("No", "Yes", "Unknown")),
    
    hypertension = case_when(
      hypertension == "0" ~ "No", hypertension == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    hypertension = factor(hypertension, levels = c("No", "Yes", "Unknown")),
    
    hyperlipidemia = case_when(
      hyperlipidemia == "0" ~ "No", hyperlipidemia == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    hyperlipidemia = factor(hyperlipidemia, levels = c("No", "Yes", "Unknown")),
    
    previous_mi = case_when(
      previous_mi == "0" ~ "No", previous_mi == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    previous_mi = factor(previous_mi, levels = c("No", "Yes", "Unknown")),
    
    previous_pci = case_when(
      previous_pci == "0" ~ "No", previous_pci == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    previous_pci = factor(previous_pci, levels = c("No", "Yes", "Unknown")),
    
    previous_cabg = case_when(
      previous_cabg == "0" ~ "No", previous_cabg == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    previous_cabg = factor(previous_cabg, levels = c("No", "Yes", "Unknown")),
    
    previous_stroke = case_when(
      previous_stroke == "0" ~ "No", previous_stroke == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    previous_stroke = factor(previous_stroke, levels = c("No", "Yes", "Unknown")),
    
    cpr_before_hospital = case_when(
      cpr_before_hospital == "0" ~ "No", cpr_before_hospital == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    cpr_before_hospital = factor(cpr_before_hospital, levels = c("No", "Yes", "Unknown")),
    
    # Killip class (binary)
    killip_class = case_when(
      killip_binary == "0"  ~ "Killip Class I",
      killip_binary == "1"  ~ "Killip Class II, III, or IV",
      is.na(killip_binary)  ~ "Missing",
      TRUE                  ~ "Missing"
    ),
    killip_class = factor(killip_class, levels = c("Killip Class I",
                                                   "Killip Class II, III, or IV",
                                                   "Missing")),
    
    # Systolic blood pressure
    systolic_bp_cat = case_when(
      is.na(systolic_blood_pressure)                                        ~ "Missing",
      systolic_blood_pressure >= 30  & systolic_blood_pressure < 120        ~ "30 \u2264 Sys BP < 120",
      systolic_blood_pressure >= 120 & systolic_blood_pressure < 140        ~ "120 \u2264 Sys BP < 140",
      systolic_blood_pressure >= 140 & systolic_blood_pressure < 160        ~ "140 \u2264 Sys BP < 160",
      systolic_blood_pressure >= 160                                        ~ "Sys BP \u2265 160",
      TRUE                                                                  ~ "Missing"
    ),
    systolic_bp_cat = factor(systolic_bp_cat,
                             levels = c("30 \u2264 Sys BP < 120", "120 \u2264 Sys BP < 140",
                                        "140 \u2264 Sys BP < 160", "Sys BP \u2265 160", "Missing")),
    
    # Diastolic blood pressure
    diastolic_bp_cat = case_when(
      is.na(diastolic_blood_pressure)                                         ~ "Missing",
      diastolic_blood_pressure >= 25 & diastolic_blood_pressure < 60          ~ "25 \u2264 Dia BP < 60",
      diastolic_blood_pressure >= 60 & diastolic_blood_pressure < 80          ~ "60 \u2264 Dia BP < 80",
      diastolic_blood_pressure >= 80 & diastolic_blood_pressure < 100         ~ "80 \u2264 Dia BP < 100",
      diastolic_blood_pressure >= 100                                         ~ "Dia BP \u2265 100",
      TRUE                                                                    ~ "Missing"
    ),
    diastolic_bp_cat = factor(diastolic_bp_cat,
                              levels = c("25 \u2264 Dia BP < 60", "60 \u2264 Dia BP < 80",
                                         "80 \u2264 Dia BP < 100", "Dia BP \u2265 100", "Missing")),
    
    # Serum creatinine
    serum_creatinine_cat = case_when(
      is.na(serum_creatinine)                            ~ "Missing",
      serum_creatinine < 60                              ~ "<60",
      serum_creatinine >= 60  & serum_creatinine < 80   ~ "60 \u2264 SCr < 80",
      serum_creatinine >= 80  & serum_creatinine < 100  ~ "80 \u2264 SCr < 100",
      serum_creatinine >= 100                            ~ "SCr \u2265 100",
      TRUE                                               ~ "Missing"
    ),
    serum_creatinine_cat = factor(serum_creatinine_cat,
                                  levels = c("<60", "60 \u2264 SCr < 80",
                                             "80 \u2264 SCr < 100", "SCr \u2265 100", "Missing")),
    
    # Vascular access site
    femoral_access = case_when(
      access_site == "88" | access_site == "1" | is.na(access_site) ~ "Femoral",
      access_site == "2"                                             ~ "Radial",
      TRUE                                                           ~ "Femoral"
    ),
    femoral_access = factor(femoral_access, levels = c("Femoral", "Radial")),
    
    # Angiography finding
    angiography_finding = case_when(
      fynd_mod == 2                                    ~ "1 vessel",
      fynd_mod == 3                                    ~ "2 vessels",
      fynd_mod == 4                                    ~ "3 vessels",
      fynd_mod %in% c(5, 6, 7, 8)                     ~ "Left main",
      fynd_mod %in% c(0, 1, 9, 10) | is.na(fynd_mod) ~ "Missing",
      TRUE                                             ~ "Missing"
    ),
    angiography_finding = factor(angiography_finding,
                                 levels = c("1 vessel", "2 vessels", "3 vessels",
                                            "Left main", "Missing")),
    
    # Renal failure
    renal_failure = case_when(
      renal_failure == "0" ~ "No", renal_failure == "1" ~ "Yes",
      TRUE ~ "Unknown"
    ),
    renal_failure = factor(renal_failure, levels = c("No", "Yes", "Unknown"))
  )

cat("Clinical variables recoded.\n\n")

# =============================================================================
# STEP 2: Recode pre-PCI medications
# =============================================================================

cat("Step 2: Recoding pre-PCI medications...\n")

target_population <- target_population %>%
  mutate(
    prepci_aspirin = case_when(
      prepci_aspirin == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_aspirin = factor(prepci_aspirin, levels = c("No", "Yes")),
    
    prepci_clopidogrel = case_when(
      prepci_clopidogrel == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_clopidogrel = factor(prepci_clopidogrel, levels = c("No", "Yes")),
    
    prepci_ticagrelor = case_when(
      prepci_ticagrelor == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_ticagrelor = factor(prepci_ticagrelor, levels = c("No", "Yes")),
    
    prepci_prasugrel = case_when(
      prepci_prasugrel == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_prasugrel = factor(prepci_prasugrel, levels = c("No", "Yes")),
    
    prepci_heparin = case_when(
      prepci_heparin == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_heparin = factor(prepci_heparin, levels = c("No", "Yes")),
    
    prepci_bivalirudin = case_when(
      prepci_bivalirudin == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_bivalirudin = factor(prepci_bivalirudin, levels = c("No", "Yes")),
    
    prepci_warfarin = case_when(
      prepci_warfarin == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_warfarin = factor(prepci_warfarin, levels = c("No", "Yes")),
    
    prepci_fibrinolytic = case_when(
      prepci_fibrinolytic == "1" ~ "Yes", TRUE ~ "No"
    ),
    prepci_fibrinolytic = factor(prepci_fibrinolytic, levels = c("No", "Yes"))
  )

cat("Pre-PCI medications recoded.\n\n")

# =============================================================================
# STEP 3: Relabel SES variables
# =============================================================================

cat("Step 3: Relabeling SES variables...\n")

target_population <- target_population %>%
  mutate(
    birthcountry_cat = case_when(
      birthcountry_cat == 1 ~ "Sweden",
      birthcountry_cat == 2 ~ "Other Nordic",
      birthcountry_cat == 3 ~ "EU/Europe",
      birthcountry_cat %in% c(4, 5) ~ "Americas/Africa/Asia/Other",
      is.na(birthcountry_cat)       ~ "Missing",
      TRUE                          ~ "Missing"
    ),
    birthcountry_cat = factor(birthcountry_cat,
                              levels = c("Sweden", "Other Nordic", "EU/Europe",
                                         "Americas/Africa/Asia/Other", "Missing")),
    
    migration_3ybefore = case_when(
      migration_3ybefore == 1 ~ "Yes", TRUE ~ "No"
    ),
    migration_3ybefore = factor(migration_3ybefore, levels = c("No", "Yes")),
    
    region = case_when(
      baseline_swedish_region == 1                            ~ "Stockholm",
      baseline_swedish_region %in% c(2, 3)                   ~ "Uppsala/Linkoping",
      baseline_swedish_region == 4                           ~ "Lund/Malmo",
      baseline_swedish_region == 5                           ~ "Gothenburg",
      is.na(baseline_swedish_region)                         ~ "Missing",
      TRUE                                                   ~ "Missing"
    ),
    region = factor(region, levels = c("Stockholm", "Uppsala/Linkoping",
                                       "Lund/Malmo", "Gothenburg", "Missing")),
    
    degurba_code = case_when(
      is.na(baseline_degurba_code)    ~ "Missing",
      baseline_degurba_code == 1      ~ "Cities",
      baseline_degurba_code == 2      ~ "Towns/suburbs",
      baseline_degurba_code == 3      ~ "Rural areas",
      TRUE                            ~ as.character(baseline_degurba_code)
    ),
    degurba_code = factor(degurba_code,
                          levels = c("Cities", "Towns/suburbs", "Rural areas", "Missing")),
    
    moved_year_before = case_when(
      moved_year_before == 1 ~ "Yes",
      moved_year_before == 0 ~ "No",
      TRUE                   ~ "Missing"
    ),
    moved_year_before = factor(moved_year_before, levels = c("No", "Yes", "Missing")),
    
    family_position = case_when(
      baseline_family_position %in% c(3, 4) | is.na(baseline_family_position) ~ "Single",
      baseline_family_position == 1                                            ~ "Married/cohabiting",
      baseline_family_position == 2                                            ~ "Single w.children",
      baseline_family_position == 5                                            ~ "Unknown",
      TRUE                                                                     ~ "Single"
    ),
    family_position = factor(family_position,
                             levels = c("Married/cohabiting", "Single w.children",
                                        "Single", "Unknown")),
    
    separation_widow_3y = case_when(
      separation_widow_3y == 1 ~ "Yes", TRUE ~ "No"
    ),
    separation_widow_3y = factor(separation_widow_3y, levels = c("No", "Yes")),
    
    education_length = case_when(
      baseline_education_length == 1        ~ "Primary \u29c49yrs",
      baseline_education_length == 2        ~ "Secondary \u29c43yrs",
      baseline_education_length == 3        ~ "Tertiary",
      is.na(baseline_education_length)      ~ "Missing",
      TRUE                                  ~ "Missing"
    ),
    education_length = factor(education_length,
                              levels = c("Primary \u29c49yrs", "Secondary \u29c43yrs",
                                         "Tertiary", "Missing")),
    
    education_area = case_when(
      baseline_education_cat == 1        ~ "General",
      baseline_education_cat == 2        ~ "Social science",
      baseline_education_cat == 3        ~ "Science/Tech/Health",
      baseline_education_cat == 4        ~ "Services",
      is.na(baseline_education_cat)      ~ "Missing",
      TRUE                               ~ "Missing"
    ),
    education_area = factor(education_area,
                            levels = c("General", "Social science",
                                       "Science/Tech/Health", "Services", "Missing")),
    
    occupational_status = case_when(
      baseline_yseg_cat == 1             ~ "Senior Officials/Managers/Professionals",
      baseline_yseg_cat == 2             ~ "Small Business Owners/Farmers",
      baseline_yseg_cat == 3             ~ "Supervisors/Technicians/Skilled Manual",
      baseline_yseg_cat == 4             ~ "Office/Trade/Services/Care",
      baseline_yseg_cat == 5             ~ "Other Manual Worker",
      baseline_yseg_cat == 6             ~ "No Occupational Status/Unemployed",
      is.na(baseline_yseg_cat)           ~ "Missing",
      TRUE                               ~ "Missing"
    ),
    occupational_status = factor(occupational_status,
                                 levels = c("Senior Officials/Managers/Professionals",
                                            "Small Business Owners/Farmers",
                                            "Supervisors/Technicians/Skilled Manual",
                                            "Office/Trade/Services/Care",
                                            "Other Manual Worker",
                                            "No Occupational Status/Unemployed",
                                            "Missing")),
    
    unemployed_3ybefore = case_when(
      unemployed_3ybefore == 1 ~ "Yes", TRUE ~ "No"
    ),
    unemployed_3ybefore = factor(unemployed_3ybefore, levels = c("No", "Yes")),
    
    income_disposable_q5 = case_when(
      income_disposable_3ymean_q5 == 1 | is.na(income_disposable_3ymean_q5) ~ "Q1 (lowest)",
      income_disposable_3ymean_q5 == 2 ~ "Q2",
      income_disposable_3ymean_q5 == 3 ~ "Q3",
      income_disposable_3ymean_q5 == 4 ~ "Q4",
      income_disposable_3ymean_q5 == 5 ~ "Q5 (highest)",
      TRUE                             ~ "Q1 (lowest)"
    ),
    income_disposable_q5 = factor(income_disposable_q5,
                                  levels = c("Q1 (lowest)", "Q2", "Q3", "Q4", "Q5 (highest)"))
  )

cat("SES variables relabeled.\n\n")

# =============================================================================
# STEP 3B: Check sparse categories in observational data
# =============================================================================

cat("Step 3B: Checking sparse categories in observational data...\n")

temp_covariates <- c(
  "stemi", "female", "age", "bmi_cat", "smoking_status",
  "diabetes", "hypertension", "hyperlipidemia",
  "previous_mi", "previous_pci", "previous_cabg", "previous_stroke",
  "cpr_before_hospital", "killip_class",
  "systolic_bp_cat", "diastolic_bp_cat", "serum_creatinine_cat",
  "femoral_access", "angiography_finding", "renal_failure",
  "prepci_aspirin", "prepci_clopidogrel", "prepci_ticagrelor",
  "prepci_prasugrel", "prepci_heparin", "prepci_warfarin",
  "birthcountry_cat", "region", "family_position",
  "education_area", "education_length", "occupational_status",
  "unemployed_3ybefore", "income_disposable_q5"
)

check_sparse <- function(data, enrol_val, treat_val) {
  data %>%
    filter(enrol == enrol_val, treat_received == treat_val) %>%
    select(all_of(temp_covariates)) %>%
    select(where(is.factor)) %>%
    pivot_longer(everything()) %>%
    count(name, value) %>%
    filter(n < 30) %>%
    arrange(name, n)
}

cat("\nControl sparse categories (<30 obs):\n")
print(check_sparse(target_population, "Not randomised", "Control"))
cat("\nIntervention sparse categories (<30 obs):\n")
print(check_sparse(target_population, "Not randomised", "Intervention"))

rm(temp_covariates)
cat("Sparse category check complete.\n\n")

# =============================================================================
# STEP 3C: Collapse rare and sparse categories
# =============================================================================

cat("Step 3C: Collapsing rare categories...\n")

target_population <- target_population %>%
  mutate(
    # Comorbidities: collapse Unknown to No
    across(c(diabetes, hypertension, hyperlipidemia, previous_mi,
             previous_stroke, previous_cabg),
           ~ factor(ifelse(. == "Yes", "Yes", "No"), levels = c("No", "Yes"))),
    
    # BMI: collapse Underweight into Normal
    bmi_cat = factor(case_when(
      bmi_cat == "Underweight" ~ "Normal",
      TRUE                     ~ as.character(bmi_cat)
    ), levels = c("Normal", "Overweight", "Obesity", "Unknown")),
    
    # Systolic BP: collapse Missing into lowest category
    systolic_bp_cat = factor(case_when(
      systolic_bp_cat == "Missing" ~ "30 \u2264 Sys BP < 120",
      TRUE                        ~ as.character(systolic_bp_cat)
    ), levels = c("30 \u2264 Sys BP < 120", "120 \u2264 Sys BP < 140",
                  "140 \u2264 Sys BP < 160", "Sys BP \u2265 160")),
    
    # Medications: collapse Unknown to No
    across(c(prepci_aspirin, prepci_clopidogrel, prepci_ticagrelor, prepci_heparin),
           ~ factor(ifelse(. == "Yes", "Yes", "No"), levels = c("No", "Yes"))),
    
    # Education length: collapse Missing into Primary
    education_length = factor(case_when(
      education_length == "Missing" ~ "Primary \u29c49yrs",
      TRUE                         ~ as.character(education_length)
    ), levels = c("Primary \u29c49yrs", "Secondary \u29c43yrs", "Tertiary")),
    
    # Occupational status: collapse Missing into No Occupational Status/Unemployed
    occupational_status = factor(case_when(
      occupational_status == "Missing" ~ "No Occupational Status/Unemployed",
      TRUE                             ~ as.character(occupational_status)
    ), levels = c("Senior Officials/Managers/Professionals",
                  "Small Business Owners/Farmers",
                  "Supervisors/Technicians/Skilled Manual",
                  "Office/Trade/Services/Care",
                  "Other Manual Worker",
                  "No Occupational Status/Unemployed"))
  )

cat("Rare categories collapsed.\n\n")

# =============================================================================
# STEP 4: Define covariates
# =============================================================================

covariates <- c(
  "stemi", "female", "age", "bmi_cat", "smoking_status",
  "diabetes", "hypertension", "hyperlipidemia",
  "previous_mi", "previous_pci", "previous_cabg", "previous_stroke",
  "killip_class", "systolic_bp_cat", "diastolic_bp_cat",
  "serum_creatinine_cat", "femoral_access", "angiography_finding",
  "renal_failure", "prepci_aspirin", "prepci_clopidogrel",
  "prepci_ticagrelor", "prepci_heparin",
  "birthcountry_cat", "region", "family_position",
  "education_area", "education_length", "occupational_status",
  "unemployed_3ybefore", "income_disposable_q5"
)

cat("Covariates defined:", length(covariates), "total\n")
cat("Excluded: cpr_before_hospital, prepci_prasugrel, prepci_warfarin\n\n")

# =============================================================================
# STEP 5: Select final variables
# =============================================================================

cat("Step 5: Selecting final variables...\n")

target_population <- target_population %>%
  select(lopnr, enrol, treat_received, primary_ev180, primary_time,
         all_of(covariates))

cat("Final dataset: ", nrow(target_population), "observations,",
    ncol(target_population), "variables\n\n")

# =============================================================================
# FINAL CHECK
# =============================================================================

na_check <- target_population %>%
  select(all_of(covariates)) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_na") %>%
  filter(n_na > 0)

if (nrow(na_check) > 0) {
  cat("Variables with NA values:\n")
  print(na_check)
} else {
  cat("No NA values in covariates.\n")
}

# Clean up environment
rm(list = setdiff(ls(), c("target_population", "covariates")))

cat("\n=== Preparation complete ===\n")
cat("target_population ready:", nrow(target_population), "observations\n")
cat("covariates ready:", length(covariates), "variables\n")