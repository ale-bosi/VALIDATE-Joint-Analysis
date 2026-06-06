# VALIDATE-Joint-Analysis
Joint analysis of the VALIDATE-SWEDEHEART trial

## Overview

This repository contains the R code for a joint analysis of the VALIDATE-SWEDEHEART trial, a registry-based randomized controlled trial nested in SWEDEHEART evaluating bivalirudin versus unfractionated heparin during percutaneous coronary intervention in patients with myocardial infarction.

The analysis combines data from trial participants and eligible non-participants to estimate the per-protocol effect of bivalirudin versus heparin on all-cause mortality in the full target population. Three identification strategies are considered — generalizability analysis, stratified joint analysis, and pooled joint analysis — each resting on different assumptions about the comparability of trial participants and non-participants. Treatment effects are estimated using standardization (g-formula), inverse probability weighting (IPW), and augmented inverse probability weighting (AIPW).
