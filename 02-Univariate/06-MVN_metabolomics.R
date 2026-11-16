# =============================================================================
# 06-MVN_metabolomics.R
#
# Metabolome-wide association study. Each denoised metabolic feature is
# regressed on each TRAP exposure (change from annual average) with a
# multivariate normal model (see run_mvn() in functions.R), adjusting for age,
# sex, BMI and health group. The five exposures are run in parallel.
#
# Input:  New_data/metabolites_denoised_295.rds, New_data/covars_expo_295.rds
# Output: Results/MVN/metab/<exposure>_metab_mvn.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(nlme)

source("Scripts/functions.R")

metabolomics <- readRDS("New_data/metabolites_denoised_295.rds")
covars_expo <- readRDS("New_data/covars_expo_295.rds")

dir.create("Results/MVN/metab", showWarnings = FALSE, recursive = TRUE)

mvn_by_exposure(metabolomics, covars_expo,
  rhs = "Expo.delta + Age + sex + BMI + Group",
  out_file = "Results/MVN/metab/%s_metab_mvn.rds"
)
