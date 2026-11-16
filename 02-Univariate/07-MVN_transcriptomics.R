# =============================================================================
# 07-MVN_transcriptomics.R
#
# Transcriptome-wide association study, using the same multivariate normal
# models as for metabolomics (see 06-MVN_metabolomics.R). A sensitivity
# analysis additionally adjusts for temperature and relative humidity.
#
# Input:  New_data/mRNA_denoised_252.rds, New_data/covars_expo_252.rds
# Output: Results/MVN/mRNA/<exposure>_mRNA_mvn.rds
#         Results/MVN/mRNA/Sensitivity/<exposure>_mRNA_mvn_temp_rhum.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(nlme)

source("Scripts/functions.R")

mRNA <- readRDS("New_data/mRNA_denoised_252.rds")
covars_expo.252 <- readRDS("New_data/covars_expo_252.rds")

dir.create("Results/MVN/mRNA/Sensitivity", showWarnings = FALSE, recursive = TRUE)

# Main models
mvn_by_exposure(mRNA, covars_expo.252,
  rhs = "Expo.delta + Age + sex + BMI + Group",
  out_file = "Results/MVN/mRNA/%s_mRNA_mvn.rds"
)

# Sensitivity analysis: adjusting for temperature and relative humidity
mvn_by_exposure(mRNA, covars_expo.252,
  rhs = "Expo.delta + Age + sex + BMI + Group + RHUM + temp",
  out_file = "Results/MVN/mRNA/Sensitivity/%s_mRNA_mvn_temp_rhum.rds"
)
