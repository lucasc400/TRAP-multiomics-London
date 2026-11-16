# =============================================================================
# 03-Denoising.R
#
# Removes technical variation from the omics data. For each feature, a linear
# mixed model is fitted with the technical covariates as random intercepts
# (box and plate for metabolomics; array and batch for transcriptomics), and
# the estimated random intercepts are subtracted from the measured levels.
#
# Input:  New_data/metabolites_{300,295}.rds, New_data/metab_tech_confounders_{300,295}.rds
#         New_data/mRNA_252.rds, New_data/mRNA_tech_252.rds (02-Data_cleaning.R)
# Output: New_data/metabolites_denoised_{300,295}.rds, New_data/mRNA_denoised_252.rds
#         Results/Denoising/metab.lmer.295.rds, Results/Denoising/mRNA.lmer.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(lme4)

source("Scripts/functions.R")

dir.create("Results/Denoising", showWarnings = FALSE, recursive = TRUE)


# Metabolomics, all samples -----------------------------------------------
metabolites.300 <- readRDS("New_data/metabolites_300.rds")
metab_tech_confounders.300 <- readRDS("New_data/metab_tech_confounders_300.rds")

metab.lmer <- Multilmer(metabolites.300, metab_tech_confounders.300, random = c("Box_No", "Plate_No"), residuals = TRUE, tol = 1e-4)
saveRDS(metab.lmer$residuals, "New_data/metabolites_denoised_300.rds")


# Metabolomics, outliers removed ------------------------------------------
metabolites.295 <- readRDS("New_data/metabolites_295.rds")
metab_tech_confounders.295 <- readRDS("New_data/metab_tech_confounders_295.rds")

metab.lmer.295 <- Multilmer(metabolites.295, metab_tech_confounders.295, random = c("Box_No", "Plate_No"), residuals = TRUE, tol = 1e-4)
saveRDS(metab.lmer.295, "Results/Denoising/metab.lmer.295.rds")
saveRDS(metab.lmer.295$residuals, "New_data/metabolites_denoised_295.rds")


# Transcriptomics ---------------------------------------------------------
mRNA.252 <- readRDS("New_data/mRNA_252.rds")
mRNA.tech.252 <- readRDS("New_data/mRNA_tech_252.rds")

mRNA.lmer <- Multilmer(mRNA.252, mRNA.tech.252, random = c("array", "batch"), residuals = TRUE, tol = 1e-4)
saveRDS(mRNA.lmer, "Results/Denoising/mRNA.lmer.rds")
saveRDS(mRNA.lmer$residuals, "New_data/mRNA_denoised_252.rds")
