# =============================================================================
# 16-Sensitivity_outliers.R
#
# Sensitivity analysis on the five outlying metabolomic samples:
#   1. models on all 300 samples, including the outliers
#   2. models on 295 samples keeping the outliers and removing five randomly
#      chosen other samples (25 replicates), to compare the number of
#      selected features with the main analysis
# Models that give an error are set to NA.
#
# Input:  New_data/metabolites_denoised_300.rds, New_data/covars_expo_300.rds
# Output: Results/MVN/metab/Sensitivity/<exposure>_metab_res_300_.rds
#         Results/MVN/metab/Sensitivity/<exposure>_metab_res_remove_5_<replicate>.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(nlme)
library(dplyr)
library(tidyr)

source("Scripts/functions.R")

metabolomics <- readRDS("New_data/metabolites_denoised_300.rds")
covars_expo <- readRDS("New_data/covars_expo_300.rds")

dir.create("Results/MVN/metab/Sensitivity", showWarnings = FALSE, recursive = TRUE)

rhs <- "Expo.delta + Age + sex + BMI + Group"
outliers <- c("6011", "1511", "7623", "10123", "9721")
n_replicates <- 25


# 1. All 300 samples ------------------------------------------------------
mvn_by_exposure(metabolomics, covars_expo,
  rhs = rhs,
  out_file = "Results/MVN/metab/Sensitivity/%s_metab_res_300_.rds",
  on_failure = "NA"
)


# 2. Outliers kept, five other samples removed ----------------------------
for (i in 1:n_replicates) {
  set.seed(i)
  index.290 <- sample(which(!rownames(metabolomics) %in% outliers), 290)
  index.295 <- c(index.290, which(rownames(metabolomics) %in% outliers))
  mvn_by_exposure(metabolomics[index.295, ], covars_expo[index.295, ],
    rhs = rhs,
    out_file = paste0("Results/MVN/metab/Sensitivity/%s_metab_res_remove_5_", i, ".rds"),
    on_failure = "NA"
  )
}


# Results -----------------------------------------------------------------
expo.files <- c(NO2 = "NO2", BC = "CBLK", PM25 = "PM25", PM10 = "PM10", PCNT = "PCNT")
n_selected <- function(file) sum(readRDS(file)$pval[, "Expo.delta"] < alpha_metab, na.rm = TRUE)

# Number of selected features across replicates
sens.sign.number <- sapply(expo.files, function(e) {
  sapply(1:n_replicates, function(i) {
    n_selected(sprintf("Results/MVN/metab/Sensitivity/%s_metab_res_remove_5_%d.rds", e, i))
  })
})
as.data.frame(sens.sign.number) %>%
  pivot_longer(everything(), names_to = "Expo") %>%
  group_by(Expo) %>%
  summarise(
    lower_q = quantile(value, 0.25),
    median = median(value),
    upper_q = quantile(value, 0.75)
  )

# Number of selected features with all 300 samples
sapply(expo.files, function(e) n_selected(sprintf("Results/MVN/metab/Sensitivity/%s_metab_res_300_.rds", e)))
