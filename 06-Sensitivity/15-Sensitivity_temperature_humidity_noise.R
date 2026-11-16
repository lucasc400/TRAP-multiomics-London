# =============================================================================
# 15-Sensitivity_temperature_humidity_noise.R
#
# Sensitivity analysis on meteorological conditions and noise: the exposure
# models for metabolomics are additionally adjusted for
#   M1: temperature
#   M2: temperature and relative humidity
#   M3: temperature, relative humidity and noise
# and compared with the main models. The transcriptomics models adjusted for
# temperature and humidity are fitted in 07-MVN_transcriptomics.R.
#
# Input:  New_data/metabolites_denoised_295.rds, New_data/covars_expo_295.rds
#         Results/MVN/{metab,mRNA}_mvn_pval_beta_list.rds (08-Analysis_univariate.R)
#         Results/MVN/mRNA/Sensitivity/*_mRNA_mvn_temp_rhum.rds (07-MVN_transcriptomics.R)
#         Results/metab_annotation_1309.csv
# Output: Results/MVN/metab/Sensitivity/<exposure>_metab_mvn_{temp,temp_rhum,temp_rh_noise}.rds
#         Results/Figures/metab_temp_humd_sensitivity_ent.pdf
#         Results/Figures/mRNA_temp_humd_sensitivity.pdf
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(nlme)
library(dplyr)
library(pheatmap)

source("Scripts/functions.R")

metabolomics <- readRDS("New_data/metabolites_denoised_295.rds")
covars_expo <- readRDS("New_data/covars_expo_295.rds")
metab.list <- colnames(metabolomics)

dir.create("Results/MVN/metab/Sensitivity", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Figures", showWarnings = FALSE, recursive = TRUE)

# Correlation between exposures, temperature and humidity
covars_expo_temp_rhum <- covars_expo[, c("RHUM", "temp", "CBLK_delta", "NO2_delta", "PM10_delta", "PM25_delta", "PCNT_delta")]
covars_expo_temp_rhum <- covars_expo_temp_rhum[complete.cases(covars_expo_temp_rhum), ]
pheatmap(cor(covars_expo_temp_rhum, method = "spearman"),
  treeheight_row = FALSE, treeheight_col = FALSE, display_numbers = TRUE
)


# Models ------------------------------------------------------------------
models <- c(
  temp = "Expo.delta + Age + sex + BMI + Group + temp",
  temp_rhum = "Expo.delta + Age + sex + BMI + Group + temp + RHUM",
  temp_rh_noise = "Expo.delta + Age + sex + BMI + Group + temp + RHUM + NOIS"
)
for (m in names(models)) {
  mvn_by_exposure(metabolomics, covars_expo,
    rhs = models[[m]],
    out_file = paste0("Results/MVN/metab/Sensitivity/%s_metab_mvn_", m, ".rds")
  )
}


# Comparison with the main models: metabolomics ---------------------------
expo.files <- c(PM25 = "PM25", PM10 = "PM10", NO2 = "NO2", BC = "CBLK", PCNT = "PCNT")
read_pvals <- function(omic, suffix) {
  sapply(expo.files, function(e) {
    readRDS(sprintf("Results/MVN/%s/Sensitivity/%s_%s_mvn_%s.rds", omic, e, omic, suffix))$pval[, "Expo.delta"]
  })
}
metab.pval.df <- readRDS("Results/MVN/metab_mvn_pval_beta_list.rds")$pval
metab.sens.pval <- lapply(setNames(names(models), names(models)), function(m) read_pvals("metab", m))

# p-values, main model vs adjusted for temperature and humidity (M2), with
# the Bonferroni (dotted) and ENT (dashed) thresholds
pdf("Results/Figures/metab_temp_humd_sensitivity_ent.pdf", width = 9, height = 6)
par(mar = c(6, 5, 6, 2), mfrow = c(2, 3))
for (expo in names(expo.files)) {
  plot(-log10(metab.pval.df[[expo]]), -log10(metab.sens.pval$temp_rhum[, expo]),
    xlab = paste(expo, "(base model)"), ylab = paste(expo, "(+ temp + humidity)")
  )
  cor_coef <- cor(-log10(metab.pval.df[[expo]]), -log10(metab.sens.pval$temp_rhum[, expo]), use = "complete.obs")
  title(paste("Corr:", round(cor_coef, 3)), cex = 0.8)
  abline(v = -log10(alpha_metab), lty = 3)
  abline(h = -log10(alpha_metab), lty = 3)
  abline(v = -log10(alpha_metab_ent), lty = 2, col = "red")
  abline(h = -log10(alpha_metab_ent), lty = 2, col = "red")
}
mtext("Metabolomics", side = 3, line = -2, outer = TRUE)
dev.off()

# Features selected in the main model (rows) and in each adjusted model (columns)
for (m in names(models)) {
  for (expo in names(expo.files)) {
    cat("\n", m, "-", expo, "\n")
    print(table(main = metab.pval.df[[expo]] < alpha_metab, adjusted = metab.sens.pval[[m]][, expo] < alpha_metab))
  }
}

# Features selected with NO2 that are not selected after adjusting for temperature (M1)
metab.annotation <- read.csv("Results/metab_annotation_1309.csv", fileEncoding = "UTF-8-BOM")
metab.annotation %>%
  filter(metab %in% metab.list[metab.pval.df$NO2 < alpha_metab & metab.sens.pval$temp[, "NO2"] >= alpha_metab])


# Comparison with the main models: transcriptomics ------------------------
mRNA.pval.df <- readRDS("Results/MVN/mRNA_mvn_pval_beta_list.rds")$pval
mRNA.M2.pval <- read_pvals("mRNA", "temp_rhum")

pdf("Results/Figures/mRNA_temp_humd_sensitivity.pdf", width = 9, height = 6)
par(mar = c(6, 5, 6, 2), mfrow = c(2, 3))
for (expo in names(expo.files)) {
  plot(-log10(mRNA.pval.df[[expo]]), -log10(mRNA.M2.pval[, expo]),
    xlab = paste(expo, "(Base model)"), ylab = paste(expo, "(+ temp + humidity)")
  )
  cor_coef <- cor(-log10(mRNA.pval.df[[expo]]), -log10(mRNA.M2.pval[, expo]), use = "complete.obs")
  title(paste("Corr:", round(cor_coef, 3)), cex = 0.8)
  abline(v = -log10(alpha_mrna), lty = 2, col = "red")
  abline(h = -log10(alpha_mrna), lty = 2, col = "red")
}
mtext("Transcriptomics", side = 3, line = -2, outer = TRUE)
dev.off()
