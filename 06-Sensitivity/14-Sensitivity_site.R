# =============================================================================
# 14-Sensitivity_site.R
#
# Sensitivity analyses on site (Oxford Street vs Hyde Park) for metabolomics:
#   1. site as the predictor instead of the measured exposures
#   2. exposure models additionally adjusted for site
#   3. exposure models stratified by site (one residual variance per time point)
# Models that give a warning or an error are set to NA.
#
# Input:  New_data/metabolites_denoised_295.rds, New_data/covars_expo_295.rds
#         Results/MVN/metab_mvn_pval_beta_list.rds, Results/MVN/sign.metab.unique.rds
#         Results/metab_annotation_1309.csv
# Output: Results/MVN/metab/Sensitivity/Site_metab_mvn.rds
#         Results/MVN/metab/Sensitivity/<exposure>_site_metab_mvn.rds
#         Results/MVN/metab/Sensitivity/<exposure>_{ox,hp}_metab_mvn.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(nlme)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)

source("Scripts/functions.R")

metabolomics <- readRDS("New_data/metabolites_denoised_295.rds")
covars_expo <- readRDS("New_data/covars_expo_295.rds")
metab.list <- colnames(metabolomics)

dir.create("Results/MVN/metab/Sensitivity", showWarnings = FALSE, recursive = TRUE)

# Main results (08-Analysis_univariate.R)
main <- readRDS("Results/MVN/metab_mvn_pval_beta_list.rds")
metab.pval.df <- main$pval
metab.beta.df <- main$beta
sign.metab.unique <- readRDS("Results/MVN/sign.metab.unique.rds")
expo.names <- c("PM25", "PM10", "NO2", "BC", "PCNT")


# 1. Site instead of exposure ---------------------------------------------
site.res <- run_mvn(metabolomics, covars_expo, rhs = "Site + Age + sex + BMI + Group", on_failure = "NA_warning")
saveRDS(site.res, "Results/MVN/metab/Sensitivity/Site_metab_mvn.rds")

metab.pval.df$Site <- site.res$pval[, "SiteOxford Street"]
metab.beta.df$Site <- site.res$beta[, "SiteOxford Street"]

# p-values: site vs exposure (red: selected in the main analysis)
par(mfrow = c(2, 3))
for (expo in expo.names) {
  max.pval <- ceiling(max(-log10(metab.pval.df$Site), -log10(metab.pval.df[, expo]), na.rm = TRUE))
  point.col <- ifelse(metab.pval.df[, expo] < alpha_metab, 2, 1)
  plot(-log10(metab.pval.df$Site), -log10(metab.pval.df[, expo]),
    xlab = "-log10 p-value for Site (Sensitivity)",
    ylab = paste("-log10 p-value for", expo, "(Main)"), main = expo, pch = 1, cex = 0.5,
    xlim = c(1, max.pval), ylim = c(1, max.pval),
    col = point.col
  )
  cor_coef <- cor(-log10(metab.pval.df$Site), -log10(metab.pval.df[, expo]), use = "complete.obs")
  text(x = 2, y = max.pval - 1, labels = paste("Cor:", round(cor_coef, 2)), pos = 4, cex = 1)
  abline(h = -log10(alpha_metab), lty = 2)
  abline(v = -log10(alpha_metab), lty = 2)
}

# Coefficients: site vs exposure
par(mfrow = c(2, 3))
for (expo in expo.names) {
  plot(metab.beta.df$Site, metab.beta.df[[expo]], xlab = "Site", ylab = expo)
  cor_coef <- cor(metab.beta.df[[expo]], metab.beta.df$Site, use = "complete.obs")
  title(paste("Corr:", round(cor_coef, 3)), cex = 0.8)
}

# Features selected with site and with each exposure
metab.sign.site.tf.df <- metab.pval.df %>%
  filter(metab %in% c(sign.metab.unique, metab.list[which(metab.pval.df$Site < alpha_metab)])) %>%
  tibble::column_to_rownames("metab") %>%
  mutate(across(everything(), ~ .x < alpha_metab))
upset(metab.sign.site.tf.df, c("PM25", "PM10", "PCNT", "BC", "NO2", "Site"),
  base_annotations = list("Intersection size" = intersection_size(text = list(size = 3))),
  set_sizes = (upset_set_size() +
    geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = "count", col = "black", size = 4) +
    expand_limits(y = 85)),
  themes = upset_default_themes(text = element_text(size = 15)),
  name = ""
)

# Volcano plot for site
par(mfrow = c(1, 1))
plot(metab.beta.df$Site, -log10(metab.pval.df$Site),
  main = "Volcano plot for site",
  xlab = "Beta", ylab = "-log10(p-value)"
)
abline(h = -log10(alpha_metab), lty = 2)

metab.annotation <- read.csv("Results/metab_annotation_1309.csv", fileEncoding = "UTF-8-BOM")
metab.annotation %>% filter(metab == "644.5358@7.188575")


# 2. Exposure models adjusted for site ------------------------------------
mvn_by_exposure(metabolomics, covars_expo,
  rhs = "Expo.delta + Age + sex + BMI + Group + Site",
  out_file = "Results/MVN/metab/Sensitivity/%s_site_metab_mvn.rds",
  on_failure = "NA_warning"
)

# Models for the PMs had convergence issues: compare for NO2
NO2_site_metab_mvn <- readRDS("Results/MVN/metab/Sensitivity/NO2_site_metab_mvn.rds")
plot(-log10(metab.pval.df$NO2), -log10(NO2_site_metab_mvn$pval[, "Expo.delta"]),
  xlab = "Base model", ylab = "Adjusted for site",
  xlim = c(0, 15.5), ylim = c(0, 15.5)
)


# 3. Exposure models stratified by site -----------------------------------
for (Loc in c("ox", "hp")) {
  in_site <- covars_expo$Loc == Loc
  mvn_by_exposure(metabolomics[in_site, ], covars_expo[in_site, ],
    rhs = "Expo.delta + Age + sex + BMI + Group",
    out_file = paste0("Results/MVN/metab/Sensitivity/%s_", Loc, "_metab_mvn.rds"),
    var_group = "Time_Point",
    on_failure = "NA_warning"
  )
}

# p-values in each site
expo.files <- c(PM25 = "PM25", PM10 = "PM10", NO2 = "NO2", BC = "CBLK", PCNT = "PCNT")
stratified_pvals <- function(Loc) {
  data.frame(metab = metab.list, sapply(expo.files, function(e) {
    readRDS(sprintf("Results/MVN/metab/Sensitivity/%s_%s_metab_mvn.rds", e, Loc))$pval[, "Expo.delta"]
  }))
}
metab.hp.pval.df <- stratified_pvals("hp")
metab.ox.pval.df <- stratified_pvals("ox")

bind_rows(
  metab.hp.pval.df %>% mutate(Loc = "hp"),
  metab.ox.pval.df %>% mutate(Loc = "ox")
) %>%
  pivot_longer(PM25:PCNT, names_to = "Expo") %>%
  pivot_wider(values_from = value, names_from = Loc) %>%
  ggplot(aes(-log10(ox), -log10(hp))) +
  geom_point(size = 0.1) +
  facet_wrap(~Expo) +
  theme_bw()

# Hyde Park vs Oxford Street (red: selected in the main analysis)
par(mfrow = c(2, 3))
for (expo in expo.names) {
  max.pval <- ceiling(max(-log10(metab.hp.pval.df[, expo]), -log10(metab.ox.pval.df[, expo]), na.rm = TRUE))
  point.col <- ifelse(metab.pval.df[, expo] < alpha_metab, 2, 1)
  plot(-log10(metab.hp.pval.df[, expo]), -log10(metab.ox.pval.df[, expo]),
    xlab = "Hyde Park", ylab = "Oxford Street", main = expo, pch = 1, cex = 0.5,
    xlim = c(1, max.pval), ylim = c(1, max.pval),
    col = point.col
  )
  cor_coef <- cor(-log10(metab.hp.pval.df[, expo]), -log10(metab.ox.pval.df[, expo]), use = "complete.obs")
  text(x = 2, y = max.pval - 1, labels = paste("Cor:", round(cor_coef, 2)), pos = 4, cex = 1)
  abline(h = -log10(alpha_metab), lty = 2)
  abline(v = -log10(alpha_metab), lty = 2)
}

# Exposure levels by site
covars_expo %>%
  group_by(ID, Loc) %>%
  filter(Per != "before") %>%
  slice(1) %>%
  pivot_longer(CBLK_delta:PCNT_delta, names_to = "expo") %>%
  ggplot(aes(Site, value)) +
  geom_boxplot() +
  facet_wrap(~expo, scale = "free") +
  theme_bw()
