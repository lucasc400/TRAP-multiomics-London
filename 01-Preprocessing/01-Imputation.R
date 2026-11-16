# =============================================================================
# 01-Imputation.R
#
# Prepares the untargeted metabolomics data: links each sample to its subject,
# visit, site and time point (XP_ID), attaches the technical covariates (box
# and plate), log2-transforms intensities, removes features missing in 40% or
# more of the samples and imputes the remaining missing values with QRILC.
# Imputation is run with and without the five outlying samples identified by
# PCA (see 05-EDA.R).
#
# Input:  Data/OxFord_Street_2_Final_Export_Metabolomics_June_2015.xlsx
#         Data/Analysis_chart_2015.xlsx
#         Data/Exposomics sample spreadsheet- Oxford Street Study II.xlsx
# Output: New_data/metab_tech_confounders.rds
#         New_data/metabolites_imputed_full.rds  (all samples)
#         New_data/metabolites_imputed.rds       (outliers removed)
#         New_data/metabolites_logged_not_imputed.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(readxl)
library(imputeLCMD)
library(tidyr)
library(dplyr)
library(stringr)

metab <- read_excel("Data/OxFord_Street_2_Final_Export_Metabolomics_June_2015.xlsx", sheet = 1, skip = 3)


# Prepare data ------------------------------------------------------------
metab.complete <- metab[complete.cases(metab), ]
metab.data <- metab.complete[7:ncol(metab.complete)]
metab.covars <- metab.complete[1:6]
metab.data <- apply(metab.data, 2, function(x) ifelse(x == 1, NA, x)) # intensity of 1 = not detected

# Sample information: subject, visit, site and time point
metab.covars <- metab.covars %>%
  rename(Subject_ID = `Subject ID`) %>%
  mutate(
    Site = case_when(str_detect(Site, "Hyde") ~ "Hyde Park", TRUE ~ "Oxford Street"),
    Subject_ID = str_replace(Subject_ID, "O", "0")
  ) %>%
  mutate(Site = case_when(Visit == "V2" & Subject_ID == "RS010" ~ "Hyde Park", TRUE ~ Site))


# Technical covariates ----------------------------------------------------
analysis_chart <- read_excel("Data/Analysis_chart_2015.xlsx", sheet = 1, skip = 1)
sample_info.1 <- read_excel("Data/Exposomics sample spreadsheet- Oxford Street Study II.xlsx", sheet = 1)
sample_info.2 <- read_excel("Data/Exposomics sample spreadsheet- Oxford Street Study II.xlsx", sheet = 2)
sample_info.1 <- sample_info.1[-21, ]

# Site visited at each visit by each subject
site_lookup_table <- as.data.frame(mapply(c, sample_info.1[, 1:11], sample_info.1[, 12:22], sample_info.2[, 1:11], sample_info.2[, 12:22]))
colnames(site_lookup_table) <- str_replace(colnames(site_lookup_table), "\\.\\.\\..", "")
site_lookup_table$Visit <- c(
  rep("V1", nrow(sample_info.1)), rep("V2", nrow(sample_info.1)),
  rep("V1", nrow(sample_info.2)), rep("V2", nrow(sample_info.2))
)
site_lookup_table <- site_lookup_table %>%
  mutate(Site = case_when(str_detect(Site, "Hyde") ~ "Hyde Park", TRUE ~ "Oxford Street")) %>%
  dplyr::select(`Subject ID`, Visit, Site)

# Box and plate of each sample, and sample ID: XP_ID = <subject><site 1/2><time point 1/2/3>
tech_covars.long <- pivot_longer(analysis_chart, cols = 6:8, names_to = "Time_Point")
tech_covars.long.clean <- tech_covars.long %>%
  mutate(Time_Point = case_when(
    Time_Point == "Time Point 1 (8 Am)" ~ "T1 (8AM)",
    Time_Point == "Time Point 2 (4 pm)" ~ "T2 (4PM)",
    Time_Point == "Time Point 3 (24 hrs)" ~ "T3 (24hrs)"
  )) %>%
  rename(
    Box_No = "Metabolomics – Lyon",
    Plate_No = `Plate No`,
    Visit = `...2`,
    Subject_ID = `Subject ID`
  ) %>%
  left_join(., site_lookup_table %>% rename(Subject_ID = "Subject ID"), by = c("Subject_ID", "Visit")) %>%
  mutate(
    Subject_ID = str_replace(Subject_ID, "O", "0"),
    ID = as.integer(str_replace(Subject_ID, "RS", "")),
    ID.1 = case_when(str_detect(Site, "Hyde") ~ 1, TRUE ~ 2),
    ID.2 = case_when(
      Time_Point == "T1 (8AM)" ~ 1,
      Time_Point == "T2 (4PM)" ~ 2,
      TRUE ~ 3
    ),
    XP_ID = paste0(ID, ID.1, ID.2)
  ) %>%
  dplyr::select(Subject_ID, Box_No, Plate_No, Time_Point, Visit, Site, ID, XP_ID)

# Subject 10 visited Hyde Park at both visits: fix the site and XP_ID of visit 2
tech_covars.long.clean <- tech_covars.long.clean %>%
  mutate(
    XP_ID = ifelse((Visit == "V2") & (ID == 10), as.integer(XP_ID) - 10, as.integer(XP_ID)),
    Site = ifelse((Visit == "V2") & (ID == 10), "Hyde Park", Site)
  )

# Remove samples with missing site
tech_covars.long.clean <- tech_covars.long.clean[complete.cases(tech_covars.long.clean), ]

# Add XP_ID to the metabolomics data
metab.covars.xp_id <- metab.covars %>%
  left_join(., tech_covars.long.clean, by = c("Subject_ID", "Visit", "Site", "Time_Point"))
rownames(metab.data) <- metab.covars.xp_id$XP_ID

# Align technical covariates with metabolomics data
tech_covars.long.clean <- as.data.frame(tech_covars.long.clean) # avoid tibble behaviour
rownames(tech_covars.long.clean) <- tech_covars.long.clean$XP_ID
tech_covars.long.clean.aligned <- tech_covars.long.clean[rownames(metab.data), ]


# Imputation --------------------------------------------------------------
metab.data.logged <- log2(metab.data)

# Keep features missing in at most 40% of the samples
pct_missing <- apply(metab.data, 2, function(x) sum(is.na(x)) / nrow(metab.data) * 100)
metab.60 <- metab.data.logged[, pct_missing <= 40]

# All samples
set.seed(400)
metab.60.imputed <- t(impute.QRILC(t(metab.60), tune.sigma = 1)[[1]])

# Without the five outlying samples (see 05-EDA.R)
outliers <- c("6011", "1511", "7623", "10123", "9721")
metab.60.no.outliers <- metab.60[!rownames(metab.60) %in% outliers, ]
set.seed(400)
metab.60.no.outliers.imputed <- t(impute.QRILC(t(metab.60.no.outliers), tune.sigma = 1)[[1]])


# Save data ---------------------------------------------------------------
dir.create("New_data", showWarnings = FALSE)
saveRDS(tech_covars.long.clean.aligned, "New_data/metab_tech_confounders.rds")
saveRDS(metab.60.imputed, "New_data/metabolites_imputed_full.rds")
saveRDS(metab.60.no.outliers.imputed, "New_data/metabolites_imputed.rds")
saveRDS(metab.data.logged, "New_data/metabolites_logged_not_imputed.rds")
