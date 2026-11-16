# =============================================================================
# 02-Data_cleaning.R
#
# Cleans the transcriptomics data and its technical covariates (array and
# batch), computes the TRAP exposures and merges them with the covariates.
# Exposures are expressed as the change from the participant's annual average
# ("delta"): zero before the walk, and the walk measurement minus the annual
# average at 2 h and 24 h after the walk. Builds the analysis datasets:
#   - 300 metabolomic profiles (50 participants x 2 sites x 3 time points)
#   - 295 metabolomic profiles (five outlying samples removed)
#   - 252 transcriptomic profiles with complete covariates
#
# Input:  Data/OXF2_Tx_mRNA_procd_bkgrcorr_compld.txt, Data/OXF2_Tx_mRNA_LabCode.csv
#         Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt
#         Data/Exposures.csv, Data/All_participants_ID_group.xlsx
#         New_data/metabolites_imputed*.rds, New_data/metab_tech_confounders.rds (01-Imputation.R)
# Output: New_data/covars_expo_{300,295,252}.rds
#         New_data/metabolites_{300,295}.rds, New_data/metab_tech_confounders_{300,295}.rds
#         New_data/mRNA_{290,252}.rds, New_data/mRNA_tech_{290,252}.rds
#         New_data/expo.measured.100.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(readxl)
library(stringr)


# Read data ---------------------------------------------------------------
OXF2 <- read.csv("Data/OXF2_Tx_mRNA_procd_bkgrcorr_compld.txt", sep = "")
OXF2_Tx_mRNA_LabCode <- read.csv("Data/OXF2_Tx_mRNA_LabCode.csv")
Full_covar_add <- read.table("Data/Full_covar_including_additional_airpollution_measurements_ratio_offset_annual_dist_temp_caff.txt",
  sep = "\t", header = TRUE
)
expo.measured <- read.csv("Data/Exposures.csv")
metabolites.full <- readRDS("New_data/metabolites_imputed_full.rds")
metabolites.no.outliers <- readRDS("New_data/metabolites_imputed.rds")
metab.tech.confounders <- readRDS("New_data/metab_tech_confounders.rds")


# Transcriptomics ---------------------------------------------------------
colnames(OXF2)[which(colnames(OXF2) == "OXF2_RS026_Hy_24.1")] <- "OXF2_RS026_Hy_8" # mislabelled sample
OXF2 <- OXF2[2:nrow(OXF2), ]
mRNA <- as.data.frame(t(OXF2))

# Time point and sample ID (XP_ID = <subject><site 1/2><time point 1/2/3>)
mRNA$idfull <- rownames(mRNA)
mRNA$id <- as.numeric(substring(mRNA$idfull, 8, 10))
mRNA$Time_Point <- as.factor(ifelse(grepl("_4", mRNA$idfull), "T2 (4PM)",
  ifelse(grepl("_8", mRNA$idfull), "T1 (8AM)", "T3 (24hrs)")
))
numbers1 <- ifelse(grepl("Hy_8", mRNA$idfull), 11,
  ifelse(grepl("Hy_4", mRNA$idfull), 12,
    ifelse(grepl("Hy_24", mRNA$idfull), 13,
      ifelse(grepl("Ox_8", mRNA$idfull), 21,
        ifelse(grepl("Ox_4", mRNA$idfull), 22, 23)
      )
    )
  )
)
mRNA$XP_ID <- as.numeric(paste(mRNA$id, numbers1, sep = ""))
mRNA <- mRNA[, c(30924:30927, 1:30923)] # sample information first

# Technical covariates
mRNA.tech <- OXF2_Tx_mRNA_LabCode
mRNA.tech$id <- as.numeric(substring(mRNA.tech$original.code, 3, 5))
mRNA.tech <- mRNA.tech[!is.na(mRNA.tech$id), ]
numbers2 <- ifelse(grepl("1-8am", mRNA.tech$original.code), 11,
  ifelse(grepl("1-4pm", mRNA.tech$original.code), 12,
    ifelse(grepl("1-24h", mRNA.tech$original.code), 13,
      ifelse(grepl("2-8am", mRNA.tech$original.code), 21,
        ifelse(grepl("2-4pm", mRNA.tech$original.code), 22, 23)
      )
    )
  )
)
mRNA.tech$XP_ID <- as.numeric(paste(mRNA.tech$id, numbers2, sep = ""))
mRNA.tech <- mRNA.tech[, c(20:21, 1:19)]
mRNA.tech <- mRNA.tech[mRNA.tech$Array.ID != "array not okay", ] # failed arrays
mRNA.tech$array <- as.integer(substring(mRNA.tech$Array.ID, 1, 5))
mRNA.tech <- mRNA.tech[, c(1:5, 7:9, 11:12, 22)]

# Join transcripts with technical covariates
mRNA.full <- inner_join(mRNA.tech, mRNA, by = "XP_ID")


# Exposure data -----------------------------------------------------------
# Annual averages
expo.annual.300 <- Full_covar_add[, c("id", "CBLK_annual", "NO2_annual", "PM10_annual", "PM25_annual", "PCNT_annual", "Loc", "Per")]
expo.annual <- expo.annual.300[!duplicated(expo.annual.300$id), c("id", "CBLK_annual", "NO2_annual", "PM10_annual", "PM25_annual", "PCNT_annual")]
rownames(expo.annual) <- expo.annual$id

# Measurements during the walks
expo.measured.hp <- expo.measured[expo.measured$Site == "Hyde Park", c("id", "CBLK", "NO2", "PM10", "PM25", "PCNT")]
expo.measured.ox <- expo.measured[expo.measured$Site == "Oxford St", c("id", "CBLK", "NO2", "PM10", "PM25", "PCNT")]
rownames(expo.measured.hp) <- expo.measured.hp$id
rownames(expo.measured.ox) <- expo.measured.ox$id
expo.measured.hp <- expo.measured.hp[rownames(expo.annual), ]
expo.measured.ox <- expo.measured.ox[rownames(expo.annual), ]

# Change from the annual average
expo.delta.hp <- expo.measured.hp
expo.delta.ox <- expo.measured.ox
expo.delta.hp[, 2:6] <- expo.measured.hp[, 2:6] - expo.annual[, 2:6]
expo.delta.ox[, 2:6] <- expo.measured.ox[, 2:6] - expo.annual[, 2:6]
expo.delta <- rbind(expo.delta.hp, expo.delta.ox)
colnames(expo.delta) <- c("id", "CBLK_delta", "NO2_delta", "PM10_delta", "PM25_delta", "PCNT_delta")
expo.delta$Site <- c(rep("Hyde Park", 50), rep("Oxford Street", 50))

# One row per participant, site and time point; the change is zero before the walk
expo.300 <- rbind(expo.delta, expo.delta, expo.delta)
expo.300$Per <- c(rep("before", 100), rep("after_1", 100), rep("after_2", 100))
expo.300[expo.300$Per == "before", c("CBLK_delta", "NO2_delta", "PM10_delta", "PM25_delta", "PCNT_delta")] <- 0
expo.300 <- expo.300[order(expo.300$Site), ]
expo.300 <- expo.300[order(expo.300$id), ]
expo.300 <- merge(expo.300, expo.annual, all.x = TRUE, by = "id")

expo.measured.100 <- expo.measured.hp %>%
  mutate(Site = "Hyde Park") %>%
  bind_rows(expo.measured.ox %>% mutate(Site = "Oxford Street"))


# Covariates and exposures: metabolomics (300 samples) --------------------
rownames(Full_covar_add) <- Full_covar_add$XP_ID
rownames(mRNA.full) <- mRNA.full$XP_ID
metab.tech.confounders.300 <- metab.tech.confounders[rownames(Full_covar_add), ]

covars.300 <- merge(Full_covar_add[, c("ID", "XP_ID", "Group", "Loc", "Per", "Age", "sex", "BMI", "Caffeine", "NOIS", "RHUM", "temp")],
  metab.tech.confounders.300[, c("XP_ID", "Site", "Time_Point", "Box_No", "Plate_No")],
  all.x = TRUE, by = "XP_ID"
)
rownames(covars.300) <- covars.300$XP_ID
covars.300 <- covars.300[as.character(rownames(Full_covar_add)), ]

# wave: site x time point combination (one residual variance per wave in the MVN models)
covars.300$wave <- factor(as.numeric(factor(10 * as.numeric(factor(covars.300$Site)) +
  as.numeric(factor(covars.300$Time_Point)))))
covars.300$Age <- as.numeric(covars.300$Age)
covars.300$Caffeine <- as.factor(covars.300$Caffeine)
covars.300$Group <- relevel(as.factor(covars.300$Group), "HEALTHY")

covars.expo.300 <- merge(covars.300, expo.300, by.x = c("ID", "Site", "Per"), by.y = c("id", "Site", "Per"), all.x = TRUE)
rownames(covars.expo.300) <- covars.expo.300$XP_ID
covars.expo.300 <- covars.expo.300[rownames(covars.300), ]


# Covariates and exposures: transcriptomics (252 samples) -----------------
# Samples with covariates, exposures and transcripts
select.idx.252 <- intersect(rownames(Full_covar_add), rownames(mRNA.full))

mRNA.tech.290 <- mRNA.full[, c(1:12, 14)]
rownames(mRNA.tech.290) <- mRNA.tech.290$XP_ID
colnames(mRNA.tech.290)[1] <- "ID"
mRNA.290 <- mRNA.full[, 15:ncol(mRNA.full)]

covars.252 <- Full_covar_add[select.idx.252, c("ID", "XP_ID", "Group", "Loc", "Per", "Age", "sex", "BMI", "Caffeine", "NOIS", "RHUM", "temp")]
covars.252 <- merge(covars.252, metab.tech.confounders[, c("XP_ID", "Site", "Time_Point")], all.x = TRUE, by = "XP_ID")
covars.252 <- merge(covars.252, mRNA.tech[, c("XP_ID", "array", "batch")], all.x = TRUE, by = "XP_ID")
rownames(covars.252) <- covars.252$XP_ID
covars.252$wave <- factor(as.numeric(factor(10 * as.numeric(factor(covars.252$Site)) +
  as.numeric(factor(covars.252$Time_Point)))))
covars.252$Age <- as.numeric(covars.252$Age)
covars.252$Caffeine <- as.factor(covars.252$Caffeine)
covars.252$batch <- as.integer(covars.252$batch)
covars.252$Group <- relevel(as.factor(covars.252$Group), "HEALTHY")

covars.expo.252 <- merge(covars.252, expo.300, by.x = c("ID", "Site", "Per"), by.y = c("id", "Site", "Per"), all.x = TRUE)
rownames(covars.expo.252) <- covars.expo.252$XP_ID
mRNA.252 <- mRNA.290[as.character(rownames(covars.expo.252)), ]
mRNA.tech.252 <- mRNA.tech.290[as.character(rownames(covars.expo.252)), ]


# Metabolomics ------------------------------------------------------------
metabolites.300 <- metabolites.full[rownames(Full_covar_add), ]

# Without the five outlying samples
XP.ID.no.outliers <- as.character(intersect(rownames(metabolites.no.outliers), rownames(Full_covar_add)))
covars.expo.295 <- covars.expo.300[XP.ID.no.outliers, ]
metabolites.295 <- metabolites.no.outliers[XP.ID.no.outliers, ]
metab.tech.confounders.295 <- metab.tech.confounders.300[XP.ID.no.outliers, ]


# Save data ---------------------------------------------------------------
saveRDS(covars.expo.300, "New_data/covars_expo_300.rds")
saveRDS(metabolites.300, "New_data/metabolites_300.rds")
saveRDS(metab.tech.confounders.300, "New_data/metab_tech_confounders_300.rds")

saveRDS(covars.expo.295, "New_data/covars_expo_295.rds")
saveRDS(metabolites.295, "New_data/metabolites_295.rds")
saveRDS(metab.tech.confounders.295, "New_data/metab_tech_confounders_295.rds")

saveRDS(covars.expo.252, "New_data/covars_expo_252.rds")
saveRDS(mRNA.tech.290, "New_data/mRNA_tech_290.rds")
saveRDS(mRNA.290, "New_data/mRNA_290.rds")
saveRDS(mRNA.tech.252, "New_data/mRNA_tech_252.rds")
saveRDS(mRNA.252, "New_data/mRNA_252.rds")

saveRDS(expo.measured.100, "New_data/expo.measured.100.rds")


# Descriptive statistics --------------------------------------------------
# Mean exposure at each site
expo.measured.100 %>%
  group_by(Site) %>%
  summarise(across(CBLK:PCNT, function(x) mean(x, na.rm = TRUE)))

# Number of participants and of samples per health group with transcriptomics data
id_group <- read_excel("Data/All_participants_ID_group.xlsx") %>%
  rename("ID" = `Subject ID`) %>%
  mutate(ID = as.character(as.integer(str_remove(ID, "RS"))))
subject_of <- function(xp_id) as.character(sapply(xp_id, function(x) substr(x, 1, nchar(x) - 2)))

count_by_group <- function(xp_id) {
  participants <- id_group %>%
    filter(ID %in% unique(subject_of(xp_id))) %>%
    count(Group, name = "participants")
  samples <- data.frame(ID = subject_of(xp_id)) %>%
    left_join(id_group, by = "ID") %>%
    count(Group, name = "samples")
  full_join(participants, samples, by = "Group")
}
count_by_group(rownames(mRNA.252)) # analysed
count_by_group(rownames(mRNA.290)) # before filtering on covariates
count_by_group(mRNA$XP_ID) # before removing failed arrays
