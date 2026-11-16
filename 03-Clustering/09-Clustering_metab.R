# =============================================================================
# 09-Clustering_metab.R
#
# Groups the selected metabolic features into clusters of correlated
# features (typically ions of the same compound) using consensus clustering
# calibrated by stability (sharp), and represents each cluster by its medoid.
# The medoids are used in the conditional independence networks.
#
# Input:  New_data/metabolites_denoised_295.rds
#         Results/MVN/sign.metab.unique.rds, Results/sign.metab.info.beta.rds,
#         Results/MVN/sign_metab_annotation_pct_missing.csv (08-Analysis_univariate.R)
# Output: Results/Network/metab_clustering_stab.rds
#         Results/Network/Cluster_description.csv
#         Results/Clustering/metab_summarised.rds (medoids)
#         Results/sign.metab.info.beta.cluster.rds, Results/sign.metab.info.beta.cluster.missing.csv
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(sharp)
library(dplyr)

source("Scripts/functions.R")

dir.create("Results/Network", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Clustering", showWarnings = FALSE, recursive = TRUE)

metabolomics <- readRDS("New_data/metabolites_denoised_295.rds")
metab.sign.list <- readRDS("Results/MVN/sign.metab.unique.rds")
metab.sign <- metabolomics[, metab.sign.list]


# Consensus clustering ----------------------------------------------------
# Number of clusters calibrated on 500 subsamples of 80% of the observations
# (seed 11 as in the published analysis)
metab.cluster.stab <- Clustering(metab.sign, row = FALSE, nc = 1:ncol(metab.sign), tau = 0.8, K = 500, seed = 11)
saveRDS(metab.cluster.stab, "Results/Network/metab_clustering_stab.rds")


# Summarise each cluster by its medoid ------------------------------------
metab.clusters <- Clusters(metab.cluster.stab)
summarised <- GroupingSummary(x = metab.sign, group = metab.clusters, summary = "medoid")
saveRDS(summarised$xsummarised, "Results/Clustering/metab_summarised.rds")
write.csv(summarised$description, "Results/Network/Cluster_description.csv")


# Add clusters to the table of selected features --------------------------
metab.clusters.df <- summarised$description[, c("feature", "group")]
rownames(metab.clusters.df) <- NULL
colnames(metab.clusters.df) <- c("metab", "Cluster")

sign.metab.info.beta.cluster <- readRDS("Results/sign.metab.info.beta.rds") %>%
  left_join(., metab.clusters.df, by = "metab") %>%
  arrange(Cluster)
saveRDS(sign.metab.info.beta.cluster, "Results/sign.metab.info.beta.cluster.rds")

# With the proportion of missing values
missing.df <- read.csv("Results/MVN/sign_metab_annotation_pct_missing.csv", row.names = 1)
sign.metab.info.beta.cluster %>%
  left_join(., missing.df %>% dplyr::select(metab, missing_obs, missing_pct), by = "metab") %>%
  write.csv("Results/sign.metab.info.beta.cluster.missing.csv")
