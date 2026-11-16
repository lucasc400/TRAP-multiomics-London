# =============================================================================
# 10-Clustering_description.R
#
# Describes the clusters of metabolic features: calibration of the consensus
# clustering, ranges in mass and retention time within clusters, and the
# annotation of the features in each cluster.
#
# Input:  Results/Network/metab_clustering_stab.rds (09-Clustering_metab.R)
#         Results/sign.metab.info.rds (08-Analysis_univariate.R)
# Output: Results/Figures/Clustering_calibration_plot.pdf
#         Results/Figures/Sign_metab_cluster_mass_rt.pdf
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(sharp)
library(stringr)
library(ggplot2)
library(dplyr)

dir.create("Results/Figures", showWarnings = FALSE, recursive = TRUE)

metab.cluster.stab <- readRDS("Results/Network/metab_clustering_stab.rds")


# Calibration -------------------------------------------------------------
pdf("Results/Figures/Clustering_calibration_plot.pdf", width = 8, height = 6)
par(mar = c(6, 6, 2, 2))
CalibrationPlot(metab.cluster.stab, xlab = "Number of clusters")
dev.off()


# Mass and retention time ranges within clusters --------------------------
myclusters <- Clusters(metab.cluster.stab)

# Clusters with more than one feature
multiv_clusters <- sort(unique(myclusters[duplicated(myclusters)]))
mymultiv_clusters <- myclusters[myclusters %in% multiv_clusters]

# Mass and retention time are encoded in the feature ID: <mass>@<retention time>
mass <- as.numeric(str_split_fixed(names(mymultiv_clusters), "@", 2)[, 1])
time <- as.numeric(str_split_fixed(str_split_fixed(names(mymultiv_clusters), ":", 2)[, 1], "@", 2)[, 2])
sapply(split(time, f = mymultiv_clusters), function(x) diff(range(x))) # retention time range (min)
sapply(split(mass, f = mymultiv_clusters), function(x) diff(range(x))) # mass range

cluster.df <- data.frame(
  metab = names(myclusters),
  cluster = as.numeric(myclusters)
) %>%
  mutate(
    mass = as.numeric(str_split_fixed(names(myclusters), "@", 2)[, 1]),
    rt = as.numeric(str_split_fixed(names(myclusters), "@", 2)[, 2])
  ) %>%
  group_by(cluster) %>%
  mutate(n_metab_in_cluster = n()) %>%
  ungroup()

p <- cluster.df %>%
  mutate(cluster_ns = case_when(
    n_metab_in_cluster == 1 ~ "Cluster with\none feature",
    TRUE ~ paste0(cluster, " (", n_metab_in_cluster, ")")
  )) %>%
  ggplot(aes(rt, mass, col = cluster_ns)) +
  geom_point() +
  geom_line(aes(group = cluster)) +
  theme_bw() +
  scale_colour_manual(
    values = c("#996633", "#CC3333", "#FFCC00", "#006633", "#CC9999", "#e87722", "#660066", "#000000", "#0019a8", "#0099CC", "#66CCCC", "grey90"),
    name = "Cluster ID\n(number of\nfeatures)"
  ) +
  xlab("Retention time (min)") +
  ylab("Mass")
ggsave("Results/Figures/Sign_metab_cluster_mass_rt.pdf", p)


# Annotation of the features in each cluster ------------------------------
sign.metab.info <- readRDS("Results/sign.metab.info.rds")
cluster.df %>%
  left_join(., sign.metab.info %>% dplyr::select(metab, `Number of associated TRAP`, Annotation), by = "metab") %>%
  arrange(cluster) %>%
  print(n = Inf)
