# =============================================================================
# 05-EDA.R
#
# Exploratory data analysis: PCA of the metabolomics data (identifying the
# five outlying samples) and of the transcriptomics data, caffeine intake by
# exposure, the effect of denoising on the PCA, and exposure levels at the
# two sites.
#
# Input:  New_data/metabolites_{300,295}.rds, New_data/metabolites_denoised_295.rds
#         New_data/mRNA_252.rds, New_data/mRNA_denoised_252.rds
#         New_data/covars_expo_300.rds, New_data/*tech*.rds
# Output: Results/Figures/PCA_*.pdf, Results/Figures/Expo_measured.pdf
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(ggplot2)
library(dplyr)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(stringr)
library(patchwork)
library(ggrepel)

dir.create("Results/Figures", showWarnings = FALSE, recursive = TRUE)

metab.300 <- readRDS("New_data/metabolites_300.rds")
covars_expo <- readRDS("New_data/covars_expo_300.rds")


# PCA ---------------------------------------------------------------------
# Metabolomics
metab.PCA <- PCA(metab.300, graph = FALSE)
p <- fviz_pca_ind(metab.PCA, geom.ind = "text", addEllipses = TRUE) + ggtitle("")
ggsave("Results/Figures/PCA_metabolomics.pdf", p, width = 6, height = 6)

# Label the five outlying samples only
outliers <- c("6011", "1511", "7623", "10123", "9721")
metab.PCA.plot <- fviz_pca_ind(metab.PCA, geom.ind = "point", addEllipses = TRUE) + ggtitle("")
obs.labels <- as.character(metab.PCA.plot$data$name)
p <- metab.PCA.plot + geom_text_repel(aes(label = ifelse(obs.labels %in% outliers, obs.labels, "")))
ggsave("Results/Figures/PCA_metabolomics_outliers_only.pdf", p, width = 6, height = 6)

# Transcriptomics
mRNA.252 <- readRDS("New_data/mRNA_252.rds")
mRNA.PCA <- PCA(mRNA.252, graph = FALSE)
p <- fviz_pca_ind(mRNA.PCA, geom.ind = "text", addEllipses = TRUE) + ggtitle("")
ggsave("Results/Figures/PCA_mRNA.pdf", p, width = 6, height = 6)

# Scree plot for transcriptomics
mRNA.ev <- as.numeric(mRNA.PCA$eig[, "cumulative percentage of variance"])
par(mar = c(4, 6, 4, 4))
plot(mRNA.ev,
  pch = 19, col = "navy", las = 1, type = "b", ylim = c(0, 100), cex = 0.5,
  ylab = "Cumulative explained\nvariance (%)", xlab = "PCs", cex.lab = 1.5
)


# Caffeine ----------------------------------------------------------------
covars_expo %>%
  dplyr::select(Caffeine, CBLK_delta:PCNT_delta, Site) %>%
  pivot_longer(CBLK_delta:PCNT_delta, names_to = "Exposure") %>%
  ggplot(aes(Site, value, fill = Caffeine)) +
  geom_boxplot() +
  facet_wrap(~Exposure, scales = "free") +
  theme_bw() +
  theme(legend.position = "top")

# One exposure value per participant and site
covars_expo %>%
  dplyr::select(Caffeine, CBLK_delta:PCNT_delta, Site, Per) %>%
  filter(Per == "after_1") %>%
  pivot_longer(CBLK_delta:PCNT_delta, names_to = "Exposure") %>%
  ggplot(aes(Site, value, fill = Caffeine)) +
  geom_boxplot() +
  facet_wrap(~Exposure, scales = "free") +
  theme_bw() +
  theme(legend.position = "top")

# Exposure by caffeine intake: t-tests
covars_expo %>%
  dplyr::select(Caffeine, CBLK_delta:PCNT_delta, Site, Per) %>%
  filter(Per == "after_1") %>%
  pivot_longer(CBLK_delta:PCNT_delta, names_to = "Exposure") %>%
  group_by(Site, Exposure) %>%
  summarise(pval = t.test(value[Caffeine == 1], value[Caffeine == 0])$p.value)

covars_expo %>%
  dplyr::select(Caffeine, CBLK_delta:PCNT_delta, Site) %>%
  pivot_longer(CBLK_delta:PCNT_delta, names_to = "Exposure") %>%
  group_by(Exposure, Site, Caffeine) %>%
  summarise(
    mean.expo = mean(value, na.rm = TRUE),
    sd.expo = sd(value, na.rm = TRUE)
  )


# Effect of denoising -----------------------------------------------------
# PCA of the data before and after denoising, coloured by technical covariate
pca_by <- function(pca, groups, legend_title, title, legend = TRUE) {
  p <- fviz_pca_ind(pca,
    geom.ind = "text",
    col.ind = groups,
    legend.title = legend_title,
    addEllipses = TRUE
  ) + ggtitle(title)
  if (!legend) {
    p <- p + theme(legend.position = "none")
  }
  p
}

# Metabolomics
metab_tech_confounders <- readRDS("New_data/metab_tech_confounders_295.rds")
metab_tech_confounders$Plate_No <- as.character(metab_tech_confounders$Plate_No)
metabolites.pca <- PCA(readRDS("New_data/metabolites_denoised_295.rds"), graph = FALSE)
metabolites.undenoised.pca <- PCA(readRDS("New_data/metabolites_295.rds"), graph = FALSE)

denoising.metab <- (pca_by(metabolites.undenoised.pca, metab_tech_confounders$Box_No, "Box", "Undenoised", legend = FALSE) +
  pca_by(metabolites.pca, metab_tech_confounders$Box_No, "Box", "Denoised")) /
  (pca_by(metabolites.undenoised.pca, metab_tech_confounders$Plate_No, "Plate", "Undenoised", legend = FALSE) +
    pca_by(metabolites.pca, metab_tech_confounders$Plate_No, "Plate", "Denoised")) +
  plot_annotation("Metabolomics")
ggsave("Results/Figures/PCA_metabolomics_denoising.pdf", denoising.metab, width = 10, height = 10)

# Transcriptomics
mRNA_tech_confounders <- readRDS("New_data/mRNA_tech_252.rds")
mRNA_tech_confounders$array <- as.character(mRNA_tech_confounders$array)
mRNA_tech_confounders$batch <- as.character(mRNA_tech_confounders$batch)
mRNA.pca <- PCA(readRDS("New_data/mRNA_denoised_252.rds"), graph = FALSE)
mRNA.undenoised.pca <- PCA(mRNA.252, graph = FALSE)

denoising.mRNA <- (pca_by(mRNA.undenoised.pca, mRNA_tech_confounders$array, "Array", "Undenoised", legend = FALSE) +
  pca_by(mRNA.pca, mRNA_tech_confounders$array, "Array", "Denoised")) /
  (pca_by(mRNA.undenoised.pca, mRNA_tech_confounders$batch, "Batch", "Undenoised", legend = FALSE) +
    pca_by(mRNA.pca, mRNA_tech_confounders$batch, "Batch", "Denoised")) +
  plot_annotation("mRNA")
ggsave("Results/Figures/PCA_mRNA_denoising.pdf", denoising.mRNA, width = 10, height = 10)


# Exposure levels ---------------------------------------------------------
# Measured levels = change from annual average + annual average
covars_expo <- covars_expo %>%
  mutate(
    BC_measured = CBLK_delta + CBLK_annual,
    NO2_measured = NO2_delta + NO2_annual,
    PM25_measured = PM25_delta + PM25_annual,
    PM10_measured = PM10_delta + PM10_annual,
    PCNT_measured = PCNT_delta + PCNT_annual
  )

expo.p1 <- covars_expo %>%
  pivot_longer(BC_measured:PM10_measured, names_to = "TRAP_measured") %>%
  mutate(TRAP_measured = str_remove(TRAP_measured, "_measured")) %>%
  ggplot(aes(TRAP_measured, value, fill = Site)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  theme(legend.position = "none") +
  ylab(expression(paste(mu, "g m"^-3))) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C"))

expo.p2 <- covars_expo %>%
  pivot_longer(PCNT_measured, names_to = "TRAP_measured") %>%
  mutate(TRAP_measured = str_remove(TRAP_measured, "_measured")) %>%
  ggplot(aes(TRAP_measured, value, fill = Site)) +
  geom_boxplot() +
  theme_bw() +
  xlab("") +
  ylab(expression(paste("Particle number count/ cm"^3))) +
  scale_fill_manual(values = c("#377EB8", "#E41A1C"))

p <- expo.p1 + expo.p2 + plot_layout(widths = c(4, 1), ncol = 2)
ggsave("Results/Figures/Expo_measured.pdf", p)
