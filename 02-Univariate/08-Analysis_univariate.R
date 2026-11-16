# =============================================================================
# 08-Analysis_univariate.R
#
# Summarises the association studies: selects the metabolic features
# (Bonferroni) and transcripts (effective number of tests) associated with at
# least one TRAP exposure, annotates them, and draws the volcano, Manhattan
# and UpSet plots. Also writes the input files for the Mummichog pathway
# analysis, plots selected tryptophan metabolites, and computes the
# proportion of missing (imputed) values among the selected features.
#
# Input:  Results/MVN/{metab,mRNA}/<exposure>_*_mvn.rds (06-, 07-MVN_*.R)
#         Results/metab_annotation_1309.csv (feature annotation, see README)
#         Results/mRNA_annotation/annotation.matched.rds (04-mRNA_annotation.R)
#         New_data/*.rds
# Output: Results/MVN/sign.*.rds, Results/MVN/*_mvn_pval_beta_list.rds
#         Results/sign.metab.info.rds, Results/sign.metab.info.beta.rds
#         Results/MVN/sign_mRNA_table_with_genename.csv
#         Results/MVN/sign_metab_annotation_pct_missing.csv
#         Results/Mummichog/<exposure>.res.txt (input for MetaboAnalyst)
#         Results/Figures/*.pdf
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ComplexUpset)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(nlme)

source("Scripts/functions.R")

dir.create("Results/Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Mummichog", showWarnings = FALSE, recursive = TRUE)

mRNA <- readRDS("New_data/mRNA_denoised_252.rds")
metabolomics <- readRDS("New_data/metabolites_denoised_295.rds")
covars_expo <- readRDS("New_data/covars_expo_300.rds")
metab.list <- colnames(metabolomics)
mRNA.list <- colnames(mRNA)


# MVN results -------------------------------------------------------------
# Exposures as named in the paper (names) and in the result files (values)
expo.files <- c(PM25 = "PM25", PM10 = "PM10", NO2 = "NO2", BC = "CBLK", PCNT = "PCNT")
read_mvn <- function(omic) {
  lapply(expo.files, function(e) readRDS(sprintf("Results/MVN/%s/%s_%s_mvn.rds", omic, e, omic)))
}
metab.mvn <- read_mvn("metab")
mRNA.mvn <- read_mvn("mRNA")

# Estimates for the exposure: one row per feature, one column per exposure
expo_table <- function(mvn, stat) {
  sapply(mvn, function(x) x[[stat]][, "Expo.delta"])
}
metab.pval.df <- data.frame(metab = metab.list, expo_table(metab.mvn, "pval"))
metab.beta.df <- data.frame(metab = metab.list, expo_table(metab.mvn, "beta"))
mRNA.pval.df <- data.frame(mRNA = mRNA.list, expo_table(mRNA.mvn, "pval"))
mRNA.beta.df <- data.frame(mRNA = mRNA.list, expo_table(mRNA.mvn, "beta"))

saveRDS(list(pval = metab.pval.df, beta = metab.beta.df), "Results/MVN/metab_mvn_pval_beta_list.rds")
saveRDS(list(pval = mRNA.pval.df, beta = mRNA.beta.df), "Results/MVN/mRNA_mvn_pval_beta_list.rds")


# Selected features -------------------------------------------------------
# Metabolic features: Bonferroni
sign.metab.list <- lapply(c(NO2 = "NO2", BC = "BC", PCNT = "PCNT", PM10 = "PM10", PM25 = "PM25"), function(e) {
  metab.list[metab.pval.df[[e]] < alpha_metab]
})
sign.metab.unique <- unique(unlist(sign.metab.list, use.names = FALSE))

# Metabolic features significant at the ENT threshold, used to colour the
# network nodes by exposure (05-Networks)
sign.metab.ent.list <- lapply(c(NO2 = "NO2", CBLK = "BC", PCNT = "PCNT", PM10 = "PM10", PM25 = "PM25"), function(e) {
  metab.list[metab.pval.df[[e]] < alpha_metab_ent]
})

# Transcripts: ENT
sign.mRNA.list <- lapply(c(NO2 = "NO2", PM25 = "PM25", PM10 = "PM10", BC = "BC", PCNT = "PCNT"), function(e) {
  mRNA.list[mRNA.pval.df[[e]] < alpha_mrna]
})
sign.mRNA.unique <- unique(unlist(sign.mRNA.list, use.names = FALSE))

saveRDS(sign.metab.list, "Results/MVN/sign.metab.list.rds")
saveRDS(sign.metab.unique, "Results/MVN/sign.metab.unique.rds")
saveRDS(sign.metab.ent.list, "Results/MVN/sign.metab.ent.list.rds")
saveRDS(sign.mRNA.unique, "Results/MVN/sign.mRNA.unique.rds")


# Annotation of the metabolic features ------------------------------------
sign.metab.info <- metab.pval.df %>%
  mutate(
    across(PM25:PCNT, ~ .x < alpha_metab),
    `Number of associated TRAP` = PM25 + PM10 + NO2 + BC + PCNT
  ) %>%
  filter(`Number of associated TRAP` > 0) %>%
  dplyr::select(-c(PM25:PCNT)) %>%
  mutate(
    `Retention time` = as.numeric(str_split_fixed(metab, "@", 2)[, 2]),
    `Mass` = as.numeric(str_split_fixed(metab, "@", 2)[, 1])
  ) %>%
  relocate(`Number of associated TRAP`, .after = last_col()) %>%
  arrange(desc(`Number of associated TRAP`), `Retention time`)

metab.annotation <- read.csv("Results/metab_annotation_1309.csv", fileEncoding = "UTF-8-BOM")
colnames(metab.annotation)[5] <- "Annotation"

sign.metab.info <- sign.metab.info %>%
  left_join(., metab.annotation %>% dplyr::select(metab, Annotation), by = "metab") %>%
  mutate(Annotation = str_remove(Annotation, "\\*"))
sign.metab.info.beta <- sign.metab.info %>%
  left_join(., metab.beta.df, by = "metab")

saveRDS(sign.metab.info, "Results/sign.metab.info.rds")
saveRDS(sign.metab.info.beta, "Results/sign.metab.info.beta.rds")


# Annotation of the transcripts -------------------------------------------
z.table.unique <- readRDS("Results/mRNA_annotation/annotation.matched.rds")
z.table.sign <- z.table.unique[z.table.unique$NAME %in% sign.mRNA.unique, c("NAME", "GENE_SYMBOL", "GENE_NAME", "GO_ID")]
rm(z.table.unique)
saveRDS(z.table.sign, "Results/mRNA_annotation/annotation.significant.rds")

# Transcripts mapped to a gene, with the exposure(s) they are associated with
sign.mRNA.table <- z.table.sign %>%
  filter(GENE_NAME != "") %>%
  dplyr::select(NAME, GENE_SYMBOL) %>%
  mutate(
    NO2 = ifelse(NAME %in% sign.mRNA.list$NO2, "NO2", NA),
    BC = ifelse(NAME %in% sign.mRNA.list$BC, "BC", NA),
    PM25 = ifelse(NAME %in% sign.mRNA.list$PM25, "PM25", NA),
    PM10 = ifelse(NAME %in% sign.mRNA.list$PM10, "PM10", NA),
    PCNT = ifelse(NAME %in% sign.mRNA.list$PCNT, "PCNT", NA)
  ) %>%
  rename(
    "Accession number" = NAME,
    "Gene symbol" = GENE_SYMBOL
  ) %>%
  unite(., col = "Associated exposure", NO2:PCNT, sep = "+", na.rm = TRUE) %>%
  left_join(., mRNA.beta.df, by = c("Accession number" = "mRNA")) %>%
  arrange(`Associated exposure`)
write.csv(sign.mRNA.table, "Results/MVN/sign_mRNA_table_with_genename.csv", row.names = FALSE)


# Volcano plots -----------------------------------------------------------
metab.pval.beta.df <- bind_rows(
  metab.pval.df %>% pivot_longer(PM25:PCNT, names_to = "expo") %>% mutate(type = "pval"),
  metab.beta.df %>% pivot_longer(PM25:PCNT, names_to = "expo") %>% mutate(type = "beta")
)
mRNA.pval.beta.df <- bind_rows(
  mRNA.pval.df %>% pivot_longer(PM25:PCNT, names_to = "expo") %>% mutate(type = "pval"),
  mRNA.beta.df %>% pivot_longer(PM25:PCNT, names_to = "expo") %>% mutate(type = "beta")
)

volcano <- function(df, id, alpha, label) {
  df %>%
    pivot_wider(names_from = type, values_from = value, id_cols = all_of(c(id, "expo"))) %>%
    mutate(Significance = case_when(pval < alpha ~ label, TRUE ~ "Insignificant")) %>%
    ggplot(aes(beta, -log10(pval), col = Significance, alpha = Significance)) +
    geom_point(size = 0.5) +
    facet_wrap(~expo, scale = "free") +
    theme_bw() +
    xlab("Beta coefficient for change in TRAP levels") +
    ylab(expression(paste("-log"[10] * "(", italic("p"), "-value)"))) +
    scale_colour_manual(values = c("red", "darkgrey")) +
    scale_alpha_manual(values = c(1, 0.4)) +
    theme(legend.position = "top") +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed")
}
ggsave("Results/Figures/Volcano_metabolomics_col_significance.pdf",
  volcano(metab.pval.beta.df, "metab", alpha_metab, "Bonferroni significant")
)
ggsave("Results/Figures/Volcano_mRNA_col_significance.pdf",
  volcano(mRNA.pval.beta.df, "mRNA", alpha_mrna, "ENT significant")
)


# Manhattan plots ---------------------------------------------------------
# Mass and retention time are encoded in the feature ID: <mass>@<retention time>
add_mass_rt <- function(df) {
  df %>% mutate(
    mass = as.numeric(str_split_fixed(str_split_fixed(metab, ":", 2)[, 1], "@", 2)[, 1]),
    rt = as.numeric(str_split_fixed(str_split_fixed(metab, ":", 2)[, 1], "@", 2)[, 2])
  )
}

p <- metab.pval.beta.df %>%
  filter(type == "pval") %>%
  add_mass_rt() %>%
  mutate(sign = case_when(value < alpha_metab ~ "Bonferroni significant", TRUE ~ "Insignificant")) %>%
  ggplot(aes(rt, -log10(value), col = sign)) +
  geom_point(size = 0.5) +
  facet_grid(rows = vars(expo)) +
  theme_bw() +
  scale_color_manual(values = c("red", "grey90"), name = "Significance") +
  theme(legend.position = "top") +
  xlab("Retention time (min)") +
  ylab(expression(paste("-log"[10] * "(", italic("p"), "-value)"))) +
  scale_x_continuous(breaks = seq(0, 11, 1)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))
ggsave("Results/Figures/Manhattan_metabolomics_bonf.pdf", p, height = 6, width = 11)

# NO2: features selected with NO2, labelled with their annotation
fix_annotation <- function(x) {
  case_when(
    x == "Stearoylcarnitine (C18:0)?" ~ "Stearoylcarnitine (C18:0)",
    x == "Undecanoylcarnitine (C11:0)?" ~ "Undecanoylcarnitine (C11:0)",
    x == "phenylacetylcarnitine?" ~ "Phenylacetylcarnitine",
    TRUE ~ x
  )
}
metab.pval.beta.df.annotation <- metab.pval.beta.df %>%
  filter(expo == "NO2") %>%
  pivot_wider(id_cols = c(expo, metab), names_from = type, values_from = value) %>%
  add_mass_rt() %>%
  filter(pval < alpha_metab) %>%
  left_join(., sign.metab.info %>% dplyr::select(metab, Annotation), by = "metab") %>%
  mutate(Annotation = fix_annotation(Annotation))

# Version 1: compounds highlighted in colour, others labelled
highlighted <- c("Caffeine", "Phenylalanine", "Cortisol", "Indolelactic acid", "Isoleucine", "Indole-3-propionic acid", "Paraxanthine", "Unknown")
p <- metab.pval.beta.df.annotation %>%
  mutate(
    feature_color = factor(ifelse(Annotation %in% highlighted, Annotation, "Other"), levels = c(highlighted, "Other")),
    Annotation = ifelse(Annotation %in% highlighted, "", Annotation)
  ) %>%
  ggplot(aes(rt, -log10(pval), col = feature_color)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_color_manual(values = c(brewer.pal(8, "Set1"), "grey80"), name = "") +
  xlab("Retention time (min)") +
  ylab(expression(paste("-log"[10] * "(", italic("p"), "-value)"))) +
  scale_x_continuous(breaks = seq(0, 11, 1)) +
  geom_label_repel(aes(label = Annotation),
    color = "black", size = 3.5, box.padding = 0.5, point.padding = NA,
    max.overlaps = 13
  )
ggsave("Results/Figures/Manhattan_metabolomics_bonf_annotated_NO2.pdf", p, height = 6, width = 11)

# Version 2: one letter per annotated compound
mapping <- setNames(
  c("", letters[1:length(unique(metab.pval.beta.df.annotation$Annotation)) - 1]),
  unique(metab.pval.beta.df.annotation$Annotation)
)
colours.16 <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#BEAED4", "#FDBF6F",
  "#FF7F00", "#7FC97F", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FDC086", "grey50"
)
p <- metab.pval.beta.df.annotation %>%
  mutate(
    label = mapping[Annotation],
    label_annotation = case_when(Annotation == "Unknown" ~ Annotation, TRUE ~ paste0(label, " - ", Annotation)),
    Sign = factor(ifelse(beta < 0, "Negative", "Positive"), levels = c("Positive", "Negative"))
  ) %>%
  ggplot(aes(rt, -log10(pval))) +
  geom_point(aes(col = label_annotation)) +
  geom_text_repel(aes(label = label, fontface = 1), size = 3) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank()) +
  xlab("Retention time (min)") +
  ylab(expression(paste("-log"[10] * "(", italic("p"), "-value)"))) +
  scale_colour_manual(values = colours.16) +
  scale_x_continuous(breaks = seq(0, 11, 1))
ggsave("Results/Figures/Manhattan_metabolomics_NO2_annotation.pdf", p)
ggsave("Results/Figures/Manhattan_metabolomics_NO2_annotation_posneg.pdf", p + facet_grid(rows = vars(Sign)))


# UpSet plots -------------------------------------------------------------
expo_order <- c("PM25", "PM10", "PCNT", "BC", "NO2")

metab.sign.tf.df <- metab.pval.df %>%
  filter(metab %in% sign.metab.unique) %>%
  column_to_rownames("metab") %>%
  mutate(across(everything(), ~ .x < alpha_metab))
mRNA.sign.tf.df <- mRNA.pval.df %>%
  filter(mRNA %in% sign.mRNA.unique) %>%
  column_to_rownames("mRNA") %>%
  mutate(across(everything(), ~ .x < alpha_mrna))

upset_plot <- function(df, ymax, ...) {
  upset(df, expo_order,
    base_annotations = list("Intersection size" = intersection_size(text = list(size = 3))),
    set_sizes = (upset_set_size() +
      geom_text(aes(label = after_stat(count)), hjust = 1.1, stat = "count", col = "black", size = 4) +
      expand_limits(y = ymax)),
    themes = upset_default_themes(text = element_text(size = 15)),
    name = "", ...
  )
}
metab.upset <- upset_plot(metab.sign.tf.df, ymax = 85)
mrna.upset <- upset_plot(mRNA.sign.tf.df[, expo_order], ymax = 37, sort_sets = FALSE)
p <- (metab.upset | mrna.upset) + plot_layout(guides = "collect")
ggsave("Results/Figures/upset_combined.pdf", p, width = 12, height = 6)


# Input for Mummichog -----------------------------------------------------
# m/z of [M+H]+ ions, p-values and t-scores of all features, for the
# "functional analysis" module of MetaboAnalyst (positive mode)
mummichog_input <- function(mvn) {
  data.frame(
    m.z = as.numeric(str_split_fixed(metab.list, "@", 2)[, 1]) + 1.00728,
    `p-value` = mvn$pval[, "Expo.delta"],
    `t-score` = mvn$tval[, "Expo.delta"],
    check.names = FALSE
  )
}
for (e in names(metab.mvn)) {
  write.table(mummichog_input(metab.mvn[[e]]), paste0("Results/Mummichog/", e, ".res.txt"), quote = FALSE, row.names = FALSE)
}


# Exposure levels over time -----------------------------------------------
covars_expo %>%
  pivot_longer(CBLK_delta:PCNT_delta, names_to = "Expo_delta") %>%
  ggplot(aes(Time_Point, value)) +
  geom_point() +
  geom_line(aes(group = ID)) +
  facet_grid(cols = vars(Site), rows = vars(Expo_delta), scale = "free") +
  theme_bw()


# Tryptophan metabolites --------------------------------------------------
# Feature levels by site and time point, one line per participant
spaghetti <- function(features, title) {
  metabolomics %>%
    as.data.frame() %>%
    dplyr::select(all_of(features)) %>%
    rownames_to_column(var = "XP_ID") %>%
    left_join(., covars_expo %>% mutate(XP_ID = as.character(XP_ID)) %>% dplyr::select(XP_ID, ID, Site, Time_Point), by = "XP_ID") %>%
    mutate(Time_Point = case_when(
      Time_Point == "T1 (8AM)" ~ "Time point 1",
      Time_Point == "T2 (4PM)" ~ "Time point 2",
      TRUE ~ "Time point 3"
    )) %>%
    pivot_longer(cols = all_of(features), names_to = "metab") %>%
    mutate(Site_Time = paste(Site, Time_Point, sep = "\n")) %>%
    ggplot(aes(Site_Time, value)) +
    geom_point() +
    geom_line(aes(group = ID)) +
    facet_wrap(~metab) +
    theme_bw() +
    ggtitle(title) +
    xlab("") +
    ylab("log2(Peak intensity)")
}

# Indolelactic acid and indolepropionic acid
ia.p1 <- spaghetti(c("227.0561@3.8431869", "205.0735@3.8432918"), "Indolelactic acid")
ia.p2 <- spaghetti(c("211.0606@4.565889", "189.0787@4.5659895"), "Indolepropionic acid") +
  scale_y_continuous(breaks = seq(10, 20, 2))
ggsave("Results/Figures/Indolelactate_indoproprionate_spaghetti.pdf", ia.p1 / ia.p2, height = 6, width = 12)

# Kynurenine and tryptophan: features within 7 ppm of their masses
k_t <- c(208.0847, 204.0898)
metab.mass.6040 <- as.numeric(str_split_fixed(colnames(metabolomics), "@", 2)[, 1])
for (mass.a in k_t) {
  delta <- mass.a * 7e-06
  matched <- which(mass.a > metab.mass.6040 - delta & mass.a < metab.mass.6040 + delta)
  print(paste(mass.a, metab.mass.6040[matched], matched, "Match"))
}
metab.pval.beta.df %>%
  filter(metab %in% colnames(metabolomics)[c(3872, 4181, 4991)], expo == "NO2") %>%
  pivot_wider(id_cols = c("metab", "expo"), values_from = value, names_from = type)

ggsave("Results/Figures/Kynurenine_spaghetti.pdf", spaghetti(c("208.0842@1.8398161", "208.0847@1.8400313"), "Kynurenine"), height = 3, width = 12)
ggsave("Results/Figures/Tryptophan_spaghetti.pdf", spaghetti("204.0902@2.490191", "Tryptophan"), height = 3, width = 6)

# Kynurenine/tryptophan ratio and NO2
metab.ratios <- metabolomics %>%
  as.data.frame() %>%
  dplyr::select(c("208.0842@1.8398161", "208.0847@1.8400313", "204.0902@2.490191")) %>%
  rownames_to_column(var = "XP_ID") %>%
  left_join(., covars_expo %>% mutate(XP_ID = as.character(XP_ID)) %>%
    dplyr::select(XP_ID, ID, Site, Time_Point, NO2_delta, Group, Age, sex, BMI, wave), by = "XP_ID") %>%
  mutate(
    ratio1 = `208.0842@1.8398161` / `204.0902@2.490191`,
    ratio2 = `208.0847@1.8400313` / `204.0902@2.490191`
  )
for (ratio in c("ratio1", "ratio2")) {
  model.ratio <- gls(as.formula(paste(ratio, "~ NO2_delta + Group + Age + sex + BMI")),
    data = metab.ratios[complete.cases(metab.ratios), ],
    correlation = corSymm(form = ~ 1 | ID),
    weights = varIdent(form = ~ 1 | wave),
    na.action = na.pass,
    control = list(singular.ok = TRUE)
  )
  print(summary(model.ratio)$tTable[, c("Value", "p-value")])
}


# Missing values among the selected features ------------------------------
metab.not.imputed <- readRDS("New_data/metabolites_logged_not_imputed.rds")
metab.not.imputed.295 <- metab.not.imputed[rownames(metab.not.imputed) %in% rownames(metabolomics), ]

missing.df <- data.frame(
  metab = sign.metab.unique,
  missing_obs = as.numeric(apply(metab.not.imputed.295[, sign.metab.unique], 2, function(x) sum(is.na(x))))
)
missing.df <- sign.metab.info %>%
  left_join(., missing.df, by = "metab") %>%
  mutate(missing_pct = round(100 * missing_obs / nrow(metabolomics), 2)) %>%
  arrange(desc(missing_pct))
write.csv(missing.df, "Results/MVN/sign_metab_annotation_pct_missing.csv")

hist(missing.df$missing_pct,
  breaks = 10, xlab = "Percentage missing",
  main = paste("Percentage missing amongst", length(sign.metab.unique), "significant metabolic features")
)
