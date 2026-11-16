# =============================================================================
# 11-Mummichog_pathways.R
#
# Summarises the Mummichog pathway analysis and maps the selected metabolic
# features to the enriched pathways.
#
# Mummichog (v2.0) is run outside R, with the "functional analysis" module of
# MetaboAnalyst (https://www.metaboanalyst.ca), positive mode, 8 ppm mass
# tolerance and the default database. The input files are written by
# 08-Analysis_univariate.R (Results/Mummichog/<exposure>.res.txt). For each
# exposure, save the pathway and matched compound tables from MetaboAnalyst as
# Results/Mummichog/Output/Default/<exposure>_pathway.csv and
# Results/Mummichog/Output/Default/<exposure>_compound.csv.
#
# Input:  Results/Mummichog/Output/Default/*_{pathway,compound}.csv
#         Results/MVN/sign.metab.unique.rds (08-Analysis_univariate.R)
# Output: Results/Figures/Mummichog_heatmap_default.pdf
#         Results/Mummichog/metab.78.compounds.pathways.default.csv
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

dir.create("Results/Figures", showWarnings = FALSE, recursive = TRUE)

expo.files <- c(NO2 = "NO2", PCNT = "PCNT", BC = "BC", PM25 = "PM25", PM10 = "PM10")
read_mummichog <- function(e, table) {
  read.csv(sprintf("Results/Mummichog/Output/Default/%s_%s.csv", e, table))
}


# Enriched pathways -------------------------------------------------------
pathway_df <- bind_rows(lapply(names(expo.files), function(e) {
  read_mummichog(expo.files[[e]], "pathway") %>% mutate(TRAP = e)
})) %>%
  dplyr::rename("Pathway" = X, FET = `P.Fisher.`)

# Pathways enriched (Fisher exact test p < 0.05) for at least one exposure
sign_pathway_0.05 <- pathway_df %>%
  group_by(Pathway) %>%
  summarise(n_TRAP_sign = sum(FET <= 0.05)) %>%
  arrange(desc(n_TRAP_sign)) %>%
  filter(n_TRAP_sign > 0) %>%
  pull(Pathway)

sign_pathway_table_wide <- pathway_df %>%
  filter(Pathway %in% sign_pathway_0.05) %>%
  pivot_wider(names_from = TRAP, id_cols = c(Pathway, Pathway.total, Hits.total), values_from = FET) %>%
  replace(is.na(.), 1) %>%
  arrange(desc(Pathway.total))


# Heatmap -----------------------------------------------------------------
# Three blocks: pathway names, numbers of metabolites, and p-values by exposure
pathway.matrix <- sign_pathway_table_wide$Pathway %>% as.matrix()
pathway.col.matrix <- rep(c(0, 1), length.out = nrow(pathway.matrix)) # alternating row shading
pathway.col.df <- as.data.frame(pathway.col.matrix)
colnames(pathway.col.df) <- "Molecular Pathways"

metab.no.matrix <- sign_pathway_table_wide %>%
  dplyr::select(Pathway.total, Hits.total) %>%
  as.matrix()
metab.no.col.matrix <- cbind(pathway.col.matrix, pathway.col.matrix)
colnames(metab.no.col.matrix) <- c("Number of\nmetabolites\nin pathway", "Number of\noverlapping\nfeatures")

pval.matrix <- sign_pathway_table_wide %>%
  dplyr::select(NO2:PM10) %>%
  as.matrix()
pval.matrix.rounded <- round(pval.matrix, digits = 2)
pval.col.matrix <- as.data.frame(pval.matrix) %>%
  mutate(across(everything(), ~ case_when(
    .x < 0.05 ~ "< 0.05",
    .x < 0.1 ~ "< 0.1",
    TRUE ~ ">= 0.1"
  ))) %>%
  as.matrix()

color_mapping_table <- colorRamp2(c(0, 1), c("grey95", "white"))
text_cells <- function(values) {
  function(j, i, x, y, width, height, fill) {
    grid.text(values[i, j], x = x, y = y, gp = gpar(fontsize = 8))
  }
}

ht_opt(DIMNAME_PADDING = unit(5, unit = "mm"))
heatmap_columns <- list(
  height = unit(9, "cm"),
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 8)
)
p1 <- do.call(Heatmap, c(list(pathway.col.df,
  col = color_mapping_table,
  cluster_rows = FALSE, show_heatmap_legend = FALSE,
  width = unit(7, "cm"),
  cell_fun = text_cells(pathway.matrix)
), heatmap_columns))
p2 <- do.call(Heatmap, c(list(metab.no.col.matrix,
  col = color_mapping_table,
  cluster_rows = FALSE, cluster_columns = FALSE, show_heatmap_legend = FALSE,
  width = unit(4, "cm"),
  cell_fun = text_cells(metab.no.matrix)
), heatmap_columns))
p3 <- do.call(Heatmap, c(list(pval.col.matrix,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = structure(c("orange", "skyblue", "white"), names = c("< 0.05", "< 0.1", ">= 0.1")),
  width = unit(6, "cm"),
  name = "p-value",
  cell_fun = text_cells(pval.matrix.rounded)
), heatmap_columns))

pdf("Results/Figures/Mummichog_heatmap_default.pdf", height = 4.2, width = 8)
draw(p1 + p2 + p3, ht_gap = unit(c(0, 0), "cm"), padding = unit(c(1, 0, 1, 0), unit = "cm"))
dev.off()


# Pathways of the selected metabolic features -----------------------------
metab.sign.list <- readRDS("Results/MVN/sign.metab.unique.rds")
metab.78.df <- data.frame(metab = metab.sign.list) %>%
  mutate(
    mass = as.numeric(str_split_fixed(metab, "@", 2)[, 1]),
    rt = as.numeric(str_split_fixed(metab, "@", 2)[, 2]),
    mz = mass + 1.00728 # [M+H]+
  )

# Compounds matched to each feature (within 8 ppm). All features were used as
# input for every exposure, so the matched compounds are the same for all
# exposures.
compounds <- read_mummichog("BC", "compound")
match_compounds <- function(features, compounds, ppm = 8e-06) {
  matched.list <- lapply(seq_len(nrow(features)), function(i) {
    mz <- features$mz[i]
    delta <- mz * ppm
    hits <- compounds[(compounds$Query.Mass < mz + delta) & (compounds$Query.Mass > mz - delta), ]
    if (nrow(hits) > 0) cbind(metab = features$metab[i], hits)
  })
  names(matched.list) <- features$metab
  return(matched.list)
}
matched.78.list <- match_compounds(metab.78.df, compounds)
matched.78.unique.df <- do.call(rbind.data.frame, matched.78.list) %>%
  distinct(metab, Matched.Compound, Matched.Form) %>%
  remove_rownames()

# Compounds in each pathway (across exposures)
pathways.all <- bind_rows(lapply(c("BC", "NO2", "PCNT", "PM10", "PM25"), function(e) {
  read_mummichog(e, "pathway") %>% mutate(Expo = e)
})) %>%
  rename(Pathway = X, FET = `P.Fisher.`, Cpd.Hits = cpd.hits)
pathways.FET.005 <- pathways.all %>%
  filter(FET < 0.05) %>%
  pull(Pathway) %>%
  unique()
pathways.FET.010 <- pathways.all %>%
  filter(FET < 0.1) %>%
  pull(Pathway) %>%
  unique()

pathways.hits <- pathways.all %>% dplyr::select(Pathway, Pathway.total, Cpd.Hits)
pathways.compound.list <- list()
for (i in 1:nrow(pathways.hits)) {
  p <- pathways.hits[i, "Pathway"]
  cpds <- str_split(pathways.hits[i, "Cpd.Hits"], ";")[[1]]
  pathways.compound.list[[p]] <- unique(c(pathways.compound.list[[p]], cpds))
}

# Pathways (among `pathways.to.use`) containing each compound
compound2pathway <- function(compounds, pathways.list, pathways.to.use = "everything") {
  if (!identical(pathways.to.use, "everything")) {
    pathways.list <- pathways.list[pathways.to.use]
  }
  compound.match.pathway.list <- list()
  for (cpd in compounds) {
    for (p in names(pathways.list)) {
      if (cpd %in% pathways.list[[p]]) {
        compound.match.pathway.list[[cpd]] <- unique(c(compound.match.pathway.list[[cpd]], p))
      }
    }
    if (is.null(compound.match.pathway.list[[cpd]])) {
      compound.match.pathway.list[[cpd]] <- NA
    }
  }
  return(compound.match.pathway.list)
}

# Pathways of each feature, through its matched compounds
metab2pathway <- function(metabs, compounds.pathways) {
  metab.pathways.list <- list()
  for (i in 1:nrow(metabs)) {
    m <- metabs[i, "metab"]
    cmp <- metabs[i, "Matched.Compound"]
    if (!all(is.na(compounds.pathways[[cmp]]))) {
      metab.pathways.list[[m]] <- unique(c(metab.pathways.list[[m]], compounds.pathways[[cmp]]))
    }
  }
  return(metab.pathways.list)
}

compounds.matched.005 <- compound2pathway(matched.78.unique.df$Matched.Compound, pathways.compound.list, pathways.FET.005)
compounds.matched.010 <- compound2pathway(matched.78.unique.df$Matched.Compound, pathways.compound.list, pathways.FET.010)
metab.pathways.FET.005.list <- metab2pathway(matched.78.unique.df, compounds.matched.005)
metab.pathways.FET.010.list <- metab2pathway(matched.78.unique.df, compounds.matched.010)

metab.pathway.df <- metab.78.df %>%
  dplyr::select(metab, mass, rt) %>%
  mutate(matched = metab %in% matched.78.unique.df$metab) %>%
  arrange(rt)
metab.pathway.df$pathway.005 <- as.character(sapply(metab.pathways.FET.005.list[metab.pathway.df$metab], FUN = paste, collapse = "; "))
metab.pathway.df$pathway.010 <- as.character(sapply(metab.pathways.FET.010.list[metab.pathway.df$metab], FUN = paste, collapse = "; "))
metab.pathway.df$Matched.compounds <- as.character(sapply(metab.pathway.df$metab, function(x) {
  paste(matched.78.list[[x]]$Matched.Compound, collapse = "; ")
}))
write.csv(metab.pathway.df, "Results/Mummichog/metab.78.compounds.pathways.default.csv")
