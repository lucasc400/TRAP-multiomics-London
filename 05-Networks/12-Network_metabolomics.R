# =============================================================================
# 12-Network_metabolomics.R
#
# Conditional independence networks of the selected metabolic features
# (represented by their cluster medoids), estimated separately at the three
# time points (2 h before, 2 h and 24 h after the walk) with a graphical LASSO
# calibrated by stability (sharp). Nodes are coloured by enriched pathway and
# by associated exposure, and labelled with their annotation (SI Table S5).
#
# Input:  New_data/metabolites_denoised_295.rds, New_data/covars_expo_295.rds
#         Results/MVN/sign.metab.unique.rds, Results/MVN/sign.metab.ent.list.rds
#         Results/Clustering/metab_summarised.rds (09-Clustering_metab.R)
#         Results/Mummichog/metab.78.compounds.pathways.default.csv (11-Mummichog_pathways.R)
#         Supplementary_Materials/SI_Table_5.xlsx (annotation, SI Table S5 of the paper)
# Output: Results/Network/metab_63_tp{1,2,3}.rds
#         Results/Network/metab.63.pathway.default.compound.pdf
#         Results/Network/metab.63.exposure.compound.pdf
#         Results/Network/Calibration_metabolomics.pdf
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(igraph)
library(sharp)
library(dplyr)
library(readxl)
library(ggvenn)

dir.create("Results/Network", showWarnings = FALSE, recursive = TRUE)

covars_expo <- readRDS("New_data/covars_expo_295.rds")
metab.sign.list <- readRDS("Results/MVN/sign.metab.unique.rds")
metab_summarised <- readRDS("Results/Clustering/metab_summarised.rds")
n.nodes <- ncol(metab_summarised)


# Data at each time point -------------------------------------------------
time.points <- c(tp1 = "before", tp2 = "after_1", tp3 = "after_2")
time.titles <- c(tp1 = "2 h before", tp2 = "+2 h", tp3 = "+24 h")
metab.summarised.tp <- lapply(time.points, function(tp) metab_summarised[covars_expo$Per == tp, ])


# Pathways ----------------------------------------------------------------
metab_78_compounds_pathways <- read.csv("Results/Mummichog/metab.78.compounds.pathways.default.csv")
metab_78_compounds_pathways <- metab_78_compounds_pathways[, colnames(metab_78_compounds_pathways) != "X"]
metab_78_compounds_pathways$pathway.005 <- with(
  metab_78_compounds_pathways,
  ifelse(pathway.005 == "", "No matched pathways", pathway.005)
)

# Pathway names with line breaks for the legend
metab_78_compounds_pathways_clean <- metab_78_compounds_pathways %>%
  mutate(pathway.clean = case_when(
    pathway.005 == "Prostaglandin formation from arachidonate" ~ "Prostaglandin formation\nfrom arachidonate",
    pathway.005 == "Biopterin metabolism; Tryptophan metabolism" ~ "Biopterin metabolism;\nTryptophan metabolism",
    pathway.005 == "Tryptophan metabolism; Chondroitin sulfate degradation; Heparan sulfate degradation; Caffeine metabolism; Biopterin metabolism" ~
      "Tryptophan metabolism;\nChondroitin sulfate degradation;\nHeparan sulfate degradation;\nCaffeine metabolism;\nBiopterin metabolism",
    TRUE ~ pathway.005
  )) %>%
  pull(pathway.clean)
names(metab_78_compounds_pathways_clean) <- metab_78_compounds_pathways$metab

pathway.palette <- c("#56B4E9", "#E69F00", "#999999", "#009E73", "#CC3333", "#0072B2")
pathway.levels <- levels(as.factor(metab_78_compounds_pathways_clean))
colours.pathway.list <- pathway.palette[as.integer(as.factor(metab_78_compounds_pathways_clean))]
names(colours.pathway.list) <- metab_78_compounds_pathways$metab

# Number of features per pathway after summarising the clusters
metab_78_compounds_pathways %>%
  filter(metab %in% colnames(metab_summarised)) %>%
  count(pathway.005)


# Exposures ---------------------------------------------------------------
# Exposure(s) each feature is associated with (ENT threshold)
sign.metab.ent.list <- readRDS("Results/MVN/sign.metab.ent.list.rds")
metab.sign.by.expo <- sapply(metab.sign.list, function(m) {
  paste(names(sign.metab.ent.list)[sapply(sign.metab.ent.list, function(l) m %in% l)], collapse = "_")
})
expo.levels <- c("CBLK", "NO2", "PCNT", "PM10", "PM25", "CBLK_PCNT", "NO2_CBLK", "PM10_PM25", "CBLK_PM10")
expo.palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999", "#D55E00", "#CC79A7", "#CC3333")
colours.expo.list <- expo.palette[as.integer(factor(metab.sign.by.expo, levels = expo.levels))]
names(colours.expo.list) <- metab.sign.list
expo.in.legend <- c(1, 2, 3, 4, 6, 7, 8) # exposure groups present in the networks


# Annotation --------------------------------------------------------------
si.table.5 <- read_excel("Supplementary_Materials/SI_Table_5.xlsx")
si.table.5$`Metabolite name` <- ifelse(is.na(si.table.5$`Metabolite name`), "Unknown", si.table.5$`Metabolite name`)
metab.annotation.match <- structure(si.table.5$`Metabolite name`, names = si.table.5$`Metabolite ID`)


# Stability selection (seed 11 as in the published analysis) --------------
stab <- lapply(metab.summarised.tp, function(x) {
  GraphicalModel(x, PFER_thr = 50, K = 500, n_cat = 3, pi_list = seq(0.6, 0.9, by = 0.01), seed = 11)
})
for (tp in names(stab)) {
  saveRDS(stab[[tp]], paste0("Results/Network/metab_", n.nodes, "_", tp, ".rds"))
}

pdf("Results/Network/Calibration_metabolomics.pdf", width = 15, height = 4)
par(mar = c(7, 6, 6, 6), mfrow = c(1, 3))
for (tp in names(stab)) CalibrationPlot(stab[[tp]])
dev.off()


# Networks ----------------------------------------------------------------
# Node size and label size increase with degree; nodes numbered 1 to n.nodes,
# or labelled with their annotation when known
clean.graph <- function(grp) {
  grp <- set_vertex_attr(grp, "size", value = degree(grp) * 0.35 + 7)
  grp <- set_vertex_attr(grp, "label.cex", value = (degree(grp) * 0.6 + 6) / 13)
  grp <- set_vertex_attr(grp, "frame.color", value = rep("grey", vcount(grp)))
  grp <- set_vertex_attr(grp, "label", value = 1:vcount(grp))
  annotation <- metab.annotation.match[V(grp)$name]
  grp <- set_vertex_attr(grp, "label", value = ifelse(annotation == "Unknown", V(grp)$label, annotation))
  return(grp)
}
graphs <- lapply(stab, function(s) clean.graph(Graph(s, satellites = TRUE)))

# Only plot the features with at least one edge at one of the time points
network.seed <- 65
set.seed(network.seed)
layout <- layout_with_fr(graphs$tp1)
remov <- which(degree(graphs$tp1) == 0 & degree(graphs$tp2) == 0 & degree(graphs$tp3) == 0)
kept <- setdiff(seq_len(n.nodes), remov)

# Three networks side by side, with a legend on the right
plot_networks <- function(file, colours, legend_labels, legend_colours, legend_y, legend_spacing) {
  pdf(file, width = 12, height = 4.5)
  par(mfrow = c(1, 3), mar = c(1, 1, 3, 1), oma = c(0, 0, 0, 15), xpd = NA)
  for (tp in names(graphs)) {
    grp.sub <- induced_subgraph(graphs[[tp]], V(graphs[[tp]])[kept])
    plot(grp.sub, layout = layout[kept, ], vertex.color = colours[V(grp.sub)$name])
    title(time.titles[tp], line = -2)
    box()
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 15, 0, 0), mar = c(4, 15, 0, 4), new = TRUE)
  plot(0, 0, type = "l", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend(0.65, legend_y,
    legend = legend_labels, col = legend_colours,
    lwd = 5, xpd = TRUE, horiz = FALSE, cex = 1, seg.len = 1, bty = "n", y.intersp = legend_spacing
  )
  dev.off()
}

# Coloured by pathway (pathways of the plotted features only)
pathways.used <- data.frame(pathway = pathway.levels, col = pathway.palette[seq_along(pathway.levels)]) %>%
  filter(pathway %in% metab_78_compounds_pathways_clean[V(graphs$tp1)$name[kept]])
plot_networks(paste0("Results/Network/metab.", n.nodes, ".pathway.default.compound.pdf"),
  colours = colours.pathway.list,
  legend_labels = pathways.used$pathway, legend_colours = pathways.used$col,
  legend_y = 1, legend_spacing = 3
)

# Coloured by exposure
plot_networks(paste0("Results/Network/metab.", n.nodes, ".exposure.compound.pdf"),
  colours = colours.expo.list,
  legend_labels = expo.levels[expo.in.legend], legend_colours = expo.palette[expo.in.legend],
  legend_y = 0.9, legend_spacing = 4
)


# Edges by time point -----------------------------------------------------
edges <- lapply(graphs, function(g) {
  e <- as.data.frame(as_edgelist(g, names = TRUE))
  colnames(e) <- c("Start", "End")
  e$Edge <- paste(e$Start, e$End, sep = " - ")
  e
})
ggvenn(lapply(edges, function(e) e$Edge), digits = 0)

# Edges present only 2 h after the walk ("acute" edges), with annotations
acute <- edges$tp2[!edges$tp2$Edge %in% c(edges$tp1$Edge, edges$tp3$Edge), ]
metab.annotation.match.cluster <- ifelse(metab.annotation.match[colnames(metab_summarised)] != "Unknown",
  metab.annotation.match[colnames(metab_summarised)],
  1:n.nodes
)
names(metab.annotation.match.cluster) <- colnames(metab_summarised)
paste(metab.annotation.match.cluster[acute$Start], metab.annotation.match.cluster[acute$End], sep = " - ")
