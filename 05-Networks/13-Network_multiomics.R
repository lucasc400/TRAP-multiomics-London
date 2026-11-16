# =============================================================================
# 13-Network_multiomics.R
#
# Multi-omic conditional independence networks of the selected metabolic
# features (cluster medoids) and the selected transcripts mapped to a known
# gene, estimated at the three time points.
# The graphical LASSO is calibrated by stability with block-specific
# penalties (sharp). The main networks use a PFER threshold of 50;
# thresholds of 100 and 150 are used as a sensitivity analysis.
#
# Input:  New_data/covars_expo_295.rds, New_data/mRNA_denoised_252.rds
#         Results/mRNA_annotation/z.table.matched.rds (04-mRNA_annotation.R)
#         Results/MVN/sign.mRNA.unique.rds (08-Analysis_univariate.R)
#         Results/Clustering/metab_summarised.rds (09-Clustering_metab.R)
#         Supplementary_Materials/SI_Table_5.xlsx (annotation, SI Table S5 of the paper)
# Output: Results/Network/mRNA_38_metab_63_tp{1,2,3}[_100,_150].rds
#         Results/Network/metab.63.mRNA.38.annotation.compound[.pfer100,.pfer150].pdf
#         Results/Network/Summary_multiomic_network_by_PFER.pdf
#         Results/Network/Calibration_multiomic.pdf
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(igraph)
library(sharp)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

dir.create("Results/Network", showWarnings = FALSE, recursive = TRUE)

covars_expo <- readRDS("New_data/covars_expo_295.rds")
mRNA <- readRDS("New_data/mRNA_denoised_252.rds")
mRNA.annotation <- readRDS("Results/mRNA_annotation/z.table.matched.rds")
mRNA.sign.list <- readRDS("Results/MVN/sign.mRNA.unique.rds")
metab_summarised <- readRDS("Results/Clustering/metab_summarised.rds")


# Data --------------------------------------------------------------------
# Selected transcripts mapped to a gene, named by gene symbol
mRNA.w.gene <- mRNA.annotation[mRNA.annotation$GENE_NAME != "", "NAME"]
mRNA.sign.list.g <- intersect(mRNA.w.gene, mRNA.sign.list)

# Samples with both metabolomics and transcriptomics
common.id <- intersect(rownames(mRNA), rownames(metab_summarised))
covars_expo.common <- covars_expo[common.id, ]
mRNA.sign.g <- mRNA[common.id, mRNA.sign.list.g]
colnames(mRNA.sign.g) <- mRNA.annotation[mRNA.sign.list.g, "GENE_SYMBOL"]

n.metab <- ncol(metab_summarised)
n.mrna <- ncol(mRNA.sign.g)
file.prefix <- paste0("mRNA_", n.mrna, "_metab_", n.metab)
plot.prefix <- paste0("Results/Network/metab.", n.metab, ".mRNA.", n.mrna, ".annotation.compound")

# Data at each time point: metabolic features first, then transcripts (row
# selection as in the published analysis)
time.points <- c(tp1 = "before", tp2 = "after_1", tp3 = "after_2")
time.titles <- c(tp1 = "2 h before", tp2 = "+2 h", tp3 = "+24 h")
metab.mrna.tp <- lapply(time.points, function(tp) {
  in_tp <- which(covars_expo.common$Per == tp)
  cbind(metab_summarised[in_tp, ], mRNA.sign.g[in_tp, ])
})


# Annotation --------------------------------------------------------------
si.table.5 <- read_excel("Supplementary_Materials/SI_Table_5.xlsx")
si.table.5$`Metabolite name` <- ifelse(is.na(si.table.5$`Metabolite name`), "Unknown", si.table.5$`Metabolite name`)
metab.annotation.match <- structure(si.table.5$`Metabolite name`, names = si.table.5$`Metabolite ID`)


# Stability selection -----------------------------------------------------
pfer.thresholds <- c(50, 100, 150)
stab <- list()
for (pfer in pfer.thresholds) {
  stab[[as.character(pfer)]] <- lapply(metab.mrna.tp, function(x) {
    GraphicalModel(x, pk = c(n.metab, n.mrna), PFER_thr = pfer, K = 500, pi_list = seq(0.6, 0.9, by = 0.01), n_cat = 3, seed = 11)
  })
  for (tp in names(time.points)) {
    suffix <- ifelse(pfer == 50, "", paste0("_", pfer))
    saveRDS(stab[[as.character(pfer)]][[tp]], paste0("Results/Network/", file.prefix, "_", tp, suffix, ".rds"))
  }
}

# Calibration of the main networks (one panel per block)
pdf("Results/Network/Calibration_multiomic.pdf", width = 15, height = 10)
par(mar = c(7, 6, 7, 6), mfrow = c(3, 3), oma = c(0, 5, 5, 0))
for (tp in names(time.points)) CalibrationPlot(stab[["50"]][[tp]])
dev.off()


# Networks ----------------------------------------------------------------
# Metabolic features (blue) are numbered 1 to n.metab, or labelled with their
# annotation when known; transcripts (orange) are labelled with their gene
clean.graph <- function(grp) {
  is.metab <- seq_len(vcount(grp)) <= n.metab
  grp <- set_vertex_attr(grp, "size", value = degree(grp) * 0.35 + 7)
  grp <- set_vertex_attr(grp, "label.cex", value = rep(0.5, vcount(grp)))
  grp <- set_vertex_attr(grp, "color", value = ifelse(is.metab, "#56B4E9", "#E69F00"))
  grp <- set_vertex_attr(grp, "frame.color", value = rep("grey", vcount(grp)))
  annotation <- metab.annotation.match[V(grp)$name[is.metab]]
  metab.labels <- ifelse(annotation == "Unknown", seq_len(n.metab), annotation)
  grp <- set_vertex_attr(grp, "label", value = c(metab.labels, V(grp)$label[!is.metab]))
  return(grp)
}

# Three networks side by side, showing the nodes with at least one edge at
# one of the time points
plot_networks <- function(graphs, file, network.seed = 111) {
  remov <- which(degree(graphs$tp1) == 0 & degree(graphs$tp2) == 0 & degree(graphs$tp3) == 0)
  kept <- setdiff(seq_len(vcount(graphs$tp1)), remov)
  set.seed(network.seed)
  layout <- layout_with_fr(graphs$tp1)

  pdf(file, width = 12, height = 5)
  par(mfrow = c(1, 3), mar = c(1, 1, 3, 1), oma = c(0, 0, 0, 0), xpd = NA)
  for (tp in names(graphs)) {
    grp.sub <- induced_subgraph(graphs[[tp]], V(graphs[[tp]])[kept])
    plot(grp.sub, layout = layout[kept, ])
    title(time.titles[tp], line = -2)
    box()
    if (tp == "tp1") {
      legend(-1.1, 1.35,
        legend = c("metab", "mRNA"), col = c("#56B4E9", "#E69F00"),
        lwd = 5, xpd = TRUE, horiz = FALSE, seg.len = 1, cex = 1.4, bty = "n", y.intersp = 1
      )
    }
  }
  dev.off()
}

graphs <- lapply(stab, function(s) lapply(s, function(x) clean.graph(Graph(x, satellites = TRUE))))
plot_networks(graphs[["50"]], paste0(plot.prefix, ".pdf"))
for (pfer in c("100", "150")) {
  plot_networks(graphs[[pfer]], paste0(plot.prefix, ".pfer", pfer, ".pdf"))
}


# Number of edges by PFER threshold ---------------------------------------
# Blocks in the order used by sharp: within metabolomics, between the two
# omics, within transcriptomics
count_edges <- function(s) {
  A <- Adjacency(s)
  m <- seq_len(n.metab)
  r <- n.metab + seq_len(n.mrna)
  c(Metabolomics = sum(A[m, m]) / 2, `Cross-omics` = sum(A[m, r]), Transcriptomics = sum(A[r, r]) / 2)
}
network.summary <- bind_rows(lapply(names(stab), function(pfer) {
  bind_rows(lapply(names(time.points), function(tp) {
    data.frame(PFER = as.numeric(pfer), Timepoint = which(names(time.points) == tp), t(count_edges(stab[[pfer]][[tp]])), check.names = FALSE)
  }))
})) %>%
  pivot_longer(Metabolomics:Transcriptomics, names_to = "block", values_to = "value")
print(network.summary, n = Inf)

p <- network.summary %>%
  ggplot(aes(as.factor(Timepoint), value, col = as.factor(PFER))) +
  geom_point() +
  facet_wrap(~ as.factor(block)) +
  theme_bw() +
  scale_colour_brewer(palette = "Set1", name = "PFER threshold") +
  theme(legend.position = "top") +
  xlab("Time point") +
  ylab("Number of selected edges")
ggsave("Results/Network/Summary_multiomic_network_by_PFER.pdf", p)
