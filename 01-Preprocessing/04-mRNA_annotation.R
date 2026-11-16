# =============================================================================
# 04-mRNA_annotation.R
#
# Maps the microarray probes to gene symbols and names using the Agilent
# SurePrint G3 Human Gene Expression platform annotation from GEO (GPL16699).
#
# Input:  New_data/mRNA_denoised_252.rds (03-Denoising.R)
# Output: Results/mRNA_annotation/annotation.matched.rds (one row per probe and gene)
#         Results/mRNA_annotation/z.table.matched.rds    (aligned with the transcripts)
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("OXFORD_DIR", unset = ".")) # project root (see README)

library(GEOquery)

mRNA <- readRDS("New_data/mRNA_denoised_252.rds")
probes <- colnames(mRNA)

# Platform annotation
z <- getGEO("GPL16699")
z.table <- Table(z)
z.table.unique <- z.table[!duplicated(z.table[, c("NAME", "GENE_SYMBOL")]), ]

# Annotation of the probes in the data, in the same order as the transcripts
z.table.matched <- z.table.unique[z.table.unique$NAME %in% probes, c("NAME", "GENE_SYMBOL", "GENE_NAME", "GO_ID")]
rownames(z.table.matched) <- z.table.matched$NAME
z.table.matched <- z.table.matched[colnames(mRNA), ]
rownames(z.table.matched) <- z.table.matched$NAME

dir.create("Results/mRNA_annotation", showWarnings = FALSE, recursive = TRUE)
saveRDS(z.table.unique, "Results/mRNA_annotation/annotation.matched.rds")
saveRDS(z.table.matched, "Results/mRNA_annotation/z.table.matched.rds")
