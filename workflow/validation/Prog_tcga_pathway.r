#--------------------------------------------------------------------------------
## Description:
#   This script performs pathway enrichment on gene-level survival statistics from
#   TCGA-SARC using both:
#     Late integration (LI): per-histology Cox results pooled by random-effects 
#     meta-analysis (pooled log-HR and p-values)
#
#   For each integration type, the script runs:
#     • GSEA (fgsea) against MSigDB Hallmark and GO:BP collections using a
#       preranked list: rank = sign(log-HR) × (−log10 p-value).
#       - LI ranks use columns: Coef (pooled log-HR), Pval
#     • ORA (fora) on the set of genes with BH-FDR < thr_ora
#       - Universe = all genes available in the respective LI result table
#       - Pathway size filters: 15 ≤ size ≤ 500
#
# Required inputs
#   • Gene-level results (produced by a separate script):
#       data/results/validation/TCGA/meta_cox_tcga_histo.rda
#         - Data frame `meta.cox` with: Gene, Coef (pooled log-HR), SE,
#           CI_lower, CI_upper, Pval, padj, I2, Q_Pval
#   • MSigDB GMT files in data/rawdata/:
#       h.all.v2025.1.Hs.symbols.gmt        (Hallmark)
#       c5.go.bp.v2025.1.Hs.symbols.gmt     (GO: Biological Process)
#
# Notes
#   • Ensure pathway gene symbols match those used in the LI result tables.
#   • Output directories must exist prior to running.
#--------------------------------------------------------------------------------
########################################################################
## Load library
########################################################################

library(qs)
library(data.table)
library(magrittr)
library(fgsea)
library(reshape2)
library(SummarizedExperiment)

########################################################################
## Set up directroy
########################################################################

dir_in <- 'data/procdata/validation'
dir_out <- 'data/results/validation/TCGA'
dir_pathway <-  'data/rawdata' 

thr_ora <- 0.15
thr_gsea <- 0.05

############################################################
## load GO and HALLMARK pathways downloaded from msigdb
############################################################

hallmark_pathway <- gmtPathways(file.path(dir_pathway, "h.all.v2025.1.Hs.symbols.gmt"))
go_pathway <- gmtPathways(file.path(dir_pathway, "c5.go.bp.v2025.1.Hs.symbols.gmt"))

########################################################################
## Pathway analysis (late integration)
## GSEA and HALLMARK
########################################################################
load(file.path(dir_out, 'meta_cox_tcga_histo.rda'))

dat <- meta.cox
dat <- dat[!is.na(dat$Coef), ]

ranks <- sign(dat$Coef) * (-log10(dat$Pval))
names(ranks) <- dat$Gene
ranks <- sort(ranks, decreasing = T)
  
set.seed(123)
fgseaRes <- fgsea(hallmark_pathway, 
                  ranks, 
                  minSize=15, 
                  maxSize = 500)
  
fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
fgseaRes <- fgseaRes[order(fgseaRes$padj, fgseaRes$NES, decreasing = c(FALSE, TRUE)), ]  
sig <- fgseaRes[fgseaRes$padj < 0.10, ]

qsave(fgseaRes, file= file.path(dir_out, "hallmark_gsea_pathway.qs"))
write.csv(fgseaRes, file = file.path(dir_out, "hallmark_gsea_pathway.csv"), row.names = FALSE)
write.csv(sig, file = file.path(dir_out, "hallmark_gsea_pathway_0.1.csv"), row.names = FALSE)

###################################################
## GSEA and GO:BP
###################################################

dat <- meta.cox
dat <- dat[!is.na(dat$Coef), ]

ranks <- sign(dat$Coef) * (-log10(dat$Pval))
names(ranks) <- dat$Gene
ranks <- sort(ranks, decreasing = T)
  
set.seed(123)
fgseaRes <- fgsea(go_pathway, 
                  ranks, 
                  minSize=15, 
                  maxSize = 500)
  
fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
fgseaRes <- fgseaRes[order(fgseaRes$padj, fgseaRes$NES, decreasing = c(FALSE, TRUE)), ]  
sig <- fgseaRes[fgseaRes$padj < 0.1, ]

qsave(fgseaRes, file= file.path(dir_out, "go_gsea_pathway.qs"))
write.csv(fgseaRes, file = file.path(dir_out, "go_gsea_pathway.csv"), row.names = FALSE)
write.csv(sig, file = file.path(dir_out, "go_gsea_pathway_0.1.csv"), row.names = FALSE)

################################################################################
## ORA and HALLMARK
################################################################################

dat <- meta.cox
dat <- dat[!is.na(dat$Coef), ]
universe <- dat$Gene
genes <- dat[dat$padj < thr_ora, "Gene"]
  
ora_res <- fora(pathways = hallmark_pathway,
                genes    = genes,
                universe = universe,
                minSize = 15,
                maxSize = 500)
  
ora_res  <- as.data.frame(ora_res)
ora_res  <- ora_res[!is.na(ora_res$padj), ]
ora_res$overlapGenes <- sapply(ora_res$overlapGenes, paste, collapse = ",")
ora_res <- ora_res[order(ora_res$padj, ora_res$pval), ]
sig <- ora_res[ora_res$padj < 0.15, ]

qsave(ora_res, file= file.path(dir_out, "hallmark_ora_pathway.qs"))
write.csv(ora_res, file = file.path(dir_out, "hallmark_ora_pathway.csv"), row.names = FALSE)
write.csv(sig, file = file.path(dir_out, "hallmark_ora_pathway_0.15.csv"), row.names = FALSE)

################################################################################
## ORA and Go
################################################################################

dat <- meta.cox
dat <- dat[!is.na(dat$Coef), ]
universe <- dat$Gene 
genes <- dat[dat$padj < thr_ora, "Gene"]
  
ora_res <- fora(pathways = go_pathway,
                genes    = genes,
                universe = universe,
                minSize = 15,
                maxSize = 500)
  
ora_res  <- as.data.frame(ora_res)
ora_res  <- ora_res[!is.na(ora_res$padj), ]
ora_res$overlapGenes <- sapply(ora_res$overlapGenes, paste, collapse = ",")
ora_res <- ora_res[order(ora_res$padj, ora_res$pval), ]
sig <- ora_res[ora_res$padj < 0.15, ]

qsave(ora_res, file= file.path(dir_out, "go_ora_pathway.qs"))
write.csv(ora_res, file = file.path(dir_out, "go_ora_pathway.csv"), row.names = FALSE)
write.csv(sig, file = file.path(dir_out, "go_ora_pathway_0.15.csv"), row.names = FALSE)
