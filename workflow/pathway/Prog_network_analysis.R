##--------------------------------------------------------------------------
## Pathway Enrichment Analysis for Gene–Drug Associations in STS
##
## This script analyzes gene–drug association results from a soft tissue 
## sarcoma (STS) meta-analysis using:
##   1. Gene Set Enrichment Analysis (GSEA) with ranked correlations
##   2. Over-Representation Analysis (ORA) with significant genes (FDR < 0.1)
##
## Databases: MSigDB HALLMARK and GO Biological Processes.
##
## Outputs:
## - GSEA and ORA results per drug (.qs, .csv)
## - Ranked gene lists (.rnk) and enrichment tables for Cytoscape/EnrichmentMap
## - Pathway label files for network visualization
##
## Dependencies: qs, dplyr, purrr, fgsea, GSEABase, stringr
##--------------------------------------------------------------------------
########################################
## Load libraries
########################################

library(qs)
library(dplyr)
library(readr)
library(purrr)
library(GSEABase)
library(stringr)

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/results/drug' 
dir_pathway <- 'data/rawdata' 
dir_out <- 'data/results' 

thr_gsea <- 0.001

###################################################
## Load pathway data
###################################################

hallmark_pathway <- gmtPathways(file.path(dir_pathway, "h.all.v2025.1.Hs.symbols.gmt"))
go_pathway <- gmtPathways(file.path(dir_pathway, "c5.go.bp.v2025.1.Hs.symbols.gmt"))

###################################################
## GSEA: GO:BP result
###################################################
# load selected drugs
selected_drug <- read.csv(file.path('data/results/drug', 'selected_drugs.csv'))
selected_drug_sts <- selected_drug[selected_drug$sts == 'Yes', 'drug']

# load GSEA results using GO:BP 
dat <- qread(file.path(dir_out, 'pathway', 'go_gsea_pathway_drug.qs'))
dat <- dat[!is.na(dat$pval), ]
dat$pathway_update <- tools::toTitleCase(tolower(substr(dat$pathway, 6, nchar(dat$pathway))))
lab <- data.frame(
  Name = dat$pathway,  # must match the EM node "name"
  pretty_name = dat$pathway_update
)
write.table(lab, file = file.path(dir_out, 'network/GO', "labels.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

sig <- dat[dat$drug %in% selected_drug_sts & dat$padj < thr_gsea, ]

drug <- unique(sig$drug)
for(i in 1:length(drug)){

df <-  sig[sig$drug == drug[i], c("pathway", "pval", "padj", "NES")]

write.table(df, file = file.path(dir_out, 'network/GO/sig_0.001', paste(drug[i], 'gsea_go.txt', sep=".")), 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

df <- dat[dat$drug == drug[i], c("pathway", "pval", "padj", "NES")]

write.table(df, file = file.path(dir_out, 'network/GO/GSEA', paste(drug[i], 'gsea_go.txt', sep=".")), 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

}

###################################################
## GSEA: STS gene association result and GO
###################################################
# Ranked list (.rnk)
dat <- qread(file= file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
dat <- dat[!is.na(dat$r), ]
drug <- intersect(unique(dat$drug), selected_drug_sts)

for(i in 1:length(drug)){

df <- dat[dat$drug == drug[i], ] 
ranks <- sign(df$r) * (-log10(df$pval))
names(ranks) <- df$gene_name
ranks <- sort(ranks, decreasing = T)

# Two columns: Gene <tab> Score. Keep header (EM can handle it).
rnk <- tibble(
  Gene  = names(ranks),
  Score = as.numeric(ranks)
) %>%
  arrange(desc(Score))

write_tsv(rnk,  file = file.path(dir_out, 'network/ranks', paste(drug[i], 'metaCorRanks.rnk', sep='_')))  

}





