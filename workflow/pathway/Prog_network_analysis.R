##--------------------------------------------------------------------------
## Pathway Enrichment Analysis for Gene–Drug Associations (STS)
##
## Description:
## This script performs comprehensive pathway enrichment analyses using 
## gene–drug association results (e.g., correlation of gene expression with 
## drug sensitivity) from a soft tissue sarcoma (STS) meta-analysis.
##
## Specifically, the script carries out:
## 1. Gene Set Enrichment Analysis (GSEA) using ranked correlation values (r) and p-value
##    for each drug against three curated pathway databases (MSigDB):
##      - HALLMARK 
##      - GO Biological Processes
##
## 2. Over-Representation Analysis (ORA) using significantly associated genes 
##    (FDR < 0.1) from each drug's profile against three curated pathway databases (MSigDB):
##      - HALLMARK 
##      - GO Biological Processes
##
## Output:
## - GSEA and ORA results per drug stored as `.qs` and `.csv` files
## - Pathway data object saved for future reuse
##
## Dependencies:
## - Libraries: qs, data.table, magrittr, fgsea, reshape2
## - Pathway definitions: MSigDB GMT files (HALLMARK, GO)
##---------------------------------------------------------------------------
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






















  set.seed(123)
  fgseaRes <- fgsea(hallmark_pathway, 
                    ranks, 
                    minSize=15, 
                    maxSize = 500)
  
  fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
  data.frame(drug = drugs[i], 
             fgseaRes[order(fgseaRes$padj, fgseaRes$NES, decreasing = c(FALSE, TRUE)), ])
  
})

gsea_hallmark <- do.call(rbind, gsea_hallmark)
gsea_hallmark <- gsea_hallmark[!is.na(gsea_hallmark$pval), ]

qsave(gsea_hallmark, file= file.path(dir_out, "hallmark_gsea_pathway_drug.qs"))
write.csv(gsea_hallmark, file = file.path(dir_out, "hallmark_gsea_pathway_drug.csv"), row.names = FALSE)

###################################################
## GSEA: STS gene association result and GO
###################################################

gsea_go <- lapply(1:length(drugs), function(i){

  df <- dat[dat$drug == drugs[i], ] 
  #ranks <- df$r
  ranks <- sign(df$r) * (-log10(df$pval))
  names(ranks) <- df$gene_name
  ranks <- sort(ranks, decreasing = T)
  
  set.seed(123)
  fgseaRes <- fgsea(go_pathway, 
                    ranks, 
                    minSize=15, 
                    maxSize = 500)
  
  fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
  fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
  data.frame(drug = drugs[i], 
             fgseaRes[order(fgseaRes$padj, fgseaRes$NES, decreasing = c(FALSE, TRUE)), ])
  
})

gsea_go <- do.call(rbind, gsea_go)
gsea_go <- gsea_go[!is.na(gsea_go$pval), ]

qsave(gsea_go, file= file.path(dir_out, "go_gsea_pathway_drug.qs"))
write.csv(gsea_go, file = file.path(dir_out, "go_gsea_pathway_drug.csv"), row.names = FALSE)


################################################################################
## ORA: STS gene association result and HALLMARK
################################################################################
dat <- qread(file= file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
dat <- dat[!is.na(dat$r), ]
drugs <- unique(dat$drug)

ora_hallmark <- lapply(1:length(drugs), function(i){
  
  print(i)
  df <- dat[dat$drug == drugs[i], ] 
  df <- df[!is.na(df$r), ]
  universe <- df$gene_name 
  genes <- df[df$padj < thr_ora, "gene_name"]
  
  ora_res <- fora(pathways = hallmark_pathway,
                 genes    = genes,
                 universe = universe,
                 minSize = 15,
                 maxSize = 500)
  
  ora_res  <- as.data.frame(ora_res)

  if (!nrow(ora_res)){

   update_ora_res <- as.data.frame(t(rep(NA, length(colnames(ora_res)))))
   colnames(update_ora_res) <- colnames(ora_res)
   data.frame(drug = drugs[i], update_ora_res)

  }else{

   ora_res  <- ora_res[!is.na(ora_res$padj), ]
   ora_res$overlapGenes <- sapply(ora_res$overlapGenes, paste, collapse = ",")
   data.frame(drug = drugs[i],  ora_res[order(ora_res$padj, ora_res$pval), ])

  }
  
  
})

ora_hallmark <- dplyr::bind_rows(ora_hallmark)
ora_hallmark <- ora_hallmark[!is.na(ora_hallmark$pval), ]

qsave(ora_hallmark, file= file.path(dir_out, "hallmark_ora_pathway_drug.qs"))
write.csv(ora_hallmark, file = file.path(dir_out, "hallmark_ora_pathway_drug.csv"), row.names = FALSE)

################################################################################
## ORA: STS gene association result and Go
################################################################################

ora_go <- lapply(1:length(drugs), function(i){
  
  print(i)
  df <- dat[dat$drug == drugs[i], ] 
  df <- df[!is.na(df$r), ]
  universe <- df$gene_name 
  genes <- df[df$padj < thr_ora, "gene_name"]
  
  ora_res <- fora(pathways = go_pathway,
                 genes    = genes,
                 universe = universe,
                 minSize = 15,
                 maxSize = 500)
  
  ora_res  <- as.data.frame(ora_res)

  if (!nrow(ora_res)){

   update_ora_res <- as.data.frame(t(rep(NA, length(colnames(ora_res)))))
   colnames(update_ora_res) <- colnames(ora_res)
   data.frame(drug = drugs[i], update_ora_res)

  }else{

   ora_res  <- ora_res[!is.na(ora_res$padj), ]
   ora_res$overlapGenes <- sapply(ora_res$overlapGenes, paste, collapse = ",")
   data.frame(drug = drugs[i],  ora_res[order(ora_res$padj, ora_res$pval), ])

  }
  
})

ora_go <- dplyr::bind_rows(ora_go)
ora_go <- ora_go[!is.na(ora_go$pval), ]

qsave(ora_go, file= file.path(dir_out, "go_ora_pathway_drug.qs"))
write.csv(ora_go, file = file.path(dir_out, "go_ora_pathway_drug.csv"), row.names = FALSE)

