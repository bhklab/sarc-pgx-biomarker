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
library(data.table)
library(magrittr)
library(fgsea)
library(reshape2)
library(clusterProfiler)

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/results/drug' 
dir_pathway <- 'data/rawdata' 
dir_out <- 'data/results/pathway' 

thr_ora <- 0.15
###################################################################
## load GO and HALLMARK pathways downloaded from msigdb
##################################################################

hallmark_pathway <- gmtPathways(file.path(dir_pathway, "h.all.v2025.1.Hs.symbols.gmt"))
go_pathway <- gmtPathways(file.path(dir_pathway, "c5.go.bp.v2025.1.Hs.symbols.gmt"))

dat <- list(HALLMARK = hallmark_pathway, GO = go_pathway)
qsave(dat, file=  file.path(dir_pathway, "pathway_data.qs"))

###################################################
## GSEA: STS gene association result and HALLMARK
###################################################

dat <- qread(file= file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
dat <- dat[!is.na(dat$r), ]
drugs <- unique(dat$drug)

gsea_hallmark <- lapply(1:length(drugs), function(i){
  
  df <- dat[dat$drug == drugs[i], ] 
  #ranks <- df$r 
  ranks <- sign(df$r) * (-log10(df$pval))
  names(ranks) <- df$gene_name
  ranks <- sort(ranks, decreasing = T)
  
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

# Convert list to TERM2GENE data frame
TERM2GENE <- stack(hallmark_pathway)         
colnames(TERM2GENE) <- c("gene", "term")      
TERM2GENE <- TERM2GENE[, c("term", "gene")]   

ora_hallmark <- lapply(1:length(drugs), function(i){
  
  print(i)
  df <- dat[dat$drug == drugs[i], ] 
  df <- df[!is.na(df$r), ]
  universe <- df$gene_name 
  genes <- df[df$padj < thr_ora, "gene_name"]
  
  ora_res <- enricher(gene = genes,
                      TERM2GENE = TERM2GENE,
                      universe = universe,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      minGSSize = 10,
                      maxGSSize = 500)
  
  ora_res  <- as.data.frame(ora_res)

  if (!nrow(ora_res)){

   update_ora_res <- as.data.frame(t(rep(NA, length(colnames(ora_res)))))
   colnames(update_ora_res) <- colnames(ora_res)
   data.frame(drug = drugs[i], update_ora_res)

  }else{

   ora_res  <- ora_res[!is.na(ora_res $p.adjust), ]
   data.frame(drug = drugs[i], ora_res)

  }
  
  
})

ora_hallmark <- dplyr::bind_rows(ora_hallmark)
ora_hallmark <- ora_hallmark[!is.na(ora_hallmark$pvalu), ]

qsave(ora_hallmark, file= file.path(dir_out, "hallmark_ora_pathway_drug.qs"))
write.csv(ora_hallmark, file = file.path(dir_out, "hallmark_ora_pathway_drug.csv"), row.names = FALSE)

################################################################################
## ORA: STS gene association result and KEGG
################################################################################
# Convert list to TERM2GENE data frame
TERM2GENE <- stack(go_pathway)         
colnames(TERM2GENE) <- c("gene", "term")      
TERM2GENE <- TERM2GENE[, c("term", "gene")]   

ora_go <- lapply(1:length(drugs), function(i){
  
  df <- dat[dat$drug == drugs[i], ] 
  universe <- df$gene_name 
  genes <- df[df$padj < thr_ora, "gene_name"]
  
  ora_res <- enricher(gene = genes,
                      TERM2GENE = TERM2GENE,
                      universe = universe,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1,
                      minGSSize = 10,
                      maxGSSize = 500)
  
  ora_res  <- as.data.frame(ora_res)

  if (!nrow(ora_res)){

   update_ora_res <- as.data.frame(t(rep(NA, length(colnames(ora_res)))))
   colnames(update_ora_res) <- colnames(ora_res)
   data.frame(drug = drugs[i], update_ora_res)

  }else{

   ora_res  <- ora_res[!is.na(ora_res $p.adjust), ]
   data.frame(drug = drugs[i], ora_res)

  }
  
})

ora_go <- dplyr::bind_rows(ora_go)
ora_go <- ora_go[!is.na(ora_go$pval), ]

qsave(ora_go, file= file.path(dir_out, "go_ora_pathway_drug.qs"))
write.csv(ora_go, file = file.path(dir_out, "go_ora_pathway_drug.csv"), row.names = FALSE)

