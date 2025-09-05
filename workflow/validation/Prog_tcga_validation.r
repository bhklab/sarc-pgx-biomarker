#-----------------------------------------------------------------------
#
#
#
#-----------------------------------------------------------------------
########################################################################
## Load library
########################################################################

library(qs)
library(survival)
library(data.table)
library(magrittr)
library(fgsea)
library(reshape2)

########################################################################
## Set up directroy
########################################################################

dir_in <- 'data/procdata/validation'
dir_out <- 'data/results/validation/TCGA'
dir_pathway <-  'data/rawdata' 

time.censor <- 36
thr_ora <- 0.15

#######################################################################
## Load TCGA-SARC data
#######################################################################
dat <- qread(file.path(dir_in, 'curated_TCGA_se.qs'))

expr <- assay(dat)
clin <- data.frame(colData(dat))
annot <- data.frame(rowData(dat))

########################################################################
## Fit proportional hazard model
########################################################################

res <- lapply(1:nrow(expr), function(k){

data <- data.frame( status=clin$os_event , time=clin$os_time , variable=expr[k, ] )
data <- data[!is.na(data$variable), ]
data$time <- as.numeric(as.character(data$time))
data$variable <- as.numeric( as.character(data$variable) )
  
for(i in 1:nrow(data)){
    
    if( !is.na(as.numeric(as.character(data[ i , "time" ]))) && as.numeric(as.character(data[ i , "time" ])) > time.censor ){
      
      data[ i , "time" ] <- time.censor
      data[ i , "status" ] <- 0
      
    }
  }
  
  cox <- coxph( Surv( time , status ) ~ variable, data=data )
  
  hr <- summary(cox)$coefficients[, "coef"]
  se <- summary(cox)$coefficients[, "se(coef)"]
  n <- round(summary(cox)$n)
  low <- summary(cox)$conf.int[, "lower .95"]
  up <- summary(cox)$conf.int[, "upper .95"]
  pval <- summary(cox)$coefficients[, "Pr(>|z|)"]
  
data.frame(gene_name = rownames(expr)[k],
           HR = hr,
           SE = se,
           N = n,
           Low = low,
           Up = up,
           pval = pval)

})

res.cox <- do.call(rbind, res)
res.cox <- res.cox[!is.na(res.cox$HR), ]
res.cox$padj <- p.adjust(res.cox$pval, method = 'BH')

save(res.cox, file = file.path(dir_out, 'cox_tcga.rda'))

########################################################################
## ORA pathway analysis
########################################################################
## load GO and HALLMARK pathways downloaded from msigdb

hallmark_pathway <- gmtPathways(file.path(dir_pathway, "h.all.v2025.1.Hs.symbols.gmt"))
go_pathway <- gmtPathways(file.path(dir_pathway, "c5.go.bp.v2025.1.Hs.symbols.gmt"))

###################################################
## GSEA and HALLMARK
###################################################

dat <- res.cox
dat <- dat[!is.na(dat$HR), ]

ranks <- sign(dat$HR) * (-log10(dat$pval))
names(ranks) <- dat$gene_name
ranks <- sort(ranks, decreasing = T)
  
set.seed(123)
fgseaRes <- fgsea(hallmark_pathway, 
                  ranks, 
                  minSize=15, 
                  maxSize = 500)
  
fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
fgseaRes <- fgseaRes[order(fgseaRes$padj, fgseaRes$NES, decreasing = c(FALSE, TRUE)), ]  

qsave(fgseaRes, file= file.path(dir_out, "hallmark_gsea_pathway.qs"))
write.csv(fgseaRes, file = file.path(dir_out, "hallmark_gsea_pathway.csv"), row.names = FALSE)

###################################################
## GSEA and GO:BP
###################################################

dat <- res.cox
dat <- dat[!is.na(dat$HR), ]

ranks <- sign(dat$HR) * (-log10(dat$pval))
names(ranks) <- dat$gene_name
ranks <- sort(ranks, decreasing = T)
  
set.seed(123)
fgseaRes <- fgsea(go_pathway, 
                  ranks, 
                  minSize=15, 
                  maxSize = 500)
  
fgseaRes <- fgseaRes[!is.na(fgseaRes$padj), ]
fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
fgseaRes <- fgseaRes[order(fgseaRes$padj, fgseaRes$NES, decreasing = c(FALSE, TRUE)), ]  

qsave(fgseaRes, file= file.path(dir_out, "go_gsea_pathway.qs"))
write.csv(fgseaRes, file = file.path(dir_out, "go_gsea_pathway.csv"), row.names = FALSE)

################################################################################
## ORA and HALLMARK
################################################################################

dat <- res.cox
dat <- dat[!is.na(dat$HR), ]
universe <- dat$gene_name 
genes <- dat[dat$padj < thr_ora, "gene_name"]
  
ora_res <- fora(pathways = hallmark_pathway,
                genes    = genes,
                universe = universe,
                minSize = 15,
                maxSize = 500)
  
ora_res  <- as.data.frame(ora_res)
ora_res  <- ora_res[!is.na(ora_res$padj), ]
ora_res$overlapGenes <- sapply(ora_res$overlapGenes, paste, collapse = ",")
ora_res <- ora_res[order(ora_res$padj, ora_res$pval), ]

qsave(ora_res, file= file.path(dir_out, "hallmark_ora_pathway.qs"))
write.csv(ora_res, file = file.path(dir_out, "hallmark_ora_pathway.csv"), row.names = FALSE)

################################################################################
## ORA and Go
################################################################################

dat <- res.cox
dat <- dat[!is.na(dat$HR), ]
universe <- dat$gene_name 
genes <- dat[dat$padj < thr_ora, "gene_name"]
  
ora_res <- fora(pathways = go_pathway,
                genes    = genes,
                universe = universe,
                minSize = 15,
                maxSize = 500)
  
ora_res  <- as.data.frame(ora_res)
ora_res  <- ora_res[!is.na(ora_res$padj), ]
ora_res$overlapGenes <- sapply(ora_res$overlapGenes, paste, collapse = ",")
ora_res <- ora_res[order(ora_res$padj, ora_res$pval), ]

qsave(ora_res, file= file.path(dir_out, "go_ora_pathway.qs"))
write.csv(ora_res, file = file.path(dir_out, "go_ora_pathway.csv"), row.names = FALSE)


