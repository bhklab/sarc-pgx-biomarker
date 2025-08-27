# -----------------------------------------------------------
# Gene-Drug Association and Meta-Analysis for STS Cell Lines
#
# This script evaluates gene–drug associations across 
# soft-tissue sarcoma (STS) cell lines by integrating 
# gene expression and drug response profiles.
#
#   Data Sources:
#     - CCLE/CTRP
#     - GDSCv2
#     - NCI Sarcoma Panel
#
#   Includes:
#     * Correlation analysis between drug response (AAC) 
#       and aligned RNA-seq expression data (per study)
#     * Filtering of low-quality drugs (missingness > 70%)
#     * Multiple testing correction (BH-adjusted p-values)
#     * Meta-analysis across all three datasets using 
#       random-effects model (via `metacor`)
#
#   Outputs:
#     - Study-specific gene–drug correlation results
#     - Meta-analysis of gene–drug associations
#     - CSV tables with significant (FDR < 0.15) results
#     - R objects (qs format) for downstream use
#
#   Assumes aligned gene expression and drug response data 
#   has already been preprocessed and saved in .qs format.
# -----------------------------------------------------------
##################################################################
## Load libraries
##################################################################

library(qs)
library(PharmacoGx)
library(data.table)
library(magrittr)
library(meta)
library(metafor)

##################################################################
## Setup directory
##################################################################

dir_data <- 'data/procdata'  
dir_drug <- 'data/results/drug' 
dir_aligned <- 'data/results/aligned' 
dir_out <- 'data/results/drug'

#######################################################
## Load gene-drug data (after alignment)
#######################################################

dat <- qread( file.path(dir_drug, "drug_rna_aligned_sts.qs"))
pset_aac <- dat$pset_aac 
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

#######################################################
## Gene-drug association analysis (CCLE/CTRP) 
#######################################################
id <- colnames(pset_aac )[grep("-ccle", colnames(pset_aac ))]
pset_aac <- pset_aac[, id]
pset_aligned <- pset_aligned[, id]

# remove drugs with more than 70% missing values across CLs
missing_aac <- sapply(1:nrow(pset_aac), function(k){
  length(which(is.na(pset_aac[k,])))
})

missing_df <- data.frame(perc = round((missing_aac/ncol(pset_aac)) * 100), drug= rownames(pset_aac)) 
missing_df <- missing_df[missing_df$perc < 70, ]
pset_aac <- pset_aac[rownames(pset_aac) %in% missing_df$drug, ] # 69 drugs


gene_drug_res <- lapply(1:nrow(pset_aac), function(k){ 
  
  print(k)
  dat_aac <- pset_aac[k, ]
  
  res <- lapply(1:nrow(pset_aligned), function(i){
    
    fit <- cor.test(as.numeric(pset_aligned[i,]),as.numeric( dat_aac))

    data.frame(Ensembl_ID = gene_ann[gene_ann$Symbol == rownames(pset_aligned)[i], "EnsemblGeneId"],
               gene_name = rownames(pset_aligned)[i], 
               drug = rownames(pset_aac)[k],
               estimate = fit$estimate,
               df = fit$parameter + 2,
               lower = fit$conf.int[1],
               upper = fit$conf.int[2],
               pvalue = fit$p.value)
    
  })  
  
  cor_res <- do.call(rbind, res)
  cor_res$padj <- p.adjust(cor_res$pvalue, method = "BH")
  cor_res
  
})

gene_drug_assoc <- do.call(rbind, gene_drug_res)
rownames(gene_drug_assoc) <- NULL

qsave(gene_drug_assoc, file= file.path(dir_out, "gene_drug_assoc_sts_ccle_ctrp.qs"))
write.csv(gene_drug_assoc, file = file.path(dir_out, "gene_drug_assoc_sts_ccle_ctrp.csv"), row.names = FALSE)

sig_res <- gene_drug_assoc[gene_drug_assoc$padj < 0.15, ]
write.csv(sig_res, file=file.path(dir_out, "gene_drug_assoc_sts_ccle_ctrp_sigFDR.csv"))

#######################################################
## Gene-drug association analysis (GDSCv2) 
#######################################################
pset_aac <- dat$pset_aac 
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

# GDSC cellines
id <- colnames(pset_aac )[grep("-gdsc", colnames(pset_aac ))]
pset_aac <- pset_aac[, id]
pset_aligned <- pset_aligned[, id]

# remove drugs with more than 70% missing
missing_aac <- sapply(1:nrow(pset_aac), function(k){
  length(which(is.na(pset_aac[k,])))
})

missing_df <- data.frame(perc = round((missing_aac/ncol(pset_aac)) * 100), drug= rownames(pset_aac)) 
missing_df <- missing_df[missing_df$perc < 70, ]

pset_aac <- pset_aac[rownames(pset_aac) %in% missing_df$drug, ] 

gene_drug_res <- lapply(1:nrow(pset_aac), function(k){ 
  
  print(k)
  dat_aac <- pset_aac[k, ]
  
  res <- lapply(1:nrow(pset_aligned), function(i){
    
    fit <- cor.test(as.numeric(pset_aligned[i,]), as.numeric(dat_aac))
    data.frame(Ensembl_ID = gene_ann[gene_ann$Symbol == rownames(pset_aligned)[i], "EnsemblGeneId"],
               gene_name = rownames(pset_aligned)[i], 
               drug = rownames(pset_aac)[k],
               estimate = fit$estimate,
               df = fit$parameter + 2,
               lower = fit$conf.int[1],
               upper = fit$conf.int[2],
               pvalue = fit$p.value)
    
  })  
  
  cor_res <- do.call(rbind, res)
  cor_res$padj <- p.adjust(cor_res$pvalue, method = "BH")
  cor_res
  
})

gene_drug_assoc <- do.call(rbind, gene_drug_res)
rownames(gene_drug_assoc) <- NULL

qsave(gene_drug_assoc, file=file.path(dir_out, "gene_drug_assoc_sts_gdsc.qs"))
write.csv(gene_drug_assoc, file = file.path(dir_out, "gene_drug_assoc_sts_gdsc.csv"), row.names = FALSE)

sig_res <- gene_drug_assoc[gene_drug_assoc$padj < 0.15, ]
write.csv(sig_res, file=file.path(dir_out, "gene_drug_assoc_sts_gdsc_sigFDR.csv"))

#######################################################
## Gene-drug association analysis (NCI-Sarcoma) 
#######################################################

pset_aac <- dat$pset_aac 
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

# NCI-sarcoma cellines
id <- colnames(pset_aac )[grep("-nci", colnames(pset_aac ))]
pset_aac <- pset_aac[, id]
pset_aligned <- pset_aligned[, id]

# remove drugs with more than 70% missing
missing_aac <- sapply(1:nrow(pset_aac), function(k){
  length(which(is.na(pset_aac[k,])))
})

missing_df <- data.frame(perc = round((missing_aac/ncol(pset_aac)) * 100), drug= rownames(pset_aac)) 
missing_df <- missing_df[missing_df$perc < 70, ]

pset_aac <- pset_aac[rownames(pset_aac) %in% missing_df$drug, ] 

gene_drug_res <- lapply(1:nrow(pset_aac), function(k){ 
  
  print(k)
  dat_aac <- pset_aac[k, ]
  
  res <- lapply(1:nrow(pset_aligned), function(i){
    
    fit <- cor.test(as.numeric(pset_aligned[i,]), as.numeric(dat_aac))
    data.frame(Ensembl_ID = gene_ann[gene_ann$Symbol == rownames(pset_aligned)[i], "EnsemblGeneId"],
               gene_name = rownames(pset_aligned)[i], 
               drug = rownames(pset_aac)[k],
               estimate = fit$estimate,
               df = fit$parameter + 2,
               lower = fit$conf.int[1],
               upper = fit$conf.int[2],
               pvalue = fit$p.value)
    
  })  
  
  cor_res <- do.call(rbind, res)
  cor_res$padj <- p.adjust(cor_res$pvalue, method = "BH")
  cor_res
  
})

gene_drug_assoc <- do.call(rbind, gene_drug_res)
rownames(gene_drug_assoc) <- NULL

qsave(gene_drug_assoc, file=file.path(dir_out, "gene_drug_assoc_sts_nci.qs"))
write.csv(gene_drug_assoc, file = file.path(dir_out, "gene_drug_assoc_sts_nci.csv"), row.names = FALSE)

sig_res <- gene_drug_assoc[gene_drug_assoc$padj < 0.15, ]
write.csv(sig_res, file=file.path(dir_out, "gene_drug_assoc_sts_nci_sigFDR.csv"))

####################################################################
## Integration analysis
####################################################################
res.ccle <- qread(file=file.path(dir_out, "gene_drug_assoc_sts_ccle_ctrp.qs"))   
res.gdsc <- qread(file=file.path(dir_out, "gene_drug_assoc_sts_gdsc.qs"))  
res.nci <- qread(file=file.path(dir_out, "gene_drug_assoc_sts_nci.qs"))  

res.ccle$study <- "CCLE/CTRP"
res.gdsc$study <- "GDSC"
res.nci$study <- "NCI"
res <- rbind(res.ccle, res.gdsc, res.nci)

drugs <- intersect(intersect(res.ccle$drug, res.gdsc$drug),
                   res.nci$drug)
genes <- intersect(intersect(res.ccle$gene_name, res.gdsc$gene_name),
                   res.nci$gene_name)

meta.res.drug <- lapply(1:length(drugs), function(k){ 
  
  print(k)
  df <- res[res$drug == drugs[k], ]
  meta.res <- lapply(1:length(genes), function(i){
    
    print(i)
    df.genes <- df[df$gene_name == genes[i], ]
    
    fit <- metacor(df.genes$estimate, 
                   df.genes$df, 
                   sm = "cor",
                   control=list(stepadj=0.5, maxiter=1000))
    
   data.frame(gene_name = unique(df.genes$gene_name),
              Ensembl_ID = unique(df.genes$Ensembl_ID),
              r = fit$TE.random,
              se = fit$seTE.random,
              pval = fit$pval.random,
              I2 = fit$I2)
    
  })
  
  meta.res <- do.call(rbind, meta.res)
  meta.res <-  meta.res[!is.na(meta.res$r), ]
  meta.res$padj <- p.adjust(meta.res$pval, method = "BH")
  meta.res$drug <- drugs[k]
  meta.res
  
})

meta.res.drug <- do.call(rbind, meta.res.drug)
sig_res <- meta.res.drug[meta.res.drug$padj < 0.15, ]
write.csv(sig_res, file= file.path(dir_out, "gene_drug_assoc_sts_meta_sigFDR.csv"), row.names = FALSE)
qsave(meta.res.drug, file= file.path(dir_out, "gene_drug_assoc_sts_meta.qs"))


