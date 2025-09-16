#--------------------------------------------------------------------------------
## Description:
#   This script performs gene-level survival analysis on TCGA-SARC RNA-seq with
#   two strategies and saves results for downstream pathway enrichment:
#
#   (1) Early integration (EI): Pooled Cox model across all histologies
#       without a histology term (i.e., no subtype correction).
#
#   (2) Late integration (LI): Per-histology Cox models, followed by a
#       random-effects meta-analysis across histologies (DerSimonian–Laird
#       with Hartung–Knapp adjustment) to obtain a pooled log-HR and heterogeneity.
#
# Required input:
#   - SummarizedExperiment (qs): data/procdata/validation/curated_TCGA_se.qs
#       * assay(dat): gene expression matrix (genes × samples)
#       * colData(dat): clinical data with os_time (months) and os_event (0/1),
#                       histological (histology labels), patientid (sample ID)
#       * rowData(dat): gene annotations (row names must match expression genes)
#
# Outputs (written to data/results/validation/TCGA/):
#   - cox_tcga.rda / cox_tcga.csv
#       * EI (pooled) per-gene results across all histologies (no adjustment).
#       * Columns: gene_name, HR (log-HR; Cox coef), SE, N,
#                  Low/Up (95% CI on HR scale = exp(log-HR ± 1.96*SE)), pval, padj.
#
#   - cox_tcga_histo.rda / cox_tcga_histo.csv
#       * Per-histology (subtype) per-gene results.
#       * Columns: histo, gene_name, HR (log-HR), SE, N,
#                  Low/Up (HR scale), pval, padj.
#
#   - meta_cox_tcga_histo.rda / meta_cox_tcga_histo.csv
#       * Late-integration random-effects meta-analysis across histologies.
#       * Columns: Gene, Coef (pooled log-HR), SE, CI_lower/CI_upper (log-HR scale),
#                  Pval, padj, I2, Q_Pval.
#       * Requires ≥3 histology strata with estimable coefficients to run meta-analysis.
#
# Key analysis settings:
#   - Administrative censoring at 36 months (time.censor = 36):
#       times > 36 are truncated to 36 and status set to 0 at censoring.
#   - Model form: Surv(time, status) ~ expression (single-gene models; no other covariates).
#   - No histology correction in EI; LI operates via per-subtype fits → meta.
#
# Dependencies:
#   - qs, survival, data.table, magrittr, SummarizedExperiment, meta
#   - (fgsea/ORA packages are used only in the companion pathway script)
#
# Conventions / scale:
#   - “HR” columns in EI/per-histology tables are actually log-HR (Cox coefficients).
#     Low/Up are on the HR scale (exponentiated). In the meta table, Coef and CIs
#     are on the log-HR scale; exponentiate if you need HRs.
#--------------------------------------------------------------------------------
########################################################################
## Load library
########################################################################

library(qs)
library(survival)
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

time.censor <- 36
#######################################################################
## Load TCGA-SARC data
#######################################################################
dat <- qread(file.path(dir_in, 'curated_TCGA_se.qs'))

expr <- assay(dat)
clin <- data.frame(colData(dat))
annot <- data.frame(rowData(dat))
########################################################################
## Type 1: Fit proportional hazard model
########################################################################
## Integrate all subtypes and find association (early integration)
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
write.csv(res.cox, file = file.path(dir_out, "cox_tcga.csv"), row.names = FALSE)

########################################################################
## Type 2: Fit proportional hazard model 
########################################################################
## Step 1: Find association per subtype
group <- unique(clin$histological)

res.cox <- lapply(1:length(group), function(j){

print(j)
df.clin <- clin[clin$histological == group[j], ]
df.expr <- expr[, colnames(expr) %in% df.clin$patientid]

res <- lapply(1:nrow(expr), function(k){

data <- data.frame( status=df.clin$os_event , time=df.clin$os_time , variable=df.expr[k, ] )
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
  
data.frame(histo = group[j],
           gene_name = rownames(expr)[k],
           HR = hr,
           SE = se,
           N = n,
           Low = low,
           Up = up,
           pval = pval)

})

res <- do.call(rbind, res)
res <- res[!is.na(res$HR), ]
res$padj <- p.adjust(res$pval, method = 'BH')
res

})

res.cox <- do.call(rbind, res.cox)
save(res.cox, file = file.path(dir_out, 'cox_tcga_histo.rda'))
write.csv(res.cox, file = file.path(dir_out, "cox_tcga_histo.csv"), row.names = FALSE)

## Step 2: Integrate across subtypes
load(file.path(dir_out, 'cox_tcga_histo.rda'))
genes <- unique(res.cox$gene_name)

meta.cox <- laooly(1:length(genes), function(k){

data <- res.cox[res.cox$gene_name == genes[k], ]
data <- data.frame( Gene = genes[k],
                    Study = as.character( data$histo ),
                    N = as.numeric(as.character( data$N )),
                    Coef = as.numeric(as.character( data$HR )),
                    SE = as.numeric(as.character( data$SE )),
                    Pval = as.numeric(as.character( data$pval )))
  
  data <- data[ order( data$Coef ) , ]
  
  ## at least 3 studies needed to do the random effect meta-analyses
  if(nrow(data) >= 3){
    
    meta <- metagen( TE = Coef,
                     seTE = SE ,
                     data = data ,
                     studlab = Study ,
                     fixed = FALSE ,
                     random = TRUE ,
                     control = list( maxiter = 10000 , stepadj=0.5 ) )
    
    meta_res <- data.frame(Gene = genes[k],
                           Coef = meta$TE.random ,
                           SE = meta$seTE.random ,
                           CI_lower = meta$lower.random ,
                           CI_upper = meta$upper.random ,
                           Pval = meta$pval.random ,
                           I2 = meta$I2 ,
                           Q_Pval = meta$pval.Q )
  }else{
    
    print("not enough studies to do meta-analysis")
    meta_res <- data.frame(Gene = genes[k],
                           Coef = NA,
                           SE =  NA,
                           CI_lower = NA,
                           CI_upper = NA,
                           Pval = NA,
                           I2= NA,
                           Q_Pval = NA)

  }
  
meta_res

})

meta.cox <- do.call(rbind, meta.cox)
meta.cox <- meta.cox[!is.na(meta.cox$Coef), ]
meta.cox$padj <- p.adjust(meta.cox$Pval, method = "BH")

save(meta.cox, file = file.path(dir_out, "meta_cox_tcga_histo.rda"))
write.csv(meta.cox, file = file.path(dir_out, "meta_cox_tcga_histo.csv"), row.names = FALSE)










