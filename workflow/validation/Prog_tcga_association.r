#--------------------------------------------------------------------------------
## Description:
#   This script performs gene-level survival analysis on TCGA-SARC RNA-seq data,
#   using a **late integration (LI)** strategy to account for histological subtype heterogeneity.
#
#   Strategy:
#   (1) Late integration (LI): Per-histology Cox models are first computed independently.
#       These are then integrated using a random-effects meta-analysis (DerSimonian–Laird
#       with Hartung–Knapp adjustment) to obtain a pooled log-HR and assess heterogeneity.
#
# Required input:
#   - SummarizedExperiment (qs): data/procdata/validation/curated_TCGA_se.qs
#       * assay(dat): gene expression matrix (genes × samples)
#       * colData(dat): clinical data with os_time (months), os_event (0/1),
#                       histological subtype labels (histological), and sample ID (patientid)
#       * rowData(dat): gene annotations (row names must match expression genes)
#
# Output files (written to data/results/validation/TCGA/):
#   - cox_tcga_histo.rda / cox_tcga_histo.csv
#       * Per-histology (subtype-specific) Cox model results per gene.
#       * Columns: histo, gene_name, HR (log-HR), SE, N, Low/Up (HR scale),
#                  pval, padj (BH-corrected).
#
#   - meta_cox_tcga_histo.rda / meta_cox_tcga_histo.csv
#       * Late-integration meta-analysis across histologies.
#       * Columns: Gene, Coef (pooled log-HR), SE, CI_lower/CI_upper (log-HR scale),
#                  Pval, padj (BH), I2 (heterogeneity), Q_Pval (heterogeneity p).
#       * Meta-analysis is only performed for genes with results in ≥3 histologies.
#
# Key analysis settings:
#   - Remove low (or zero) expressed genes (zero expressed across 70% of samples)
#   - Administrative censoring at 36 months (time.censor = 36).
#       * Survival times > 36 months are right-censored at 36.
#   - Model form: Surv(time, status) ~ expression
#       * Per-gene univariate Cox models (no covariates other than gene expression).
#
# Dependencies:
#   - qs, survival, data.table, magrittr, SummarizedExperiment, meta
#
# Notes:
#   - "HR" refers to log hazard ratio (i.e., Cox model coefficient).
#   - CI bounds (Low/Up) are exponentiated to reflect HR scale.
#   - Meta-analysis outputs (Coef, CI) are in log-HR scale; exponentiate if needed.
#--------------------------------------------------------------------------------
########################################################################
## Load library
########################################################################

library(qs)
library(survival)
library(data.table)
library(magrittr)
library(meta)
library(reshape2)
library(SummarizedExperiment)

########################################################################
## Function to Remove low-expressed genes
########################################################################

 rem <- function (x, missing.perc, const.int){
    x <- as.matrix(x)
    x <- t(apply(x, 1, as.numeric))
    r <- as.numeric(apply(x, 1, function(i) sum(round(i, 6) == 
        round(log2(const.int), 6))))
    remove <- which(r > dim(x)[2] * missing.perc)
    return(remove)
}

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

expr <- assay(dat) # 16522 genes
clin <- data.frame(colData(dat))
annot <- data.frame(rowData(dat))

########################################################################
## Fit proportional hazard model 
########################################################################
## Step 1: Find association per subtype
group <- unique(clin$histological)

res.cox <- lapply(1:length(group), function(j){

print(j)

df.clin <- clin[clin$histological == group[j], ]
df.expr <- expr[, colnames(expr) %in% df.clin$patientid]

# filter out low-zero expressed genes
remove <- rem(df.expr, 0.70, 1)
df.expr <- df.expr[-remove, ] 

res <- lapply(1:nrow(df.expr), function(k){

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
           gene_name = rownames(df.expr)[k],
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

########################################################################
## Integration analysis 
########################################################################
## Step 2: Integrate across subtypes
load(file.path(dir_out, 'cox_tcga_histo.rda'))
genes <- unique(res.cox$gene_name) # 15481 genes

meta.cox <- lapply(1:length(genes), function(k){

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










