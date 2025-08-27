# -----------------------------------------------------------
# Pharmacological Class-Level (PCL) Similarity Analysis Script
#
# This script performs correlation-based similarity analysis 
# of drug response–associated gene signatures at the 
# pharmacological class level (PCL), using meta-association 
# results across soft-tissue sarcoma cell lines.
#
#   Includes:
#     * Mapping drugs to pharmacological classes (target pathways)
#     * Filtering for classes with ≥2 drugs
#     * Computing within-class Pearson correlations across gene–drug associations
#     * Visualizing PCL similarity using:
#         - Boxplots of correlation distributions per class
#         - Correlation matrices via corrplot
#
#   Inputs:
#     - Meta-analysis results of gene–drug correlations (FDR-adjusted)
#     - Drug annotation and pathway information from GDSCv2
#
#   Outputs:
#     - Boxplot of within-class predictive correlations
#     - Correlation heatmaps for each drug class and overall
#     - Flattened correlation tables for downstream reporting
#
#   Assumes meta-analysis results and drug annotations are preprocessed
#   and available in `.qs` format under defined input directories.
# -----------------------------------------------------------
#############################################################
## Load libraries
#############################################################

library(qs)
library(PharmacoGx)
library(data.table)
library(magrittr)
library(corrplot)
library(ggplot2)
library(Hmisc)
library(paletteer)

##################################################################
## Define function
##################################################################

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

##################################################################
## Setup directory
##################################################################
source("scripts/Prog_data_function.R")

dir_data <- 'data/rawdata'   
dir_in <- 'data/results/drug'  
dir_out <- 'data/results/MOA' 

#########################################################################################################################
#########################################################################################################################
########################################## Load Gene Drug response (AAC) Association ####################################
#########################################################################################################################
#########################################################################################################################
dat <- qread(file = file.path(dir_data, "sarcsets.qs"))
names(dat) <- c("CCLE", "CTRP", "gCSI",  "GDSCv2", "GDSCv1", "NCI60", "PRISM", "NCI_Sarcoma")

for(i in 1:length(dat)){
  dat[[i]] <- updateObject(dat[[i]])
}

## class of drugs
dat_drug  <- drugInfo(dat$GDSCv2)
dat_drug_prism <- drugInfo(dat$PRISM)

## meta association result
dat <- qread(file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
dat_drug <- dat_drug[rownames(dat_drug) %in% unique(dat$drug), ]
dat_drug_prism <- dat_drug_prism[rownames(dat_drug_prism) %in% unique(dat$drug), ]

## get the class of drugs
moa_dat <- lapply(1:nrow(dat_drug), function(k){

  sub_dat_drug <- dat_drug[rownames(dat_drug) == rownames(dat_drug)[k], ]
  sub_dat_drug_prism <- dat_drug_prism[rownames(dat_drug_prism) == rownames(dat_drug)[k], ]

  if(nrow(sub_dat_drug_prism) > 0){

     data.frame(treatmentid = rownames(dat_drug)[k],
             synonym = sub_dat_drug$SYNONYMS, 
             drug_target = sub_dat_drug$TARGET,
             TARGET_PATHWAY = sub_dat_drug$TARGET_PATHWAY,
             inchikey = sub_dat_drug$inchikey,
             target = sub_dat_drug_prism$target,
             moa = sub_dat_drug_prism$moa,
             FDA =sub_dat_drug$FDA )
  }else{

   data.frame(treatmentid = rownames(dat_drug)[k],
             synonym = sub_dat_drug$SYNONYMS, 
             drug_target = sub_dat_drug$TARGET,
             TARGET_PATHWAY = sub_dat_drug$TARGET_PATHWAY,
             inchikey = sub_dat_drug$inchikey,
             target = NA,
             moa = NA,
             FDA =sub_dat_drug$FDA )

  }
   
})

moa_dat <- do.call(rbind, moa_dat)
moa_dat$TARGET_PATHWAY_updated <- ifelse(moa_dat$TARGET_PATHWAY %in% c("Other, kinases", "Other"), 
                                  "Other", moa_dat$TARGET_PATHWAY)
moa_dat$TARGET_PATHWAY_updated <- ifelse(moa_dat$TARGET_PATHWAY_updated %in% c("PI3K/MTOR signaling"), 
                                  "PI3K MTOR signaling", moa_dat$TARGET_PATHWAY_updated)
moa_dat$TARGET_PATHWAY_updated <- ifelse(moa_dat$TARGET_PATHWAY_updated %in% c("RTK signaling///IGF1R signaling"), 
                                  "IGF1R signaling", moa_dat$TARGET_PATHWAY_updated)


target_pathway <- unique(moa_dat$TARGET_PATHWAY_updated)
freq_target_pathway <- lapply(1:length(target_pathway), function(k){
  
 data.frame(target_pathway = target_pathway[k], 
            freq = nrow( moa_dat[moa_dat$TARGET_PATHWAY_updated == target_pathway[k], ] ))
  
})

freq_target_pathway  <- do.call(rbind, freq_target_pathway)
freq_target_pathway <- freq_target_pathway[order(freq_target_pathway$freq, decreasing = TRUE), ]
#freq_target_pathway <- freq_target_pathway[freq_target_pathway$target_pathway != "Other", ]

freq_target_pathway_included <- freq_target_pathway[freq_target_pathway$freq >= 2, ]
target_pathway_included <- freq_target_pathway_included$target_pathway

dat_drug_class <- lapply(1:length(target_pathway_included), function(k){
  
  moa_dat[moa_dat$TARGET_PATHWAY_updated == target_pathway_included[k], "treatmentid"]
  
})

names(dat_drug_class) <- target_pathway_included
save(dat_drug_class, file = file.path(dir_out, 'moa_info.rda'))
write.csv(moa_dat, file = file.path(dir_out, 'moa_dat.csv'), row.names = FALSE)

################################################################################
## Pearson correlation as similarity metric (69 drugs) 
################################################################################
load(file.path(dir_out, 'moa_info.rda'))
cor_res <- lapply(1:length(dat_drug_class), function(k){
  
  sub_res <- dat[dat$drug %in% dat_drug_class[[k]], ]
  
  sub_res_matrix <- lapply(1:length(dat_drug_class[[k]]), function(j){
    
    sub_res[sub_res$drug == dat_drug_class[[k]][j], "r"]
    
  })
  
  sub_res_matrix <- do.call(cbind, sub_res_matrix)
  colnames(sub_res_matrix) <- dat_drug_class[[k]]
  res <- rcorr(as.matrix(sub_res_matrix))
  data.frame(target_pathway = names(dat_drug_class)[k],
             flattenCorrMatrix(res$r, res$P) )
  
  
})

cor_res <- do.call(rbind, cor_res)
write.csv(cor_res, file=file.path(dir_out, "cor_pcl.csv"), row.names = FALSE)  

################################################################################
## Boxplot for all the pharmacological classes
################################################################################

pdf(file= file.path(dir_out, "boxplot_pcl.pdf"), 
     width = 4, height = 4)

p <- ggplot(cor_res, aes(x=reorder(target_pathway, -cor), y=cor)) + 
  geom_boxplot(width = 0.5, fill= "#D8D8D8FF", color = "#181830FF") +
  coord_flip() + 
  ylim(c(-0.5, 1)) +
  xlab("") +
  ylab("predictive correlation") +
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size=9),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.text = element_text(size = 6, face="bold"),
        legend.title = element_blank()) 

p

dev.off()

################################################################################
## Upper triangles
################################################################################

for(k in 1:length(dat_drug_class)){
  
  print(k)
  sub_res <- dat[dat$drug %in% dat_drug_class[[k]], ]
  sub_res_matrix <- lapply(1:length(dat_drug_class[[k]]), function(j){
    
    sub_res[sub_res$drug == dat_drug_class[[k]][j], "r"]
    
  })
  
  sub_res_matrix <- do.call(cbind, sub_res_matrix)
  rem <- which(dat_drug_class[[k]] == "Unii-40E3azg1MX")
  if(length(rem) > 0){
    
    dat_drug_class[[k]][rem] <- "BMS-536924"
    
  }
  colnames(sub_res_matrix) <- dat_drug_class[[k]]
 
  res <- cor(sub_res_matrix)


  pdf(file=file.path(dir_out, 'corPlot/all', paste(names(dat_drug_class)[k], ".pdf", sep="")),
       width = 6, height = 6)
  
  p <- corrplot(res, type = "upper", order = "hclust", 
           tl.col = "black", tl.srt = 45, tl.cex = 0.8)
  print(p)
  
  dev.off()
  
}

##########################################################
## Drugs with significant meta-association results 
##########################################################
# Consider cut-offs for correlation (no cut-off, 30%)
dat <- qread(file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
sig <- dat[dat$padj < 0.05 & abs(dat$r) >= 0.3, ]

drug <- unique(sig$drug)
res <- lapply(1:length(drug), function(k){
  
  df <- dat[dat$drug == drug[k], "r"]
  df
  
})

res <- do.call(cbind, res)
colnames(res) <- drug
cor_res <- cor(res)

pdf(file=file.path(dir_out, paste("cor_pcl_30PercCutoffCor", ".pdf", sep="")),
     width = 10, height = 10)

p <- corrplot(cor_res, type = "upper", order = "hclust", 
              tl.col = "black", tl.srt = 90, tl.cex = 0.8)
p

dev.off()

#####################################################################################################################
################################################ clincial drugs #####################################################
#####################################################################################################################
## load selected drugs
selected_drugs <- read.csv(file.path(dir_in, 'selected_drugs.csv'))
clin_selected_drugs <- selected_drugs[selected_drugs$sts == 'Yes', 'drug']

dat_drug <- dat_drug[rownames(dat_drug) %in% clin_selected_drugs, ]
## get the class of drugs
moa_dat <- lapply(1:nrow(dat_drug), function(k){

  sub_dat_drug <- dat_drug[rownames(dat_drug) == rownames(dat_drug)[k], ]
  sub_dat_drug_prism <- dat_drug_prism[rownames(dat_drug_prism) == rownames(dat_drug)[k], ]

  if(nrow(sub_dat_drug_prism) > 0){

     data.frame(treatmentid = rownames(dat_drug)[k],
             synonym = sub_dat_drug$SYNONYMS, 
             drug_target = sub_dat_drug$TARGET,
             TARGET_PATHWAY = sub_dat_drug$TARGET_PATHWAY,
             inchikey = sub_dat_drug$inchikey,
             target = sub_dat_drug_prism$target,
             moa = sub_dat_drug_prism$moa,
             FDA =sub_dat_drug$FDA )
  }else{

   data.frame(treatmentid = rownames(dat_drug)[k],
             synonym = sub_dat_drug$SYNONYMS, 
             drug_target = sub_dat_drug$TARGET,
             TARGET_PATHWAY = sub_dat_drug$TARGET_PATHWAY,
             inchikey = sub_dat_drug$inchikey,
             target = NA,
             moa = NA,
             FDA =sub_dat_drug$FDA )

  }
   
})

moa_dat <- do.call(rbind, moa_dat)
moa_dat$TARGET_PATHWAY_updated <- ifelse(moa_dat$TARGET_PATHWAY %in% c("Other, kinases", "Other"), 
                                  "Other", moa_dat$TARGET_PATHWAY)
moa_dat$TARGET_PATHWAY_updated <- ifelse(moa_dat$TARGET_PATHWAY_updated %in% c("PI3K/MTOR signaling"), 
                                  "PI3K MTOR signaling", moa_dat$TARGET_PATHWAY_updated)
moa_dat$TARGET_PATHWAY_updated <- ifelse(moa_dat$TARGET_PATHWAY_updated %in% c("RTK signaling///IGF1R signaling"), 
                                  "IGF1R signaling", moa_dat$TARGET_PATHWAY_updated)

target_pathway <- unique(moa_dat$TARGET_PATHWAY_updated)
freq_target_pathway <- lapply(1:length(target_pathway), function(k){
  
 data.frame(target_pathway = target_pathway[k], 
            freq = nrow( moa_dat[moa_dat$TARGET_PATHWAY_updated == target_pathway[k], ] ))
  
})

freq_target_pathway  <- do.call(rbind, freq_target_pathway)
freq_target_pathway <- freq_target_pathway[order(freq_target_pathway$freq, decreasing = TRUE), ]
#freq_target_pathway <- freq_target_pathway[freq_target_pathway$target_pathway != "Other", ]

freq_target_pathway_included <- freq_target_pathway[freq_target_pathway$freq >= 2, ]
target_pathway_included <- freq_target_pathway_included$target_pathway

dat_drug_class <- lapply(1:length(target_pathway_included), function(k){
  
  moa_dat[moa_dat$TARGET_PATHWAY_updated == target_pathway_included[k], "treatmentid"]
  
})

names(dat_drug_class) <- target_pathway_included
save(dat_drug_class, file = file.path(dir_out, 'moa_info_clinical.rda'))

################################################################################
## Pearson correlation as similarity metric (69 drugs) 
################################################################################

load(file.path(dir_out, 'moa_info_clinical.rda'))
cor_res <- lapply(1:length(dat_drug_class), function(k){
  
  sub_res <- dat[dat$drug %in% dat_drug_class[[k]], ]
  
  sub_res_matrix <- lapply(1:length(dat_drug_class[[k]]), function(j){
    
    sub_res[sub_res$drug == dat_drug_class[[k]][j], "r"]
    
  })
  
  sub_res_matrix <- do.call(cbind, sub_res_matrix)
  colnames(sub_res_matrix) <- dat_drug_class[[k]]
  res <- rcorr(as.matrix(sub_res_matrix))
  data.frame(target_pathway = names(dat_drug_class)[k],
             flattenCorrMatrix(res$r, res$P) )
  
  
})

cor_res <- do.call(rbind, cor_res)
write.csv(cor_res, file=file.path(dir_out, "cor_pcl_clinical.csv"), row.names = FALSE)  

################################################################################
## Boxplot for all the pharmacological classes
################################################################################

pdf(file= file.path(dir_out, "boxplot_pcl_clinical.pdf"), 
     width = 3, height = 3)

p <- ggplot(cor_res, aes(x=reorder(target_pathway, -cor), y=cor)) + 
  geom_boxplot(width = 0.5, fill= "#D8D8D8FF", color = "#181830FF") +
  coord_flip() + 
  ylim(c(-0.5, 1)) +
  xlab("") +
  ylab("predictive correlation") +
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size=9),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.text = element_text(size = 6, face="bold"),
        legend.title = element_blank()) 

p

dev.off()

################################################################################
## Upper triangles
################################################################################

for(k in 1:length(dat_drug_class)){
  
  print(k)
  sub_res <- dat[dat$drug %in% dat_drug_class[[k]], ]
  sub_res_matrix <- lapply(1:length(dat_drug_class[[k]]), function(j){
    
    sub_res[sub_res$drug == dat_drug_class[[k]][j], "r"]
    
  })
  
  sub_res_matrix <- do.call(cbind, sub_res_matrix)
  rem <- which(dat_drug_class[[k]] == "Unii-40E3azg1MX")
  if(length(rem) > 0){
    
    dat_drug_class[[k]][rem] <- "BMS-536924"
    
  }
  colnames(sub_res_matrix) <- dat_drug_class[[k]]
 
  res <- cor(sub_res_matrix)


  pdf(file=file.path(dir_out, 'corPlot/clinical', paste(names(dat_drug_class)[k], ".pdf", sep="")),
       width = 6, height = 6)
  
  p <- corrplot(res, type = "upper", order = "hclust", 
           tl.col = "black", tl.srt = 45, tl.cex = 0.8)
  print(p)
  
  dev.off()
  
}

##########################################################
## Drugs with significant meta-association results 
##########################################################
# Consider cut-offs for correlation (no cut-off, 30%)
dat <- qread(file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
sig <- dat[dat$padj < 0.05 & abs(dat$r) >= 0.3, ]

drug <- intersect(unique(sig$drug), clin_selected_drugs)

res <- lapply(1:length(drug), function(k){
  
  df <- dat[dat$drug == drug[k], "r"]
  df
  
})

res <- do.call(cbind, res)
colnames(res) <- drug
cor_res <- cor(res)

pdf(file=file.path(dir_out, paste("cor_pcl_30PercCutoffCor_clinical", ".pdf", sep="")),
     width = 5, height = 5)

p <- corrplot(cor_res, type = "upper", order = "hclust", 
              tl.col = "black", tl.srt = 90, tl.cex = 0.8)
p

dev.off()

