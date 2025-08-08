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
source("scripts/Prog_preprocess.R")

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

## meta association result
dat <- qread(file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
dat_drug <- dat_drug[rownames(dat_drug) %in% unique(dat$drug), ]

## get the class of drugs
dat_drug$TARGET_PATHWAY <- ifelse(dat_drug$TARGET_PATHWAY %in% c("Other, kinases", "Other"), 
                                  "Other", dat_drug$TARGET_PATHWAY)

dat_drug$TARGET_PATHWAY <- ifelse(dat_drug$TARGET_PATHWAY %in% c("PI3K/MTOR signaling"), 
                                  "PI3K MTOR signaling", dat_drug$TARGET_PATHWAY)
dat_drug$TARGET_PATHWAY <- ifelse(dat_drug$TARGET_PATHWAY %in% c("RTK signaling///IGF1R signaling"), 
                                  "IGF1R signaling", dat_drug$TARGET_PATHWAY)

target_pathway <- unique(dat_drug$TARGET_PATHWAY)
freq_target_pathway <- lapply(1:length(target_pathway), function(k){
  
 data.frame(target_pathway = target_pathway[k], 
            freq = nrow( dat_drug[dat_drug$TARGET_PATHWAY == target_pathway[k], ] ))
  
})

freq_target_pathway  <- do.call(rbind, freq_target_pathway)
freq_target_pathway <- freq_target_pathway[order(freq_target_pathway$freq, decreasing = TRUE), ]
#freq_target_pathway <- freq_target_pathway[freq_target_pathway$target_pathway != "Other", ]

freq_target_pathway_included <- freq_target_pathway[freq_target_pathway$freq >= 2, ]
target_pathway_included <- freq_target_pathway_included$target_pathway

dat_drug_class <- lapply(1:length(target_pathway_included), function(k){
  
  dat_drug[dat_drug$TARGET_PATHWAY == target_pathway_included[k], "treatmentid"]
  
})

names(dat_drug_class) <- target_pathway_included

################################################################################
## Pearson correlation as similarity metric  
################################################################################

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


  pdf(file=file.path(dir_out, 'corPlot', paste(names(dat_drug_class)[k], ".pdf", sep="")),
       width = 6, height = 6)
  
  p <- corrplot(res, type = "upper", order = "hclust", 
           tl.col = "black", tl.srt = 45, tl.cex = 0.8)
  print(p)
  
  dev.off()
  
}

##########################################################
## Drugs with significant meta-association results
##########################################################
# Consider cut-offs for correlation (no cut-off, 30%, and 50%)
sig <- dat[dat$padj < 0.05 & abs(dat$r) >= 0.5, ]
dat_drug <- dat_drug[rownames(dat_drug) %in% unique(sig$drug), ]

drug <- unique(sig$drug)
res <- lapply(1:length(drug), function(k){
  
  df <- dat[dat$drug == drug[k], "r"]
  df
  
})

res <- do.call(cbind, res)
colnames(res) <- drug
cor_res <- cor(res)

pdf(file=file.path(dir_out, paste("cor_pcl_50PercCutoffCor", ".pdf", sep="")),
     width = 12, height = 12)

p <- corrplot(cor_res, type = "upper", order = "hclust", 
              tl.col = "black", tl.srt = 90, tl.cex = 0.8)
p

dev.off()


