# -----------------------------------------------------------
# Soft-Tissue Sarcoma Batch Correction & PCA Analysis Script
# This script corrects for potential batch effects in 
# soft-tissue sarcoma gene expression data and visualizes 
# cohort/subtype variation using PCA plots.
#
#   Data sources:
#     - GEO (clinical gene expression datasets)
#     - ORCESTRA (pharmacogenomics cell line data)
#
#   Includes:
#     * PCA plots of GEO datasets before/after ComBat correction
#     * PCA plots of PSet datasets before/after ComBat correction
#     * Batch effect correction using Empirical Bayes (ComBat)
#     * Merging corrected datasets for downstream analysis
#
#   Assumes input data (qs format) has already been processed
# -----------------------------------------------------------
#############################################################
## Load libraries
#############################################################

library(qs)
library(dplyr)
library(sva)
library(ggplot2)

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/procdata'  
dir_out <- 'data/results/batchcorrection' 

########################################################
# GEO data: assess potential batches across datasets
########################################################
#--- part I: distribution of studies
dat <- qread(file.path(dir_in, "PGx_gse_rna_sts.qs"))

dat_mat <- dat$TCGA_mat
dat_ann <- dat$TCGA_ann

group <- factor(dat_ann$Type)
cols <- c("#FED789FF", "#FED789FF", "#476F84FF", "#72874EFF")

pca_results <- prcomp(dat_mat, scale = TRUE)
var_res <- pca_results$sdev^2 / sum(pca_results$sdev^2)

pdf(file= file.path(dir_out, "pca_gse_cohort_before_correction.pdf"), width = 5, height = 5)

plot(pca_results$x[,1],pca_results$x[,2],pch="+", xlab='PC1 (40%)',ylab='PC2 (4%)')
ind= group == levels(group)[1]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[1])
ind= group == levels(group)[2]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[2])
ind= group == levels(group)[3]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[3])
ind= group == levels(group)[4]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[4])
legend('center',paste(levels(group)), pch="+",col=cols, cex=0.7)

dev.off()

# part II: distribution of sts subtypes
group <- dat_ann$subtype
group <- factor(group)

cols <- c( "#AAC197FF", "#4A7169FF",  "#9A9391FF", "#BC8E7DFF")

pdf(file= file.path(dir_out, "pca_gse_subtype_before_correction.pdf"), width = 5, height = 5)

plot(pca_results$x[,1],pca_results$x[,2],pch="+", xlab='PC1 (40%)',ylab='PC2 (4%)')
ind= group == levels(group)[1]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[1])
ind= group == levels(group)[2]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[2])
ind= group == levels(group)[3]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[3])
ind= group == levels(group)[4]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[4])
legend('center',paste(levels(group)), pch="+",col=cols, cex=0.7)

dev.off()

#################################################################
## Patient data: ComBat (Empirical Bayes)
#################################################################

dat <- qread(file.path(dir_in, "PGx_gse_rna_sts.qs"))

dat_mat <- dat$TCGA_mat
dat_mat <- t(dat_mat)
dat_ann <- dat$TCGA_ann

batch <- dat_ann$Type
batch <- ifelse(batch %in% c("cohort 1-GSE21050", "cohort 2-GSE21050"), "GSE21050", "GSE21122/GSE30929")
batch <- as.factor(batch)

combat_mat <- ComBat(dat=dat_mat, batch=batch, mod=NULL)
dat_mat <- t(combat_mat)

## regenerate figures after correction

group <- factor(dat_ann$Type)
cols <- c("#FED789FF", "#FED789FF", "#476F84FF", "#72874EFF")

pca_results <- prcomp(dat_mat, scale = TRUE)
var_res <- pca_results$sdev^2 / sum(pca_results$sdev^2)

pdf(file= file.path(dir_out, "pca_gse_cohort_after_correction.pdf"), width = 5, height = 5)

plot(pca_results$x[,1],pca_results$x[,2],pch="+", xlab='PC1 (8%)',ylab='PC2 (6%)')
ind= group == levels(group)[1]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[1])
ind= group == levels(group)[2]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[2])
ind= group == levels(group)[3]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[3])
ind= group == levels(group)[4]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[4])
legend('topright',paste(levels(group)), pch="+",col=cols, cex=0.5)

dev.off()


group <- dat_ann$subtype
group <- factor(group)

cols <- c( "#AAC197FF", "#4A7169FF",  "#9A9391FF", "#BC8E7DFF")

pdf(file= file.path(dir_out, "pca_gse_subtype_after_correction.pdf"), width = 5, height = 5)


plot(pca_results$x[,1],pca_results$x[,2],pch="+", xlab='PC1 (7%)',ylab='PC2 (6%)', 
     main='GEO RNA data')
ind= group == levels(group)[1]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[1])
ind= group == levels(group)[2]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[2])
ind= group == levels(group)[3]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[3])
ind= group == levels(group)[4]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[4])
legend('topright',paste(levels(group)), pch="+",col=cols, cex=0.5)

dev.off()

########################################################################
## Merge the data as corrected GEO
########################################################################

pset_ann <- dat$pset_ann
pset_mat <- dat$pset_mat
TCGA_ann <- dat$TCGA_ann
TCGA_mat <- dat_mat

dat <- list(TCGA_mat = TCGA_mat,
            TCGA_ann = TCGA_ann,
            pset_ann = pset_ann,
            pset_mat = pset_mat)

qsave(dat, file= file.path(dir_in, "PGx_correctedGEO_gse_rna_sts.qs"))

###############################################################################################################
###############################################################################################################
############################ Correct for potential batches in PSets ###########################################
###############################################################################################################
###############################################################################################################

dat <- qread(file.path(dir_in, "PGx_correctedGEO_gse_rna_sts.qs"))

dat_mat <- dat$pset_mat
dat_ann <- dat$pset_ann

group <- dat_ann$Type
group <- factor(group)

cols <- c("#1b7837", "#542788", "#bf812d")

pca_results <- prcomp(dat_mat, scale = TRUE)
var_res <- pca_results$sdev^2 / sum(pca_results$sdev^2)

pdf(file= file.path(dir_out, "pca_pset_study_before_correction.pdf"), width = 5, height = 5)

plot(pca_results$x[,1],pca_results$x[,2],pch="+", xlab='PC1 (52%)',ylab='PC2 (19%)')
ind= group == levels(group)[1]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[1])
ind= group == levels(group)[2]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[2])
ind= group == levels(group)[3]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[3])
legend('topright',paste(levels(group)), pch="+",col=cols, cex=0.7)

dev.off()

#################################################################
## PSets: ComBat (Empirical Bayes)
#################################################################

batch <- dat_ann$Type
batch <- as.factor(batch)

combat_mat <- ComBat(dat=t(dat_mat), batch=batch, mod=NULL)
dat_mat <- t(combat_mat)

## regenerate figures after correction
group <- factor(dat_ann$Type)
cols <- c("#1b7837", "#542788", "#bf812d")

pca_results <- prcomp(dat_mat, scale = TRUE)
var_res <- pca_results$sdev^2 / sum(pca_results$sdev^2)

pdf(file= file.path(dir_out, "pca_pset_study_after_correction.pdf"), width = 5, height = 5)

plot(pca_results$x[,1],pca_results$x[,2],pch="+", xlab='PC1 (12%)',ylab='PC2 (10%)')
ind= group == levels(group)[1]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[1])
ind= group == levels(group)[2]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[2])
ind= group == levels(group)[3]
points(pca_results$x[ind,1],pca_results$x[ind,2], pch="+",col=cols[3])
legend('topleft',paste(levels(group)), pch="+",col=cols, cex=0.5)

dev.off()

########################################################################
## Merge the data as corrected
########################################################################

pset_ann <- dat$pset_ann
pset_mat <- dat_mat
TCGA_ann <- dat$TCGA_ann
TCGA_mat <- dat$TCGA_mat

dat <- list(TCGA_mat = TCGA_mat,
            TCGA_ann = TCGA_ann,
            pset_ann = pset_ann,
            pset_mat = pset_mat)

qsave(dat, file= file.path(dir_in, "PGx_corrected_gse_rna_sts.qs"))

# Note
# 1. correct for GEO data and remove potential variation between cohorts: PGx_correctedGEO_gse_rna_sts
# 2. Correct for PSets and remove potential variation between cohorts
# 3. Merge them as a corrected data: PGx_corrected_gse_rna_sts
