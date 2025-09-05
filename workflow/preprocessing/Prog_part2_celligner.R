##-------------------------------------------------------------------------------
## Celligner Alignment & Visualization Script (Part 1)
## This script performs optimization of the number of contrastive principal 
## components (cPCs) to be used in the Celligner alignment between tumor and 
## cell line expression datasets.
##
## It includes:
##   - Loading required functions for Celligner integration
##   - Reading preprocessed expression data and annotations
##   - Calculating gene-wise expression statistics
##   - Performing clustering and differential expression analysis
##   - Running cPCA to identify contrastive sources of variation
##   - Iteratively removing top k cPCs (1–10) and running MNN correction
##   - Evaluating number of MNN pairs as a function of cPCs removed
##   - Visualizing MNN pairing trend and cPCA variance contributions
##
## Input: 
##   - PGx_corrected_gse_rna_sts.qs (processed expression dataset)
##   - hgnc_complete_set.txt (gene annotations)
##
## Output:
##   - pgx_rna_celligner_sts_pcs.qs (alignment results per cPC)
##   - optimal_rna_sts_k.pdf (plot of MNN count vs cPC count)
##   - cpca_sts_variance.pdf (barplot of explained variance by cPCs)
##
## Dependencies: Seurat, limma, batchelor, qread, ggplot2, etc.
## Note: Update file paths accordingly when running on different environments
##-------------------------------------------------------------------------------
##################################################################
## Libraries and functions
##################################################################
source("scripts/Prog_celligner_function.R")

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/procdata'  
dir_out <- 'data/results/aligned/'  
dir_annot <-  'data/rawdata' 

#################################################################
## Run Celligner
#################################################################
# Download gene annotation (ensemble ids) 
hgnc.complete.set <- data.table::fread( file.path(dir_annot, "hgnc_complete_set.txt")) %>% as.data.frame() 

# load corrected cell line and tumor patient data
pgx_rna <- run_Celligner(data_dir =  dir_annot, 
                         dat =  file.path(dir_in, "PGx_corrected_gse_rna_sts.qs"), 
                         remove_cPCA_dims = 1:7, 
                         plotname = " ")

qsave(pgx_rna, file = file.path(dir_out, "pgx_rna_celligner_sts.qs"))

################################################################################
## Visualization Part 1: before alignment
################################################################################

plot_uncorrected_data_type(org_dat = qread(file.path(dir_in, "PGx_corrected_gse_rna_sts.qs")), 
                           before_plot = "pgx_tumor_uncorrected",
                           dir_output = dir_out)

plot_uncorrected_data_study(org_dat = qread(file.path(dir_in, "PGx_corrected_gse_rna_sts.qs")), 
                            before_plot = "pgx_tumor_uncorrected_study",
                            dir_output = dir_out)

################################################################################
## Visualization Part 2: after alignment
################################################################################

dat <- qread(file.path(dir_out, "pgx_rna_celligner_sts.qs"))

Celligner_alignment_plot_type(aligned_data = dat$comb_obj, 
                              after_plot = "pgx_tumor_corrected",
                              dir_output = dir_out)

Celligner_alignment_plot_study(aligned_data = dat$comb_obj, 
                               org_dat = qread(file.path(dir_in, "PGx_corrected_gse_rna_sts.qs")),
                               after_plot =  "pgx_tumor_corrected_study",
                               dir_output = dir_out)



