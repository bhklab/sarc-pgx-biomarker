# ------------------------------------------------------------------------------
# STS RNA-seq Alignment Evaluation Script
# This script evaluates the effect of data alignment (via Celligner) on gene 
# expression similarity between tumor (TCGA) and cell line (PSet) samples.
#
# The following analyses are included:
#   * Anderson-Darling and Kolmogorov-Smirnov tests on gene distributions
#     before and after alignment
#   * Visualization of the proportion of genes with significantly shifted 
#     distributions (AD test)
#   * Correlation analysis between tumors and cell lines (aligned vs raw)
#   * Distribution plots and Wilcoxon test to assess alignment impact
#
# Assumes the aligned and unaligned datasets are already preprocessed.
# ------------------------------------------------------------------------------
########################################
## Load libraries and functions
########################################

library(twosamples)
library(ggplot2)
source("scripts/Prog_celligner_function.R")

##################################################################
## Setup directory
##################################################################

dir_data <- 'data/procdata' 
dir_aligned <- 'data/results/aligned' 
dir_out <- 'data/results/aligned'  

##########################################
## STS before alignment (AD: scaled data)
##########################################

dat_sts <- qread(file.path(dir_data, "PGx_corrected_gse_rna_sts.qs"))

dat_pset <- dat_sts$pset_mat
dat_tcga <- dat_sts$TCGA_mat

int <- colnames(dat_pset)

tcga_df <- data.frame(group = "TCGA", dat_tcga)
pset_df <- data.frame(group = "PSet", dat_pset)
gene_id <- colnames(pset_df)[-1]

ad_ks_res <- lapply(1:length(gene_id), function(k){
  
  print(k)  
  dat <- rbind(tcga_df[, c("group", gene_id[k])], 
               pset_df[, c("group", gene_id[k])])
  
  colnames(dat)[2] <- "gene"
  
  fit_ad <- ad_test(dat[dat$group == "PSet", "gene"], 
                    dat[dat$group == "TCGA", "gene"],
                    nboots = 1000)
  
  fit_ks <- ks_test(dat[dat$group == "PSet", "gene"], 
                    dat[dat$group == "TCGA", "gene"],
                    nboots = 1000)
  
  data.frame(gene = gene_id[k],
             ad_effect = fit_ad[1],
             ad_pvalue = fit_ad[2],
             ks_effect = fit_ks[1],
             ks_pvalue = fit_ks[2])
})

ad_ks_res <- do.call(rbind, ad_ks_res)
ad_ks_res$ad_padj_BH <- p.adjust(ad_ks_res$ad_pvalue, method = "BH")
ad_ks_res$ad_padj_bonferroni <- p.adjust(ad_ks_res$ad_pvalue, method = "bonferroni")
ad_ks_res$ks_padj_BH <- p.adjust(ad_ks_res$ks_pvalue, method = "BH")
ad_ks_res$ks_padj_bonferroni <- p.adjust(ad_ks_res$ks_pvalue, method = "bonferroni")
rownames(ad_ks_res) <- NULL

qsave(ad_ks_res, file = file.path(dir_out, "pgx_rna_before_alignment_ad_ks.qs"))

##########################################
## STS after alignment (AD: scaled data)
##########################################

dat_sts <- qread( file.path(dir_aligned, "pgx_rna_celligner_sts.qs"))

dat_expr_sts <- t(GetAssayData(dat_sts$comb_obj))
dat_annot_sts <- data.frame(sampleID =  dat_sts$comb_obj$sampleID,
                            type = dat_sts$comb_obj$type)

dat_pset <- dat_expr_sts[rownames(dat_expr_sts) %in% dat_annot_sts[dat_annot_sts$type == "CL", "sampleID"], ]
dat_tcga <- dat_expr_sts[rownames(dat_expr_sts) %in% dat_annot_sts[dat_annot_sts$type == "tumor", "sampleID"], ]

int <- colnames(dat_pset)

tcga_df <- data.frame(group = "TCGA", dat_tcga)
pset_df <- data.frame(group = "PSet", dat_pset)
gene_id <- colnames(pset_df)[-1]
  
ad_ks_res <- lapply(1:length(gene_id), function(k){
    
    print(k)  
    dat <- rbind(tcga_df[, c("group", gene_id[k])], 
                 pset_df[, c("group", gene_id[k])])
    
    colnames(dat)[2] <- "gene"
    
    fit_ad <- ad_test(dat[dat$group == "PSet", "gene"], 
                      dat[dat$group == "TCGA", "gene"],
                      nboots = 1000)
    
    fit_ks <- ks_test(dat[dat$group == "PSet", "gene"], 
                      dat[dat$group == "TCGA", "gene"],
                      nboots = 1000)
    
    data.frame(gene = gene_id[k],
               ad_effect = fit_ad[1],
               ad_pvalue = fit_ad[2],
               ks_effect = fit_ks[1],
               ks_pvalue = fit_ks[2])
  })
  
ad_ks_res <- do.call(rbind, ad_ks_res)
ad_ks_res$ad_padj_BH <- p.adjust(ad_ks_res$ad_pvalue, method = "BH")
ad_ks_res$ad_padj_bonferroni <- p.adjust(ad_ks_res$ad_pvalue, method = "bonferroni")
ad_ks_res$ks_padj_BH <- p.adjust(ad_ks_res$ks_pvalue, method = "BH")
ad_ks_res$ks_padj_bonferroni <- p.adjust(ad_ks_res$ks_pvalue, method = "bonferroni")
rownames(ad_ks_res) <- NULL

qsave(ad_ks_res, file = file.path(dir_out, "pgx_rna_after_alignment_ad_ks.qs"))

################################################################################
## Compare findings (before and after alignment using celligner)
################################################################################

ad_ks_res_raw <- qread(file.path(dir_out, "pgx_rna_before_alignment_ad_ks.qs"))
ad_ks_res_raw <- data.frame(group = "uncorrected", ad_ks_res_raw)

df_raw <- data.frame(method = "uncorrected",
                     group = "uncorrected",
                     perc = table(ad_ks_res_raw$ad_padj_BH < 0.05)["TRUE"]/nrow(ad_ks_res_raw) *100)

ad_ks_res_corrected <- qread(file.path(dir_out, "pgx_rna_after_alignment_ad_ks.qs"))
ad_ks_res_corrected <- data.frame(group = "corrected", ad_ks_res_corrected )

df_celligner <- data.frame(method = "Celligner",
                           group = "aligned",
                           perc = table(ad_ks_res_corrected$ad_padj_BH < 0.05)["TRUE"]/nrow(ad_ks_res_corrected) *100 )

# merge results
df <- rbind(df_raw, df_celligner)

pdf(file= file.path(dir_out, "ad_test_celligner.pdf"), width = 5, height = 4)

p <- ggplot(df, aes(x = group, y = perc, fill = group)) +
  geom_bar(width = 0.3, stat = "identity", col = c("#b2afa8", "#f1c769")) +
  scale_fill_manual(values = c("#f1c769", "#b2afa8")) +
  coord_flip()+
  ylim(c(0, 100)) + 
  ylab("percentage of non-preserved genes") +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.y=element_text(size=10),
        strip.text = element_text(size=10),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.text = element_text(size = 6, face="bold"),
        legend.title = element_blank()
      )

p

dev.off()

#####################################################
#####################################################
## correlation analysis  
#####################################################
#####################################################
# load data
raw_dat_sts <- qread(file.path(dir_data, "PGx_corrected_gse_rna_sts.qs"))
aligned_dat_sts <- qread(file.path(dir_aligned, "pgx_rna_celligner_sts.qs"))

TCGA_mat <- raw_dat_sts$TCGA_mat
pset_mat <- raw_dat_sts$pset_mat
alignment <-  aligned_dat_sts$comb_obj

# plot aligned data

cor_res <- calc_corrected_tumor_CL_correlation(alignment)
cell_line_tumor_distance_distribution(alignment, cor_res, "correlation_aligned", dir_out)

# plot raw (or uncorrected) data

cor_res <- calc_uncorrected_tumor_CL_correlation(TCGA_mat, pset_mat)
plot_uncorrected_distribution_of_CL_tumor_distances(cor_res, alignment, "correlation_raw", dir_out)

#########################################################
## Merge corrected and uncorrected correlation for STS
#########################################################
raw_dat_sts <- qread(file.path(dir_data, "PGx_corrected_gse_rna_sts.qs"))
aligned_dat_sts <- qread(file.path(dir_aligned, "pgx_rna_celligner_sts.qs"))

TCGA_mat <- raw_dat_sts$TCGA_mat
pset_mat <- raw_dat_sts$pset_mat
alignment <-  aligned_dat_sts$comb_obj
cor_res_aligned <- calc_corrected_tumor_CL_correlation(alignment)
cor_res_raw <- calc_uncorrected_tumor_CL_correlation(TCGA_mat, pset_mat)

tumor_names <- unique(rownames(cor_res_raw))
CL_names <- unique(colnames(cor_res_raw))
dat_raw <- lapply(1:length(tumor_names), function(k){
  
  df <- cor_res_raw[k, ]
  data.frame(tumor_names = tumor_names[k], 
             CL_names = names(df),
             r = as.numeric(df))
})

dat_raw <- do.call(rbind, dat_raw)
dat_raw$group <- "uncorrected"


dat_aligned <- lapply(1:length(tumor_names), function(k){
  
  df <- cor_res_aligned[k, ]
  data.frame(tumor_names = tumor_names[k], 
             CL_names = names(df),
             r = as.numeric(df))
})

dat_aligned <- do.call(rbind, dat_aligned)
dat_aligned$group <- "aligned"
dat <- rbind(dat_raw, dat_aligned)

ggplot(dat, aes(x = r, y = group, fill = group)) +
  ggridges::geom_density_ridges(alpha=0.8)+ 
  ggplot2::scale_fill_manual(values = c('uncorrected' = "#ccc3e1", 'aligned' = "#48376f")) + 
  ggridges::theme_ridges() +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "none",
                 text=ggplot2::element_text(size=8),
                 axis.text = ggplot2::element_text(size=10),
                 axis.title = ggplot2::element_text(size=8)) +
  ggplot2::xlab("correlation between cell lines and tumors") +
  ggplot2::ylab('') +
  ggplot2::xlim(0.4, 1)
ggsave(file.path(dir_out, 'correlation_aligned_raw.pdf'), 
       width = 3, height = 2.5)

wilcox.test(dat[dat$group == "uncorrected", "r"], 
            dat[dat$group == "aligned", "r"], 
            paired = TRUE, alternative = "two.sided")

# wilcox p-value < 2.2e-16

