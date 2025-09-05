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
################################################################################
## Libraries and functions
################################################################################

source("scripts/Prog_celligner_function.R")

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/procdata'  # '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/data'
dir_out <- 'data/results/aligned'  # '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/aligned/'
dir_annot <-  'data/rawdata' # '/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Result/aligned/'
  
########################################################
## Part 1: find optimal number of PCs
########################################################

dat <- qread(file.path(dir_in, "PGx_corrected_gse_rna_sts.qs"))
gene_stats <- calc_gene_stats(dat, dir_annot, hgnc_file = 'hgnc_complete_set.txt') # UPDATED hgnc FILE -> EDIT BY NIKTA

comb_ann <- rbind(
  dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary.Metastasis`) %>%
    dplyr::mutate(type = 'tumor'),
  dat$pset_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary.Metastasis`) %>%
    dplyr::mutate(type = 'CL')
)

TCGA_obj <- create_Seurat_object(exp_mat=dat$TCGA_mat, ann=dat$TCGA_ann, type='tumor')
pset_obj <- create_Seurat_object(dat$pset_mat, dat$pset_ann, type='CL')

TCGA_obj <- cluster_data(seu_obj = TCGA_obj)
pset_obj <- cluster_data(pset_obj)

tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj) #  Zero sample variances detected, have been offset away from zero 
CL_DE_genes <- find_differentially_expressed_genes(pset_obj) #  Zero sample variances detected, have been offset away from zero 

DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
  mutate(
    tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
    CL_rank = dplyr::dense_rank(-gene_stat_CL),
    best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
  dplyr::left_join(gene_stats, by = 'Gene')

# take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
DE_gene_set <- DE_genes %>%
  dplyr::filter(best_rank < global$top_DE_genes_per) %>%
  .[['Gene']]


cov_diff_eig <- run_cPCA(TCGA_obj, pset_obj, global$fast_cPCA)

res <- lapply(1:10, function(k){
  
  print(k)
  if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, 1:k, drop = FALSE]
  } else {
    cur_vecs <- cov_diff_eig$rotation[, 1:k, drop = FALSE]
  }
  
  rownames(cur_vecs) <- colnames(dat$TCGA_mat)
  TCGA_cor <- resid(lm(t(dat$TCGA_mat[complete.cases(dat$TCGA_mat),]) ~ 0 + cur_vecs)) %>% t() #nikta edited
  pset_cor <- resid(lm(t(dat$pset_mat) ~ 0 + cur_vecs)) %>% t()
  
  mnn_res <- run_MNN(pset_cor, TCGA_cor,  k1 = global$mnn_k_tumor, 
                     k2 = global$mnn_k_CL, ndist = global$mnn_ndist,
                     subset_genes = DE_gene_set)
  
  combined_mat <- rbind(mnn_res$corrected, pset_cor)
  
  comb_obj <- create_Seurat_object(combined_mat, comb_ann)
  comb_obj <- cluster_data(seu_obj = comb_obj)
  
  
  list(comb_obj = comb_obj, pairs = mnn_res$pairs)
  
})

res_all <- list(celligner_res = res,
            DE_gene = DE_gene_set[c(1:1000)],
            cpca = cov_diff_eig$sdev[1:10])

qsave(res_all, file=  file.path(dir_out, "pgx_rna_celligner_sts_pcs.qs"))

################################################################################
## Visualization Part 1
################################################################################

dat <- qread(file.path(dir_out, "pgx_rna_celligner_sts_pcs.qs"))

res <- sapply(1:length(dat$celligner_res), function(k){
  
  nrow(dat$celligner_res[[k]]$pairs)
  
})

num_cPCs <- c(1:10)
num_MNN <- c(as.vector(res))

num_cPCs_used <- cbind.data.frame(num_cPCs, num_MNN)
num_cPCs_used$param_type <- ifelse(num_cPCs_used$num_cPCs==7, 'original parameter', 'modified parameter')
num_cPCs_used$num_cPCs <- factor(num_cPCs_used$num_cPCs)
num_cPCs_used$num_MNN <- factor(num_cPCs_used$num_MNN)

cPCs_MNN <- ggplot2::ggplot(num_cPCs_used, ggplot2::aes(num_cPCs, num_MNN, shape=param_type, 
                                                        group= interaction('param_type'))) +
  ggplot2::geom_line(color='gray80') + 
  ggplot2::geom_point() +
  ggplot2::xlab('number of cPCs used') +
  ggplot2::ylab('number of MNN pairs') +
  ggplot2::scale_shape_manual(values=c(`modified parameter`=16, `original parameter`=8)) +
  ggplot2::guides(shape=ggplot2::guide_legend(title="")) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position='bottom', 
                 text = ggplot2::element_text(size=8),
                 axis.title = ggplot2::element_text(size=8), 
                 axis.text = ggplot2::element_text(size=8),
                 legend.margin =ggplot2::margin(0,0,0,0),
                 legend.box.margin=ggplot2::margin(-10,-30,-10,-30)) 
ggsave(file.path(dir_out, "optimal_rna_sts_k.pdf"), 
       width = 4, height = 4)

################################################################################
## Visualization Part 2 
################################################################################

var_res <- cov_diff_eig$sdev^2 / sum(cov_diff_eig$sdev^2)

pdf(file= file.path(dir_out, "cpca_sts_variance.pdf"), width = 5, height = 5)

barplot(var_res[1:10]*100, xlab = "principal component (cPCA)", ylab = "variance explained by cPCs")

dev.off()





