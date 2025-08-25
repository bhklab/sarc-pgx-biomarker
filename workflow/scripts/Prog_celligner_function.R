# -----------------------------------------------------------
# Celligner Alignment & Visualization Script Functions
# This script implements the Celligner pipeline to align and 
# compare gene expression data between tumor samples and 
# pharmacogenomic cell lines.
#
#   Data sources:
#     - GEO (soft-tissue sarcoma tumor datasets)
#     - ORCESTRA (pharmacogenomics cell line datasets)
#
#   Includes:
#     * Seurat object creation, PCA, and UMAP visualization
#     * Batch effect correction using contrastive PCA (cPCA)
#     * Integration of tumor and cell line data using Mutual Nearest Neighbors (MNN)
#     * Identification of differentially expressed genes (DEGs)
#     * Visualization of aligned vs. unaligned datasets
#     * Post-alignment correlation and tissue-specific similarity analysis
#
#   Key methods:
#     - Seurat for normalization, dimensionality reduction, clustering
#     - limma for differential expression analysis
#     - contrastive PCA (cPCA) to remove dataset-specific variation
#     - MNN for cross-dataset alignment (modified from batchelor)
#     - Visualization using ggplot2 and ggridges
#
#   Assumes input data is preprocessed and stored in .qs format
# -----------------------------------------------------------
###################################################################
## Load libraries
###################################################################

library(qs)
require(limma)
require(tibble)
require(plyr)
require(dplyr)
library(Seurat)
library(batchelor)
library(data.table)
library(magrittr)
library(ggplot2)

# ========================================== github codes ==========================================
# ====== Global_params ======
global <- list(
  n_genes = 'all', # set to 'all' to use all protein coding genes found in both datasets 
  umap_n_neighbors = 5, # num nearest neighbors used to create UMAP plot
  umap_min_dist = 0.5, # min distance used to create UMAP plot
  mnn_k_CL = 5, # number of nearest neighbors of tumors in the cell line data #FARNOOSH changed to 5
  mnn_k_tumor = 30, # number of nearest neighbors of cell lines in the tumor data 
  top_DE_genes_per = 1000, # differentially expressed genes with a rank better than this is in the cell line or tumor data
  # are used to identify mutual nearest neighbors in the MNN alignment step
  # remove_cPCA_dims = c(1,2,3,4), # which cPCA dimensions to regress out of the data 
  distance_metric = 'euclidean', # distance metric used for the UMAP projection
  mod_clust_res = 5, # resolution parameter used for clustering the data
  mnn_ndist = 3, # ndist parameter used for MNN
  n_PC_dims = 20, # number of PCs to use for dimensionality reduction #FARNOOSH 
  reduction.use = 'umap', # 2D projection used for plotting
  fast_cPCA = 10 # to run fast cPCA (approximate the cPCA eigenvectors instead of calculating all) set this to a value >= 4
)

################################################################################
## calculate gene average expression and variance for an expression matrix
################################################################################
calc_gene_stats <- function(dat, data_dir, hgnc_file) { 
  common_genes <- intersect(colnames(dat$TCGA_mat), colnames(dat$pset_mat))
  
  hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
  hgnc.complete.set <- hgnc.complete.set %>% 
    dplyr::select(Gene = ensembl_gene_id, Symbol = symbol) %>%
    dplyr::filter(Symbol %in% common_genes)  
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$Symbol),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$Symbol
  hgnc.complete.set <- hgnc.complete.set[common_genes,]  
  
  gene_stats <- data.frame(
    Tumor_SD = apply(dat$TCGA_mat, 2, sd, na.rm=T),
    pset_SD = apply(dat$pset_mat, 2, sd, na.rm=T),
    Tumor_mean = colMeans(dat$TCGA_mat, na.rm=T),
    pset_mean = colMeans(dat$pset_mat, na.rm=T),
    Symbol = common_genes,
    stringsAsFactors = F) %>% 
    dplyr::mutate(max_SD = pmax(Tumor_SD, pset_SD, na.rm=T)) #add avg and max SD per gene
  
  gene_stats <- dplyr::left_join(hgnc.complete.set, gene_stats, by = "Symbol") 
  
  return(gene_stats)
  
}

################################################################################
## create seurat objects given an expression matrix and annotation table
################################################################################
create_Seurat_object <- function(exp_mat, ann, type = NULL) { 
  seu_obj <- Seurat::CreateSeuratObject(t(exp_mat),
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = ann %>%
                                          magrittr::set_rownames(ann$sampleID))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # normalize, mean center the data, important for PCA
  seu_obj <- Seurat::NormalizeData(seu_obj)

  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)

  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                              features = rownames(Seurat::GetAssayData(seu_obj)),
                              npcs = global$n_PC_dims, verbose = F)
  
  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                               reduction = 'pca',
                               n.neighbors = global$umap_n_neighbors,
                               min.dist =  global$umap_min_dist,
                               metric = global$distance_metric, verbose=F)
  
  return(seu_obj)
}

################################################################################
## take in a Seurat object and run default Seurat clustering algorithm
################################################################################

cluster_data <- function(seu_obj) {
  seu_obj <- Seurat::FindNeighbors(seu_obj, 
                                   reduction = 'pca',
                                   dims = 1:global$n_PC_dims,
                                   k.param = 5, # it was 20  FARNOOSH
                                   force.recalc = TRUE,
                                   verbose = FALSE)
  
  seu_obj %<>% Seurat::FindClusters(reduction = 'pca', 
                                    resolution = global$mod_clust_res)
  
  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters
  
  return(seu_obj)
  
}

################################################################################
## find differentially expressed genes between clusters within the data

# Estimate linear-model stats for a matrix of data with respect to a group of phenotype variables
# using limma with empirical Bayes moderated F-stats for p-values

################################################################################

run_lm_stats_limma_group <- function (mat, phenos, covars = NULL, weights = NULL, target_type = "Gene", 
                                      limma_trend = FALSE) {
  udata <- rownames(mat) %>% intersect(rownames(phenos))
  if (!is.null(covars)) {
    udata %<>% intersect(rownames(covars))
  }
  form <- as.formula(paste("~", paste0(colnames(phenos), collapse = " + ")))
  design <- model.matrix(form, data = phenos[udata, , drop = F])
  if (!is.null(covars)) {
    covars <- data.frame(covars)
    form <- as.formula(paste("~", paste0(colnames(covars), 
                                         collapse = " + ")))
    Cdesign <- model.matrix(form, data = covars[udata, , 
                                                drop = F])
    Cdesign <- Cdesign[, setdiff(colnames(Cdesign), "(Intercept)"), 
                       drop = FALSE]
    stopifnot(length(intersect(colnames(Cdesign), colnames(design))) == 
                0)
    design %<>% cbind(Cdesign)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata, ])
    }
    else {
      weights <- weights[udata]
    }
  }
  design <- design[, colSums(design) > 2, drop = FALSE]
  targ_coefs <- setdiff(colnames(design), "(Intercept)")
  fit <- limma::lmFit(t(mat[udata, ]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- which(colnames(design) %in% targ_coefs)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf, 
                             sort.by = "F", genelist = colnames(mat))
  results %<>% tibble::rownames_to_column(var = target_type)
  results %<>% magrittr::set_colnames(revalue(colnames(.), c(AveExpr = "Avg", 
                                                             F = "F_stat", 
                                                             P.Value = "p.value", 
                                                             adj.P.Val = "q.value"))) %>% 
    na.omit() %>% dplyr::select(-ProbeID)
  return(results)
}

# ====== find_differentially_expressed_genes ====== 

find_differentially_expressed_genes <- function(seu_obj) {
  n_clusts <- nlevels(seu_obj@meta.data$seurat_clusters)
  if (n_clusts > 2) {
    cur_DE_genes <- run_lm_stats_limma_group(
      t(Seurat::GetAssayData(seu_obj, assay='RNA', layer = 'scale.data')), # chnage slot to layer
      seu_obj@meta.data %>% dplyr::select(seurat_clusters),
      limma_trend = TRUE) %>%
      dplyr::select(Gene, gene_stat = F_stat)
  } else if (n_clusts == 2) {
    cur_DE_genes <- run_lm_stats_limma(t(Seurat::GetAssayData(seu_obj, assay='RNA', layer = 'scale.data')),
                                       seu_obj@meta.data$cluster,
                                       limma_trend = TRUE) %>%
      dplyr::mutate(gene_stat = abs(t_stat)) %>%
      dplyr::select(Gene, gene_stat)
  } else {
    cur_DE_genes <- data.frame(Gene = colnames(seu_obj), gene_stat = NA)
  }
  
  return(cur_DE_genes)
  
}

################################################################################
## run contrastive principal components analysis
################################################################################
# ====== calculate the average expression per cluster  ======

get_cluster_averages <- function(mat, cluster_df) {
  n_clusts <- nlevels(cluster_df$seurat_clusters)
  clust_avgs <- matrix(NA, nrow = n_clusts, ncol = ncol(mat)) %>% 
    magrittr::set_colnames(colnames(mat)) %>% 
    magrittr::set_rownames(levels(cluster_df$seurat_clusters))
  for (ii in levels(cluster_df$seurat_clusters)) {
    clust_avgs[ii,] <- colMeans(mat[cluster_df$seurat_clusters == ii,], na.rm=T)
  }
  return(clust_avgs)
}

# ======  contrastive principal components analysis ====== 

run_cPCA_analysis <- function(TCGA_dat, pset_dat, tumor_cluster_df, CL_cluster_df, pc_dims=NULL) {
  tumor_clust_avgs <- get_cluster_averages(TCGA_dat, tumor_cluster_df)
  CL_clust_avgs <- get_cluster_averages(pset_dat, CL_cluster_df)
  
  TCGA_subtype_ms <- TCGA_dat - tumor_clust_avgs[tumor_cluster_df$seurat_clusters,]
  pset_subtype_ms <- pset_dat - CL_clust_avgs[CL_cluster_df$seurat_clusters,]
  
  TCGA_cov <- cov(TCGA_subtype_ms)
  pset_cov <- cov(pset_subtype_ms)
  
  if(!is.null(pc_dims)) {
    cov_diff_eig <- irlba::prcomp_irlba(TCGA_cov - pset_cov, n = pc_dims)
  } else {
    cov_diff_eig <- eigen(TCGA_cov - pset_cov)
  }
  return(cov_diff_eig)
}

# ====== run_cPCA ======
run_cPCA <- function(TCGA_obj, pset_obj, pc_dims = NULL) {
  cov_diff_eig <- run_cPCA_analysis(t(Seurat::GetAssayData(TCGA_obj, assay='RNA', layer = 'scale.data')), 
                                    t(Seurat::GetAssayData(pset_obj, assay='RNA', layer = 'scale.data')), 
                                    TCGA_obj@meta.data, pset_obj@meta.data, pc_dims=pc_dims)
  return(cov_diff_eig) 
}


################################################################################
# mutual nearest neighbors object with corrected data
################################################################################
# ====== run mutual nearest neighbors batch correction  ======
# Allows for separate k values per dataset, and simplifies some of the IO and doesn't use PCA reduction
modified_mnnCorrect <- function(ref_mat, 
                                targ_mat, 
                                k1 = 20, # celligner_global$mnn_k_tumor
                                k2 = 20, # celligner_global$mnn_k_CL
                                ndist = 3, 
                                subset_genes = NULL) {
  if (is.null(subset_genes)) {
    subset_genes <- colnames(ref_mat) 
  }  
  
  sets <- batchelor::findMutualNN(ref_mat[, subset_genes], 
                                  targ_mat[, subset_genes], 
                                  k1 = k2, k2 = k1, 
                                  BPPARAM = BiocParallel::SerialParam())
  mnn_pairs <- as.data.frame(sets) %>% 
    dplyr::mutate(ref_ID = rownames(ref_mat)[first],
                  targ_ID = rownames(targ_mat)[second],
                  pair = seq(nrow(.))) %>% 
    dplyr::select(-first, -second)
  
  # Estimate the overall batch vector.
  ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  overall.batch <- colMeans(ave.out$averaged)
  
  #remove variation along the overall batch vector
  ref_mat <- .center_along_batch_vector(ref_mat, overall.batch)
  targ_mat <- .center_along_batch_vector(targ_mat, overall.batch)
  
  # Recompute correction vectors and apply them.
  re.ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  targ_mat <- .tricube_weighted_correction(targ_mat, re.ave.out$averaged, re.ave.out$second, k=k2, ndist=ndist, subset_genes, BPPARAM=BiocParallel::SerialParam())
  
  final <- list(corrected = targ_mat, 
                pairs = mnn_pairs)
  return(final)
}

# ====== run mutual nearest neighbors ======
run_MNN <- function(pset_cor, TCGA_cor,  
                    k1 = global$mnn_k_tumor, 
                    k2 = global$mnn_k_CL, 
                    ndist = global$mnn_ndist, 
                    subset_genes) {
  mnn_res <- modified_mnnCorrect(pset_cor, TCGA_cor, k1 = k1, k2 = k2, ndist = ndist, 
                                 subset_genes = subset_genes)
  
  return(mnn_res)
}

################################################################################
## correlation between cell lines and tumor
################################################################################
calc_tumor_CL_cor <- function(Celligner_aligned_data, Celligner_info) {
  tumors_samples <- dplyr::filter(Celligner_info, type=='tumor')$sampleID
  cl_samples <- dplyr::filter(Celligner_info, type=='CL')$sampleID
  tumor_CL_cor <- cor(t(Celligner_aligned_data[tumor_samples,]), t(Celligner_aligned_data[cl_samples,]),
                      use='pairwise')
  
  
  return(tumor_CL_cor)
}

################################################################################
## additional functions for Celligner method
################################################################################


# ====== .average_correction  ======
# Copied from dev version of scran (2018-10-28) with slight modifications as noted
#https://github.com/MarioniLab/scran
.average_correction <- function(refdata, mnn1, curdata, mnn2)
  # Computes correction vectors for each MNN pair, and then
  # averages them for each MNN-involved cell in the second batch.
{
  corvec <- refdata[mnn1,,drop=FALSE] - curdata[mnn2,,drop=FALSE]
  corvec <- rowsum(corvec, mnn2)
  npairs <- table(mnn2)
  stopifnot(identical(names(npairs), rownames(corvec)))
  corvec <- unname(corvec)/as.vector(npairs)
  list(averaged=corvec, second=as.integer(names(npairs)))
}

# ====== .center_along_batch_vector  ======
.center_along_batch_vector <- function(mat, batch.vec) 
  # Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
  # This removes any variation along the overall batch vector within each matrix.
{
  batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
  batch.loc <- as.vector(mat %*% batch.vec)
  central.loc <- mean(batch.loc)
  mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
  return(mat)
}

# ====== .tricube_weighted_correction  ======
#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocParallel SerialParam
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, subset_genes, BNPARAM=NULL, BPPARAM=BiocParallel::SerialParam())
  # Computing tricube-weighted correction vectors for individual cells,
  # using the nearest neighbouring cells _involved in MNN pairs_.
  # Modified to use FNN rather than queryKNN for nearest neighbor finding
{
  cur.uniq <- curdata[in.mnn,,drop=FALSE]
  safe.k <- min(k, nrow(cur.uniq))
  # closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
  closest <- FNN::get.knnx(cur.uniq[, subset_genes], query=curdata[, subset_genes], k=safe.k)
  # weighted.correction <- .compute_tricube_average(correction, closest$index, closest$distance, ndist=ndist)
  weighted.correction <- .compute_tricube_average(correction, closest$nn.index, closest$nn.dist, ndist=ndist)
  curdata + weighted.correction
}

# ====== .compute_tricube_average ======
.compute_tricube_average <- function(vals, indices, distances, bandwidth=NULL, ndist=3) 
  # Centralized function to compute tricube averages.
  # Bandwidth is set at 'ndist' times the median distance, if not specified.
{
  if (is.null(bandwidth)) {
    middle <- ceiling(ncol(indices)/2L)
    mid.dist <- distances[,middle]
    bandwidth <- mid.dist * ndist
  }
  bandwidth <- pmax(1e-8, bandwidth)
  
  rel.dist <- distances/bandwidth
  rel.dist[rel.dist > 1] <- 1 # don't use pmin(), as this destroys dimensions.
  tricube <- (1 - rel.dist^3)^3
  weight <- tricube/rowSums(tricube)
  
  output <- 0
  for (kdx in seq_len(ncol(indices))) {
    output <- output + vals[indices[,kdx],,drop=FALSE] * weight[,kdx]
  }
  
  if (is.null(dim(output))) {
    matrix(0, nrow(vals), ncol(vals))
  } else {
    output
  }
}


################################################################################
# run Celligner method 
################################################################################

run_Celligner <- function(data_dir, dat, remove_cPCA_dims = c(1,2,3,4), plotname) {
  
  dat <- qread(dat)
  gene_stats <- calc_gene_stats(dat, data_dir, hgnc_file = 'hgnc_complete_set.txt') 
  
  comb_ann <- rbind(
    dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary.Metastasis`) %>%
      dplyr::mutate(type = 'tumor'),
    dat$pset_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary.Metastasis`) %>%
      dplyr::mutate(type = 'CL')
  )
  
  TCGA_obj <- create_Seurat_object(exp_mat=dat$TCGA_mat, ann=dat$TCGA_ann, type='tumor')
  pset_obj <- create_Seurat_object(exp_mat=dat$pset_mat, ann=dat$pset_ann, type='CL')
  
  TCGA_obj <- cluster_data(seu_obj = TCGA_obj)
  pset_obj <- cluster_data(pset_obj)
  
  tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj) # Warning message: Zero sample variances detected, have been offset away from zero 
  CL_DE_genes <- find_differentially_expressed_genes(pset_obj) # Warning message: Zero sample variances detected, have been offset away from zero 
  
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
  
  if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, remove_cPCA_dims, drop = FALSE]
  } else {
    cur_vecs <- cov_diff_eig$rotation[, remove_cPCA_dims, drop = FALSE]
  }
  
  rownames(cur_vecs) <- colnames(dat$TCGA_mat)
  TCGA_cor <- resid(lm(t(dat$TCGA_mat[complete.cases(dat$TCGA_mat),]) ~ 0 + cur_vecs)) %>% t() #nikta edited
  pset_cor <- resid(lm(t(dat$pset_mat) ~ 0 + cur_vecs)) %>% t()
  
  mnn_res <- run_MNN(pset_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist,
                     subset_genes = DE_gene_set)
  
  combined_mat <- rbind(mnn_res$corrected, pset_cor)
  
  comb_obj <- create_Seurat_object(combined_mat, comb_ann)
  comb_obj <- cluster_data(seu_obj = comb_obj)
  
  res <- list(comb_obj = comb_obj,
              DE_gene = DE_gene_set[c(1:1000)],
              cpca = cov_diff_eig$sdev)
  
  return(res) 
  
}

################################################################################
################################################################################
# ============================ visualization function ==========================
################################################################################
################################################################################
# figure: uncorrected data: CL vs tumor patient 

plot_uncorrected_data_type <- function(org_dat, before_plot, dir_output) {
  
  comb_ann <- cbind.data.frame(`sampleID` = c(rownames(org_dat$TCGA_mat), rownames(org_dat$pset_mat)),
                               `type` = c(rep('tumor', nrow(org_dat$TCGA_mat)), 
                                          rep("CL", nrow(org_dat$pset_mat))),
                               `Type` = c(org_dat$TCGA_ann$Type, org_dat$pset_ann$Type))
  
  comb_ann$type <- ifelse(comb_ann$type == "tumor", "GEO", comb_ann$type) 
  
  # using Seurat object to run cPCA and UMAP
  original_combined_obj <-  Seurat::CreateSeuratObject(t(rbind(org_dat$TCGA_mat,
                                                               org_dat$pset_mat)),
                                                       min.cells = 0,
                                                       min.features = 0,
                                                       meta.data = comb_ann %>%
                                                         magrittr::set_rownames(comb_ann$sampleID))
  
  original_combined_obj <- Seurat::NormalizeData(original_combined_obj)
  original_combined_obj <- Seurat::ScaleData(original_combined_obj, 
                                             features = rownames(Seurat::GetAssayData(original_combined_obj)), 
                                             do.scale = F)
  
  original_combined_obj %<>% Seurat::RunPCA(assay='RNA',
                                            features = rownames(Seurat::GetAssayData(original_combined_obj)),
                                            npcs = global$n_PC_dims,
                                            verbose = F)
  
  original_combined_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                             reduction = 'pca',
                                             n.neighbors = global$umap_n_neighbors,
                                             min.dist = global$umap_min_dist,
                                             metric = global$distance_metric,
                                             verbose=F)
  
  uncorrected_alignment <- Seurat::Embeddings(original_combined_obj, reduction = 'umap') %>%
    as.data.frame() %>%
    set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    rownames_to_column(var = 'sampleID') %>%
    left_join(comb_ann, by = 'sampleID')
  
  ggplot2::ggplot(uncorrected_alignment, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, type=='GEO'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=type)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, type=='CL'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=type), stroke=0.5) +
    ggplot2::scale_color_manual(values=c(CL="#2166ac")) +
    ggplot2::scale_fill_manual(values=c('GEO'= "#993404")) + 
    ggplot2::xlab('UMAP 1') + ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='right',
                   legend.title = ggplot2::element_blank(),   # removes "type"
                   text = ggplot2::element_text(size=8),
                   axis.text = ggplot2::element_text(size=6),
                   axis.title = ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0), 
                   #legend.box.margin=ggplot2::margin(-10,-30,-10,-30),
                   legend.box.margin=ggplot2::margin(0,0,0,0),
                   axis.line = ggplot2::element_line(linewidth = .3))
  
  ggsave(file.path(dir_output, paste(before_plot, '.pdf', sep="")), 
         width = 4, height = 3)
  
}   

################################################################################
# figure: uncorrected data: CL vs tumor patient ---> study
################################################################################

plot_uncorrected_data_study <- function(org_dat, before_plot, dir_output) {
  
  comb_ann <- cbind.data.frame(`sampleID` = c(rownames(org_dat$TCGA_mat), 
                                              rownames(org_dat$pset_mat)),
                               `type` = c(rep('tumor', nrow(org_dat$TCGA_mat)), 
                                          rep("CL", nrow(org_dat$pset_mat))),
                               `Type` = c(org_dat$TCGA_ann$Type, org_dat$pset_ann$Type))
  
  comb_ann$Type <- ifelse(comb_ann$Type ==  "cohort 1-GSE21050" | comb_ann$Type ==  "cohort 2-GSE21050",
                          "GSE21050", comb_ann$Type)
  
  # using Seurat object to run cPCA and UMAP
  original_combined_obj <-  Seurat::CreateSeuratObject(t(rbind(org_dat$TCGA_mat,
                                                               org_dat$pset_mat)),
                                                       min.cells = 0,
                                                       min.features = 0,
                                                       meta.data = comb_ann %>%
                                                         magrittr::set_rownames(comb_ann$sampleID))
  
  original_combined_obj <- Seurat::NormalizeData(original_combined_obj)
  original_combined_obj <- Seurat::ScaleData(original_combined_obj, 
                                             features = rownames(Seurat::GetAssayData(original_combined_obj)), 
                                             do.scale = F)
  
  original_combined_obj %<>% Seurat::RunPCA(assay='RNA',
                                            features = rownames(Seurat::GetAssayData(original_combined_obj)),
                                            npcs = global$n_PC_dims,
                                            verbose = F)
  
  original_combined_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                             reduction = 'pca',
                                             n.neighbors = global$umap_n_neighbors,
                                             min.dist = global$umap_min_dist,
                                             metric = global$distance_metric,
                                             verbose=F)
  
  uncorrected_alignment <- Seurat::Embeddings(original_combined_obj, reduction = 'umap') %>%
    as.data.frame() %>%
    set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    rownames_to_column(var = 'sampleID') %>%
    left_join(comb_ann, by = 'sampleID')
  
  ggplot2::ggplot(uncorrected_alignment, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, Type =='GSE21050'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=Type)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, Type =='GSE21122'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=Type)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, Type =='GSE30929'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=Type)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, Type =='CCLE'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=Type), stroke=0.5) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, Type =='GDSC'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=Type), stroke=0.5) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, Type =='NCI-Sarcoma'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=Type), stroke=0.5) +
    ggplot2::scale_color_manual(values=c('CCLE'="#1b7837", 
                                         'GDSC' = "#542788", 
                                         'NCI-Sarcoma' = "#bf812d") )+
    ggplot2::scale_fill_manual(values=c('GSE21050' = "#4d4d4d",
                                        'GSE21122' ="#a50026",
                                        'GSE30929' = "#c994c7"
    )) +
    ggplot2::xlab('UMAP 1') + ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='right',
                   legend.title = ggplot2::element_blank(),   # removes "type"
                   text = ggplot2::element_text(size=8),
                   axis.text = ggplot2::element_text(size=6),
                   axis.title = ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0), 
                   legend.box.margin=ggplot2::margin(0,0,0,0),
                   axis.line = ggplot2::element_line(linewidth = .3))
  
  
  ggsave(file.path(dir_output, paste(before_plot, '.pdf', sep="")), 
         width = 4, height = 3)
}   

################################################################################
# figure: corrected data: CL vs tumor patient 
################################################################################

Celligner_alignment_plot_type <- function(aligned_data, after_plot, dir_output) {
  
  aligned_data %<>% Seurat::RunUMAP(assay = 'RNA', 
                                    dims = 1:global$n_PC_dims,
                                    reduction = 'pca',
                                    n.neighbors = global$umap_n_neighbors,
                                    min.dist = global$umap_min_dist,
                                    metric = global$distance_metric,
                                    verbose=F)
  
  comb_ann <- aligned_data[[]] # This contains all the metadata
  
  
  alignment <- Seurat::Embeddings(aligned_data, reduction = 'umap') %>%
    as.data.frame() %>%
    set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    rownames_to_column(var = 'sampleID') %>%
    left_join(comb_ann, by = 'sampleID') 
  
  alignment$type <- ifelse(alignment$type == "tumor", "GEO", "CL") 
  
  ggplot2::ggplot(alignment, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = filter(alignment, type=='GEO'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=type)) +
    ggplot2::geom_point(data = filter(alignment, type=='CL'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=type), stroke=0.5) +
    ggplot2::scale_color_manual(values=c(CL="#2166ac")) +
    ggplot2::scale_fill_manual(values=c('GEO'= "#993404")) +
    ggplot2::xlab('UMAP 1') + ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='right',
                   legend.title = ggplot2::element_blank(),   # removes "type"
                   text = ggplot2::element_text(size=8),
                   axis.text = ggplot2::element_text(size=6),
                   axis.title = ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0), 
                   legend.box.margin=ggplot2::margin(0,0,0,0),
                   axis.line = ggplot2::element_line(linewidth = .3))
  
  
  ggsave(file.path(dir_output, paste(after_plot, '.pdf', sep="")), 
         width = 4, height = 3)
  
}


################################################################################
# figure: corrected data: CL vs tumor patient ----> study
################################################################################

Celligner_alignment_plot_study <- function(aligned_data, org_dat, after_plot, dir_output) {
  
  aligned_data %<>% Seurat::RunUMAP(assay = 'RNA', 
                                    dims = 1:global$n_PC_dims,
                                    reduction = 'pca',
                                    n.neighbors = global$umap_n_neighbors,
                                    min.dist = global$umap_min_dist,
                                    metric = global$distance_metric,
                                    verbose=F)
  
  comb_ann <- aligned_data[[]] # This contains all the metadata
  
  alignment <- Seurat::Embeddings(aligned_data, reduction = 'umap') %>%
    as.data.frame() %>%
    set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    rownames_to_column(var = 'sampleID') %>%
    left_join(comb_ann, by = 'sampleID') 
  
  alignment$Type <-  c(org_dat$TCGA_ann$Type, org_dat$pset_ann$Type)
  alignment$Type <- ifelse(alignment$Type ==  "cohort 1-GSE21050" | alignment$Type ==  "cohort 2-GSE21050",
                           "GSE21050", alignment$Type)
  
  ggplot2::ggplot(alignment, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = filter(alignment, Type =='GSE21050'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=Type)) +
    ggplot2::geom_point(data = filter(alignment, Type =='GSE21122'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=Type)) +
    ggplot2::geom_point(data = filter(alignment, Type =='GSE30929'), 
                        alpha=0.6, size=1.7, pch=21, color='white', aes(fill=Type)) +
    ggplot2::geom_point(data = filter(alignment, Type =='CCLE'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=Type), stroke=0.5) +
    ggplot2::geom_point(data = filter(alignment, Type =='GDSC'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=Type), stroke=0.5) +
    ggplot2::geom_point(data = filter(alignment, Type =='NCI-Sarcoma'), 
                        alpha=0.6, size=1.7, pch=4, aes(color=Type), stroke=0.5) +
    ggplot2::scale_color_manual(values=c('CCLE'="#1b7837", 
                                         'GDSC' = "#542788", 
                                         'NCI-Sarcoma' = "#bf812d") )+
    ggplot2::scale_fill_manual(values=c('GSE21050' = "#4d4d4d",
                                        'GSE21122' ="#a50026",
                                        'GSE30929' = "#c994c7"
    )) +
    ggplot2::xlab('UMAP 1') + ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='right',
                   legend.title = ggplot2::element_blank(),   # removes "type"
                   text = ggplot2::element_text(size=8),
                   axis.text = ggplot2::element_text(size=6),
                   axis.title = ggplot2::element_text(size=8),
                   legend.margin =ggplot2::margin(0,0,0,0), 
                   legend.box.margin=ggplot2::margin(0,0,0,0),
                   axis.line = ggplot2::element_line(linewidth = .3))

  ggsave(file.path(dir_output, paste(after_plot, '.pdf', sep="")), 
         width = 4, height = 3)
  
}

################################################################################
# figure: post-hoc analysis (correlation)
################################################################################
# pairwise correlation between tumor and CL (uncorrected data)
calc_uncorrected_tumor_CL_correlation <- function(TCGA_mat, pset_mat) {
  
  uncorrected_tumor_CL_cor <- cor(t(TCGA_mat), t(pset_mat), use='pairwise')
  
  return(uncorrected_tumor_CL_cor)
}


# pairwise correlation between tumor and CL (corrected data)
calc_corrected_tumor_CL_correlation <- function(alignment) {
  
  alignment <- cbind.data.frame(`sampleID` = alignment$sampleID,
                                `type` = alignment$type,
                                `lineage` = alignment$lineage,
                                t(as.matrix(Seurat::GetAssayData(alignment))))
  
  TCGA_mat <- alignment[alignment$type == "tumor", -c(1,2,3)]
  pset_mat <- alignment[alignment$type == "CL", -c(1,2,3)]
  
  corrected_tumor_CL_cor <- cor(t(TCGA_mat), t(pset_mat), use='pairwise')
  
  return(corrected_tumor_CL_cor)
}

# create distance function (aligned data)
cell_line_tumor_distance_distribution <- function(alignment, tumor_CL_cor, name_fig, dir_output) {
  
  alignment <- cbind.data.frame(`sampleID` = alignment$sampleID,
                                `type` = alignment$type,
                                `lineage` = alignment$lineage,
                                t(as.matrix(Seurat::GetAssayData(alignment))))
  rownames(alignment) <- NULL
  alignment$compare_types <- alignment$lineage
  
  common_cancer_types <- intersect(dplyr::filter(alignment, type=='tumor')$compare_types, 
                                   dplyr::filter(alignment, type=='CL')$compare_types)
  
  tumor_names <- character()
  CL_names <- character()
  dist_list <- numeric()
  tissue_types <- character()
  for(cancer in common_cancer_types) {
    cur_tumors <- dplyr::filter(alignment, type=='tumor' & compare_types==cancer)$sampleID
    cur_CLs <- dplyr::filter(alignment, type=='CL' & compare_types==cancer)$sampleID
    cur_dist <- reshape2::melt(as.matrix(tumor_CL_cor[cur_tumors, cur_CLs]))
    tumor_names <- c(tumor_names, as.character(cur_dist$Var1))
    CL_names <- c(CL_names, as.character(cur_dist$Var2))
    dist_list <- c(dist_list, cur_dist$value)
    tissue_types <- c(tissue_types, rep(cancer, nrow(cur_dist)))
    
  }
  
  dist_df <- cbind.data.frame(tumor_names, CL_names, dist_list, tissue_types)
  dist_df$tissue_types <- gsub("_", " ", dist_df$tissue_types)
  mean_dist <- aggregate(dist_df$dist_list, list(dist_df$tissue_types), 
                         FUN = quantile, probs = 0.25) %>% dplyr::arrange(desc(x))
  
  mean_dist$Group.1 <- rev(mean_dist$Group.1)
  dist_df$tissue_types <- factor(dist_df$tissue_types) # levels = mean_dist$Group.1
  
  tumor_dist_spread <- ggplot2::ggplot(dplyr::filter(dist_df, tissue_types != 'all'),
                                       ggplot2::aes(x = dist_list, y = tissue_types, fill = tissue_types)) +
    ggridges::geom_density_ridges(alpha=0.8)+ 
    ggplot2::scale_fill_manual(values = c('Soft Tissue' = "#456b53ff")) + 
    ggridges::theme_ridges() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   text=ggplot2::element_text(size=10),
                   axis.text = ggplot2::element_text(size=10),
                   axis.title = ggplot2::element_text(size=8)) +
    ggplot2::xlab("correlation between cell lines and tumors") +
    ggplot2::ylab('') +
    ggplot2::xlim(0.4, 1)
  
  ggsave(file.path(dir_output, paste(name_fig, '.pdf', sep="")), 
         width = 3, height = 2.5)
  
}


# plot distribution (uncorrected data)
plot_uncorrected_distribution_of_CL_tumor_distances <- function(uncorrected_tumor_CL_cor, alignment, name_fig, dir_output) {
  
  alignment <- cbind.data.frame(`sampleID` = alignment$sampleID,
                                `type` = alignment$type,
                                `lineage` = alignment$lineage,
                                t(as.matrix(Seurat::GetAssayData(alignment))))
  rownames(alignment) <- NULL
  
  alignment$compare_types <- alignment$lineage
  
  common_cancer_types <- intersect(dplyr::filter(alignment, type=='tumor')$compare_types, 
                                   dplyr::filter(alignment, type=='CL')$compare_types)
  
  tumor_names <- character()
  CL_names <- character()
  dist_list <- numeric()
  tissue_types <- character()
  for(cancer in common_cancer_types) {
    cur_tumors <- dplyr::filter(alignment, type=='tumor' & compare_types==cancer)$sampleID
    cur_CLs <- dplyr::filter(alignment, type=='CL' & compare_types==cancer)$sampleID
    cur_dist <- reshape2::melt(as.matrix(uncorrected_tumor_CL_cor[cur_tumors, cur_CLs]))
    tumor_names <- c(tumor_names, as.character(cur_dist$Var1))
    CL_names <- c(CL_names, as.character(cur_dist$Var2))
    dist_list <- c(dist_list, cur_dist$value)
    tissue_types <- c(tissue_types, rep(cancer, nrow(cur_dist)))
    
  }
  
  dist_df <- cbind.data.frame(tumor_names, CL_names, dist_list, tissue_types)
  dist_df$tissue_types <- gsub("_", " ", dist_df$tissue_types)
  
  mean_dist <- aggregate(dist_df$dist_list, list(dist_df$tissue_types), 
                         FUN = quantile, probs = 0.25) %>% dplyr::arrange(desc(x))
  
  mean_dist$Group.1 <- rev(mean_dist$Group.1)
  dist_df$tissue_types <- as.factor(dist_df$tissue_types)
  
  ggplot2::ggplot(dplyr::filter(dist_df, tissue_types != 'all'),
                  ggplot2::aes(x = dist_list, y = tissue_types, fill = tissue_types)) + 
    ggridges::geom_density_ridges(alpha=0.8) + 
    ggplot2::scale_fill_manual(values = c('Soft Tissue' = "#a09f9bff")) + 
    ggridges::theme_ridges() + 
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   text=ggplot2::element_text(size=10),
                   axis.text = ggplot2::element_text(size=10),
                   axis.title = ggplot2::element_text(size=8)) +
    ggplot2::xlab("correlation between cell lines and tumors") + 
    ggplot2::ylab('') +
    ggplot2::xlim(0.4, 1)
  
  ggsave(file.path(dir_output, paste(name_fig, '.pdf', sep="")), 
         width = 3, height = 2.5)
  
}
