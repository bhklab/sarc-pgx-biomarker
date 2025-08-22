# -----------------------------------------------------------
# Soft-Tissue Sarcoma Data Downloader & Preprocessor
# This script downloads and preprocesses gene expression data 
# for soft-tissue sarcoma from ORCESTRA (https://orcestra.ca/; 
# pharmacogenomics cell line data).
#
#   - ORCESTRA: Retrieves pharmacogenomics profiles for sarcoma cell lines
#     * Extracts gene expression microarray and drug sensitivity data
#     * Filters for soft-tissue sarcoma cell lines
#     * Prepares curated matrices for integrative analysis
#
#   - Outputs: Processed expression and phenotype data ready 
#     for biomarker discovery and integrative modeling
# -----------------------------------------------------------
###################################################################
## Load libraries
###################################################################

library(qs)
library(dplyr)
library(PharmacoGx)
library(MultiAssayExperiment)

##############################################################################
# function to curate/extract expression and drug response data from PSets
##############################################################################
# extract rna microarray expression data
data_extractor_rna <- function(pset, pset.name){
    
    marray <- summarizeMolecularProfiles(pset, mDataType="rna") %>% assay() %>% t() 
    marray <- marray[rowSums(is.na(marray)) != ncol(marray), ] # Removing rows in which the gene expression is NA for all genes
    gene_ann_arr <- featureInfo(pset, mDataType = "rna") %>% data.frame() 
    
    if (!identical(grep("AFFX", colnames(marray)), integer(0))) { 
      marray <- marray[ , -grep("AFFX", colnames(marray))]
      gene_ann_arr <- gene_ann_arr[-grep("AFFX", rownames(gene_ann_arr)) , ]
      
    } 
    
    if(pset.name == "NCI_Sarcoma"){
      
      gene_ann_arr <- gene_ann_arr[gene_ann_arr$gene_biotype == "protein_coding", ]
      gene_ann_arr <- gene_ann_arr[!duplicated(gene_ann_arr$Symbol), ]
      
      int <- intersect(gene_ann_arr$gene_id, colnames(marray))
      marray <- marray[, colnames(marray) %in% int]
      gene_ann_arr <- gene_ann_arr[gene_ann_arr$gene_id %in% int, ]
      colnames(marray) <- gene_ann_arr$Symbol
      
    }else{
      
      gene_ann_arr <- gene_ann_arr[gene_ann_arr$GeneBioType == "protein_coding", ]
      gene_ann_arr <- gene_ann_arr[!duplicated(gene_ann_arr$Symbol), ]
      
      int <- intersect(gene_ann_arr$EnsemblGeneId, colnames(marray))
      marray <- marray[, colnames(marray) %in% int]
      gene_ann_arr <- gene_ann_arr[gene_ann_arr$EnsemblGeneId %in% int, ]
      colnames(marray) <- gene_ann_arr$Symbol
      
    }
    
    
    cell_ann_marray <- cellInfo(pset)
    
    if(pset.name == "NCI_Sarcoma"){
      cell_ann_marray$tissueid <- cell_ann_marray$dimitrios.soft_vs_bone
      cell_ann_marray <- cell_ann_marray[cell_ann_marray$tissueid %in% c("Bone", "Soft"), ]
      cell_ann_marray$tissueid <- ifelse(cell_ann_marray$tissueid == "Soft", "Soft Tissue", cell_ann_marray$tissueid)
    }
    
    cell_ann_marray$mod_tissueid <- ifelse(cell_ann_marray$tissueid == "Uterus", "Soft Tissue",  # Uterus be considered as Soft Tissue
                                        cell_ann_marray$tissueid)
    
    cell_ann_marray <- cell_ann_marray[cell_ann_marray$mod_tissueid != "Other", ]  # remove Other group
    
    marray <- marray[intersect(rownames(cell_ann_marray) , rownames(marray)), ] 
    cell_ann_marray <- cell_ann_marray[intersect(rownames(cell_ann_marray), rownames(marray)), ]
    
    pgx.dat <- list(marray_exp = marray, marray_gene_ann = gene_ann_arr, 
                    cell_ann_marray = cell_ann_marray)
    
    return(pgx.dat)
}

# extract drug response data
data_extractor_drug <- function(pset){
  
  drug_aac <- summarizeSensitivityProfiles(pset, 
                                           sensitivity.measure = "aac_recomputed")
  
  cell_ann_seq <- cellInfo(pset) 
  drug_aac <- drug_aac[, intersect(rownames(cell_ann_seq ) , colnames(drug_aac))] 
  cell_ann_seq <- cell_ann_seq[intersect(rownames(cell_ann_seq), colnames(drug_aac)), ]
  
  pgx.dat <- list( cell_ann_seq = cell_ann_seq, drug_aac = drug_aac)
  
  return(pgx.dat)
  
}

