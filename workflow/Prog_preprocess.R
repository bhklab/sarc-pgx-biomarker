# -----------------------------------------------------------
# Soft-Tissue Sarcoma Data Downloader & Preprocessor
# This script downloads and preprocesses gene expression data 
# for soft-tissue sarcoma from both GEO (clinical samples) 
# and ORCESTRA (https://orcestra.ca/; pharmacogenomics cell line data).
#
#   - GEO: Downloads and curated datasets (GSE21122, GSE21050, GSE30929)
#     * Parses expression matrices and metadata
#     * Performs basic QC and formatting
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

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/rawdata'    # '/home/bioinf/bhklab/farnoosh/SARC/Data'
dir_out <- 'data/procdata'  # '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/data'
  
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

#################################################################
## Load pre-clinical data
#################################################################

dat <- qread(file = file.path(dir_in, "sarcsets.qs"))
names(dat) <- c("CCLE", "CTRP", "gCSI",  "GDSCv2", "GDSCv1", "NCI60", "PRISM", "NCI_Sarcoma")

for(i in 1:length(dat)){
  dat[[i]] <- updateObject(dat[[i]])
}

################################################
## Extract microarray and drug response data
################################################

pset_name <- c("CCLE", "GDSCv1", "GDSCv2", "NCI_Sarcoma")

# smicroarray (rna) data
rna_data <- lapply(1:length(pset_name), function(k){
  
  data_extractor_rna(dat[[which(names(dat) == pset_name[k])]], pset.name = pset_name[k])
  
})

names(rna_data) <- pset_name


# drug response data
pset_name <- c("CCLE", "CTRP", "GDSCv1", "GDSCv2", "NCI_Sarcoma")
drug_data <-  lapply(1:length(pset_name), function(k){
  
  data_extractor_drug(dat[[which(names(dat) == pset_name[k])]])
  
})

names(drug_data) <- pset_name

# update NCI AAC data 
aac_nci <- qread(file.path(dir_in, 'NCI_aac_data.qs')) # "/home/bioinf/bhklab/farnoosh/SARC/PGx_sts/Data/NCI_aac_data.qs"
drug_data$NCI$cell_ann_seq <- aac_nci$annot_cl
drug_data$NCI$drug_aac <- aac_nci$aac

dat <- list(rna = rna_data, aac = drug_data)
qsave(dat, file= file.path(dir_out, "PGx_sarc_data.qs"))

#######################################################
## GEO STS datasets
#######################################################
### GSE21122 (microarray, STS)
dat <- qread(file= file.path(dir_in, "GSE21122_RangedSummarizedExperiment_2025-07-28.qs")) # old version: 2022-06-11
          
exp_mat <- assay(dat)
exp_mat <- t(exp_mat)

cell_ann <- data.frame(colData(dat))
cell_ann <- cell_ann[cell_ann$characteristics_ch1.1 != "disease status: Normal control", ] # remove 9 normal tissues
exp_mat <- exp_mat[rownames(exp_mat) %in% cell_ann$geo_accession, ]

gene_ann <- data.frame(rowData(dat))
gene_ann <- gene_ann[gene_ann$gene_type == "protein_coding", ]
gene_ann <- gene_ann[!duplicated(gene_ann$gene_name), ]
exp_mat <- exp_mat[, colnames(exp_mat) %in% gene_ann$gene_id]
colnames(exp_mat) <- gene_ann$gene_name

dat <- list(marray_exp = exp_mat, marray_gene_ann = gene_ann, 
            cell_ann_marray = cell_ann)

qsave(dat, file=  file.path(dir_out, "GSE21122_sarc_data.qs"))
# 149 samples, 11,245 genes 

### GSE21050 (microarray, STS)

dat <- qread(file= file.path(dir_in, "GSE21050_SummarizedExperiment_2025-07-28.qs")) # old version: 2022-06-11
exp_mat <- assay(dat)
exp_mat <- t(exp_mat)

cell_ann <- data.frame(colData(dat))
cell_ann <- cell_ann[!cell_ann$characteristics_ch1.5 %in% c("", "diagnosis: Other"), ] # remove 27 other and 1 unknown
exp_mat <- exp_mat[rownames(exp_mat) %in% cell_ann$geo_accession, ]

gene_ann <- data.frame(rowData(dat))
gene_ann <- gene_ann[gene_ann$gene_type == "protein_coding", ]
gene_ann <- gene_ann[!duplicated(gene_ann$gene_name), ]
int <- intersect(gene_ann$gene_id, colnames(exp_mat))
exp_mat <- exp_mat[, colnames(exp_mat) %in% int]
gene_ann <- gene_ann[gene_ann$gene_id %in% int, ]
colnames(exp_mat) <- gene_ann$gene_name

dat <- list(marray_exp = exp_mat, marray_gene_ann = gene_ann, 
            cell_ann_marray = cell_ann)

qsave(dat, file=  file.path(dir_out, "GSE21050_sarc_data.qs"))
# 282 samples, 16,014 genes 

### GSE30929 (microarray, STS)

dat <- qread(file= file.path(dir_in, "GSE30929_RangedSummarizedExperiment_2025-07-28.qs")) # old version: 2022-06-21
exp_mat <- assay(dat)
exp_mat <- t(exp_mat)

cell_ann <- data.frame(colData(dat))
exp_mat <- exp_mat[rownames(exp_mat) %in% cell_ann$geo_accession, ]

gene_ann <- data.frame(rowData(dat))
gene_ann <- gene_ann[gene_ann$gene_type == "protein_coding", ]
gene_ann <- gene_ann[!duplicated(gene_ann$gene_name), ]
int <- intersect(gene_ann$gene_id, colnames(exp_mat))
exp_mat <- exp_mat[, colnames(exp_mat) %in% int]
gene_ann <- gene_ann[gene_ann$gene_id %in% int, ]
colnames(exp_mat) <- gene_ann$gene_name

dat <- list(marray_exp = exp_mat, marray_gene_ann = gene_ann, 
            cell_ann_marray = cell_ann)

qsave(dat, file=  file.path(dir_out, "GSE30929_sarc_data.qs"))
# 140 samples, 11,245 genes 

######################################################################################################
## PGx + GEO (RNA) datasets 
## needs to be checked! Remove normal tissue types + other types + unknowns + MFS
######################################################################################################
# GSE21122
dat <- qread(file = file.path(dir_out, "GSE21122_sarc_data.qs"))
dat_ann <- dat$cell_ann_marray
dat_ann$mod_characteristics_ch1.1 <- substr(dat_ann$characteristics_ch1.1, 17, 
                                            nchar(dat_ann$characteristics_ch1.1))
subtype <- sapply(1:nrow(dat_ann), function(i){
  
  if(dat_ann$mod_characteristics_ch1.1[i] == "Leiomyosarcoma:N/A"){ df <- "Leiomyosarcoma" }
  
  if(dat_ann$mod_characteristics_ch1.1[i] %in% c("Liposarcoma:Dedifferentiated",
                                              "Liposarcoma:Myxoid/RC",
                                              "Liposarcoma:Pleomorphic") ){ df <- "Liposarcoma" }
  
  if(dat_ann$mod_characteristics_ch1.1[i] %in% c("MFH:Myxofibrosarcoma",
                                              "MFH:Pleomorphic") ){ df <- "MFH" }
  
  df
  
})

dat_ann$subtype <- subtype
gse21122_ann <- data.frame(sampleID = dat_ann$geo_accession,
                           lineage = "Soft Tissue",
                           subtype = dat_ann$subtype,
                           subtype_original = dat_ann$characteristics_ch1.1,
                           Primary.Metastasis = "NA",
                           type = "tumor",
                           Type = "GSE21122" )

gse21122_mat <- dat$marray_exp
gse21122_mat <- gse21122_mat[rownames(gse21122_mat) %in% gse21122_ann$sampleID, ]

# GSE21050
dat <- qread(file = file.path(dir_out, "GSE21050_sarc_data.qs"))
dat_ann <- dat$cell_ann_marray
dat_ann$mod_characteristics_ch1.5 <- substr(dat_ann$characteristics_ch1.5, 12, 
                                            nchar(dat_ann$characteristics_ch1.5))

subtype <- sapply(1:nrow(dat_ann), function(i){
  
  if(dat_ann$mod_characteristics_ch1.5[i] == "Leiomyosarcoma"){ df <- "Leiomyosarcoma" }
  
  if(dat_ann$mod_characteristics_ch1.5[i] == c("Liposarcoma - dedifferentiated") ){ df <- "Liposarcoma" }
  
  if(dat_ann$mod_characteristics_ch1.5[i] == c("Undifferentiated sarcoma") ){ df <- "UPS" }
  
  df
  
})
dat_ann$subtype <- subtype

gse21050_ann <- data.frame(sampleID = dat_ann$geo_accession,
                           lineage = "Soft Tissue",
                           subtype = dat_ann$subtype,
                           subtype_original = dat_ann$characteristics_ch1.5,
                           Primary.Metastasis = "NA",
                           type = "tumor",
                           Type = "GSE21050")

gse21050_ann$Type <- ifelse(dat_ann$characteristics_ch1 == "cohort: cohort 1",
                            paste("cohort 1", "GSE21050", sep="-"), 
                            paste("cohort 2", "GSE21050", sep="-"))

gse21050_mat <- dat$marray_exp
gse21050_mat <- gse21050_mat[rownames(gse21050_mat) %in% gse21050_ann$sampleID, ]

# GSE30929
dat <- qread(file = file.path(dir_out, "GSE30929_sarc_data.qs"))
dat_ann <- dat$cell_ann_marray

gse30929_ann <- data.frame(sampleID = dat_ann$geo_accession,
                           lineage = "Soft Tissue",
                           subtype = "Liposarcoma",
                           subtype_original = dat_ann$subtype.ch1,
                           Primary.Metastasis = "NA",
                           type = "tumor",
                           Type = "GSE30929" )

gse30929_mat <- dat$marray_exp
gse30929_mat <- gse30929_mat[rownames(gse30929_mat) %in% gse30929_ann$sampleID, ]

tcga_ann <- rbind(gse21122_ann, gse21050_ann, gse30929_ann)
int <- intersect(intersect(colnames(gse21050_mat), colnames(gse21122_mat)),
                 colnames(gse30929_mat)) # 11245

gene_ann <- dat$marray_gene_ann
gene_ann <- gene_ann[gene_ann$gene_type == "protein_coding", ]

int <- intersect(int, gene_ann$gene_name) # 11245 protein-coding genes and 571 samples

gse21122_mat <- gse21122_mat[, colnames(gse21122_mat) %in% int]
gse21050_mat <- gse21050_mat[, colnames(gse21050_mat) %in% int]
gse30929_mat <- gse30929_mat[, colnames(gse30929_mat) %in% int]

tcga_mat <- rbind(gse21122_mat,
                  gse21050_mat,
                  gse30929_mat)


#####################################
## PGx data
#####################################
dat <- qread(file.path(dir_out, "PGx_sarc_data.qs"))
rna_data <- dat$rna

ccle_ann <- rna_data$CCLE$cell_ann_marray
ccle_ann$sampleid <- paste(ccle_ann$sampleid, "ccle", sep = "-")
ccle_mat <- rna_data$CCLE$marray_exp
rownames(ccle_mat) <- ccle_ann$sampleid

gdsc_ann <- rna_data$GDSCv2$cell_ann_marray
gdsc_ann$sampleid <- paste(gdsc_ann$sampleid, "gdsc", sep="-")
gdsc_mat <- rna_data$GDSCv2$marray_exp
rownames(gdsc_mat) <- gdsc_ann$sampleid

nci_ann <- rna_data$NCI_Sarcoma$cell_ann_marray
nci_ann$sampleid <- paste(nci_ann$sampleid, "nci", sep ="-")
nci_mat <- rna_data$NCI_Sarcoma$marray_exp
rownames(nci_mat) <- nci_ann$sampleid  

pset_ann <- data.frame(sampleID = c(ccle_ann$sampleid,
                                    gdsc_ann$sampleid,
                                    nci_ann$sampleid),
                       lineage = c(ccle_ann$mod_tissueid,
                                   gdsc_ann$mod_tissueid,
                                   nci_ann$mod_tissueid),
                       subtype = c(ccle_ann$cellosaurus.cellosaurus.disease,
                                   gdsc_ann$cellosaurus.cellosaurus.disease,
                                   nci_ann$cellosaurus.cellosaurus.disease),
                       subtype_original = 'NA',
                       Primary.Metastasis = "NA",
                       type = "CL",
                       Type = c(rep("CCLE", length(ccle_ann$sampleid)),
                                rep("GDSC", length(gdsc_ann$sampleid)),
                                rep("NCI-Sarcoma", length(nci_ann$sampleid))) )

int <- intersect(intersect(colnames(gdsc_mat),
                           intersect(colnames(ccle_mat), colnames(nci_mat))), 
                           colnames(tcga_mat)) # 11049


gene_ann <- dat$rna$CCLE$marray_gene_ann
gene_ann_protein_coding <- gene_ann[gene_ann$GeneBioType == "protein_coding", ]
int <- intersect(int, gene_ann$Symbol) # 11049

tcga_mat <- tcga_mat[, int]
ccle_mat <- ccle_mat[, int]
gdsc_mat <- gdsc_mat[, int]
nci_mat <- nci_mat[, int]

pset_mat <- rbind(ccle_mat, gdsc_mat, nci_mat)

dat <- list(TCGA_mat = tcga_mat,
            TCGA_ann = tcga_ann,
            pset_ann = pset_ann,
            pset_mat = pset_mat)

qsave(dat, file= file.path(dir_out, "PGx_gse_rna.qs")) # 65 STS and 96 Bone

##################################################################
##  separate STS and Bone for alignment: Microarray 
##################################################################

dat <- qread(file.path(dir_out, "PGx_gse_rna.qs"))

pset_ann <- dat$pset_ann
pset_mat <- dat$pset_mat
TCGA_ann <- dat$TCGA_ann
TCGA_mat <- dat$TCGA_mat

pset_ann_sts <- pset_ann[pset_ann$lineage == "Soft Tissue", ]
pset_mat_sts <- pset_mat[rownames(pset_mat) %in% pset_ann_sts$sampleID, ]

TCGA_ann_sts <- TCGA_ann[TCGA_ann$lineage == "Soft Tissue", ]
TCGA_mat_sts <- TCGA_mat[rownames(TCGA_mat) %in% TCGA_ann_sts$sampleID, ]

dat <- list(TCGA_mat = TCGA_mat_sts,
            TCGA_ann = TCGA_ann_sts,
            pset_ann = pset_ann_sts,
            pset_mat = pset_mat_sts)

qsave(dat, file= file.path(dir_out, "PGx_gse_rna_sts.qs"))

# Note
# sub-type and original sub-type information for cell line data were added. 

