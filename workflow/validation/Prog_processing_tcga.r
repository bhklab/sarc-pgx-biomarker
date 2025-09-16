##-------------------------------------------------------------------------------
## Script: curate_TCGA_validation_data.R
##
## Description:
## This script processes and curates TCGA sarcoma RNA-seq and clinical data
## for downstream analysis. It performs the following steps:
## 1. Loading and cleaning clinical and gene expression data.
## 2. Filtering for relevant samples and protein-coding genes.
## 3. Creating and saving a curated SummarizedExperiment object.
##
## Output:
## - A serialized .qs file containing the curated SummarizedExperiment object.
##-------------------------------------------------------------------------------
########################################################
## load libraries
########################################################
library(qs)
library(SummarizedExperiment)

########################################################
## Set up directory
########################################################

dir_in <- 'data/rawdata/validation'
dir_out <- 'data/procdata/validation'

#######################################################
## Load clinical data
#######################################################
clinical_data <- read.csv(file.path(dir_in, 'TCGA/TCGA_clinical', 'clinical_curatedTCGA.csv'))

keyVar <- c('patientID', 'gender', 'race', 'ethnicity',
            'years_to_birth', 'days_to_birth', 'histological_type',
            'tumor_tissue_site', 'vital_status', 'days_to_death',
            'days_to_last_followup', 'residual_tumor', 
            'patient.clinical_cqcf.tumor_type', 'patient.clinical_cqcf.consent_or_death_status',
            'patient.clinical_cqcf.anatomic_neoplasm_subdivision',
            'date_of_initial_pathologic_diagnosis')

clin <- clinical_data[, colnames(clinical_data) %in% keyVar]

## standardize clinical variables
clin$age <- as.numeric(clin$years_to_birth)
clin$os_time_days <- ifelse(
  clin$vital_status == 1,
  clin$days_to_death,
  clin$days_to_last_followup
)

clin$os_time <- clin$os_time_days / 30.44
clin$os_event <- clin$vital_status
colnames(clin)[1] <- 'patientid' 
clin <- clin[order(clin$patientid), ]

clin$histological <- ifelse(clin$histological_type %in% c("undifferentiated pleomorphic sarcoma (ups)",
                                                          "pleomorphic 'mfh' / undifferentiated pleomorphic sarcoma",
                                                          "giant cell 'mfh' / undifferentiated pleomorphic sarcoma with giant cells"),
                                                          "UPS-MFH",
                                                          clin$histological_type)
clin$histological <- ifelse(clin$histological %in% c("malignant peripheral nerve sheath tumors (mpnst)",
                                                          "sarcoma; synovial; poorly differentiated",
                                                          "synovial sarcoma - biphasic",
                                                          "synovial sarcoma - monophasic"),
                                                          "synovial-mpnst",
                                                          clin$histological)
clin$histological <- ifelse(clin$histological == "dedifferentiated liposarcoma", "DDLPS", clin$histological)
clin$histological <- ifelse(clin$histological == "leiomyosarcoma (lms)", "LMS", clin$histological)
clin$histological <- ifelse(clin$histological == "myxofibrosarcoma", "MFS", clin$histological)
clin <- clin[clin$histological != "desmoid tumor", ]

########################################################
## Load TCGA RNA-seq
########################################################
dat <- readRDS(file = file.path(dir_in, 'TCGA', 'TCGA_RNA_seq.rds'))
dat_ann <- readRDS( file = file.path(dir_in, 'TCGA', 'TCGA_all_sample_annot.rds'))
dat_ann <- dat_ann[dat_ann$lineage == "SARC", ]

dat <- dat[rownames(dat) %in% dat_ann$sampleID, ]
expr <- log2(dat + 1)

patientid <- sapply(1:nrow(expr), function(k){
    id <- strsplit(rownames(expr)[k], "-")[[1]] 
    paste(id[1], id[2], id[3], sep="-")
})

rownames(expr) <- patientid
expr <- t(expr)
expr <- expr[grep("^ENSG", rownames(expr)), ]
expr <- expr[order(rownames(expr)),]
expr <- expr[, order(colnames(expr))]
expr <- expr[, !duplicated(colnames(expr))]

## gene annotation file 
load(file.path(dir_in, 'Gencode.v40.annotation.RData'))
gene_ann <- features_gene
gene_ann_pc <- gene_ann[gene_ann$gene_type == 'protein_coding', ]
gene_ann_pc <- gene_ann_pc[!duplicated(gene_ann_pc$gene_name), ]
gene_ann_pc <- gene_ann_pc[order(gene_ann_pc$gene_id), ]
gene_id <- substr(rownames(gene_ann_pc), 1, 15)
rownames(gene_ann_pc) <- gene_id

## common genes
int <- intersect(rownames(expr), rownames(gene_ann_pc)) # 16522
expr <- expr[rownames(expr) %in% int, ]
gene_ann_pc <- gene_ann_pc[rownames(gene_ann_pc) %in% int, ]
rownames(expr) <- gene_ann_pc$gene_name

## common samples
int <- intersect(colnames(expr), clin$patientid) # 257 (remove two patinets with desmoid tumor)
expr <- expr[, colnames(expr) %in% int ]
clin <- clin[clin$patientid %in% int, ] 
rownames(clin) <- clin$patientid
rownames(gene_ann_pc ) <- gene_ann_pc$gene_name

## build SE object
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(expr)), 
  colData = clin,                          
  rowData = gene_ann_pc                    
)

qsave(se, file = file.path(dir_out, 'curated_TCGA_se.qs'))
