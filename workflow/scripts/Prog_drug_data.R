# -----------------------------------------------------------
# Drug Response Integration & Visualization for STS Cell Lines
#
# This script processes and visualizes drug response data 
# across multiple pharmacogenomic datasets for soft-tissue sarcoma:
#
#   - CTRP (via CCLE annotation)
#   - GDSCv2
#   - NCI-Sarcoma
#
#   Includes:
#     * Filtering and QC of drug response profiles
#     * Venn diagram of drug overlaps across datasets
#     * Heatmap of drug overlap across studies
#     * Table summarizing drug availability
#     * Alignment of drug response and gene expression profiles
#
#   Final output is a curated and aligned dataset of drug 
#   response (AAC) and gene expression profiles for overlapping 
#   STS cell lines, saved for downstream analyses.
#
#   Assumes input data (qs files) has been preprocessed
#   and available under specified `data/` directories.
# -----------------------------------------------------------
########################################
## Load libraries
########################################

library(qs)
library(ggplot2)
library(reshape2)
library(dplyr)
library(PharmacoGx)
library(VennDiagram)

##################################################################
## Setup directory
##################################################################

dir_data <- 'data/procdata' # '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/data'
dir_aligned <- 'data/results/aligned' # '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/aligned/'
dir_out <-  'data/results/drug' # '/home/bioinf/bhklab/farnoosh/STS-PGx-Biomarker/result/drug/'

################################################################################
## Load drug response data
################################################################################
dat <- qread(file= file.path(dir_data, "PGx_sarc_data.qs"))

gene_ann <- dat$rna$CCLE$marray_gene_ann

ccle_aac <- dat$aac$CCLE$drug_aac
keep <- dat$aac$CCLE$cell_ann_seq[dat$aac$CCLE$cell_ann_seq$tissueid != 'Bone', 'sampleid']
ccle_aac <- ccle_aac[, keep]
remove <- sapply(1:nrow(ccle_aac), function(k){
  sum(is.na(ccle_aac[k, ]))/ncol(ccle_aac)
})
ccle_aac <- ccle_aac[- which(remove > 0.60), ] # 22 out of 24 

ctrp_aac <- dat$aac$CTRP$drug_aac
keep <- dat$aac$CTRP$cell_ann_seq[dat$aac$CTRP$cell_ann_seq$tissueid != 'Bone', 'sampleid']
ctrp_aac <- ctrp_aac[, keep]
remove <- sapply(1:nrow(ctrp_aac), function(k){
  sum(is.na(ctrp_aac[k, ]))/ncol(ctrp_aac)
})
ctrp_aac <- ctrp_aac[- which(remove > 0.60), ] # 475 out of 544 

gdsc_aac <- dat$aac$GDSCv2$drug_aac
keep <- dat$aac$GDSCv2$cell_ann_seq[dat$aac$GDSCv2$cell_ann_seq$tissueid != 'Bone', 'sampleid']
gdsc_aac <- gdsc_aac[, keep]
remove <- sapply(1:nrow(gdsc_aac), function(k){
  sum(is.na(gdsc_aac[k, ]))/ncol(gdsc_aac)
})
gdsc_aac <- gdsc_aac[- which(remove > 0.60), ] # 309 out of 343

nci_aac <- dat$aac$NCI$drug_aac
keep <- dat$aac$NCI$cell_ann_seq[!is.na(dat$aac$NCI$cell_ann_seq$dimitrios.soft_vs_bone), ]
keep <- keep[keep$dimitrios.soft_vs_bone %in% c("Soft", "Mixed"), 'sampleid']
nci_aac <- nci_aac[, keep] # 437

## Venn-diagram
n1 <- length(rownames(ctrp_aac))
n2 <- length(rownames(gdsc_aac))
n3 <- length(rownames(nci_aac))
n12 <- length(intersect(rownames(ctrp_aac), rownames(gdsc_aac)))
n13 <- length(intersect(rownames(ctrp_aac), rownames(nci_aac)))
n23 <- length(intersect(rownames(gdsc_aac), rownames(nci_aac)))
n123 <- length(intersect(intersect(rownames(ctrp_aac), rownames(gdsc_aac)), rownames(nci_aac)))

pdf(file = file.path(dir_out, "venn_drug.pdf"), width = 5.5, height = 5.5)

venn.plot <- draw.triple.venn(
  area1 = n1,
  area2 = n2,
  area3 = n3,
  n12 = n12,
  n13 = n13,
  n23 = n23,
  n123 = n123,
  category = c("CCLE/CTRP", "GDSC", "NCI-Sarcoma"),
  fill = c("#64894DFF", "#99B6BDFF", "#ECC9A0FF"),
  cat.col = c("#252525", "#252525", "#252525"),
  #lty = "dashed",
  cex = 2,
  cat.cex = 1
)

dev.off()

##################################################################
## Drugs overlap between data sources
##################################################################

df <- data.frame(drug = c(rownames(ccle_aac), rownames(ctrp_aac), rownames(gdsc_aac), rownames(nci_aac)),
                 group = c(rep("CCLE", nrow(ccle_aac)),
                           rep("CTRP", nrow(ctrp_aac)), 
                           rep("GDSC", nrow(gdsc_aac)),
                           rep("NCI-Sarcoma", nrow(nci_aac))) )

drug <- unique(df$drug)

df.table <- lapply(1:length(drug), function(j){
  
  sub.df <- df[df$drug == drug[j], ] 

  data.frame(drug = drug[j],
             CCLE = sum(sub.df$group == "CCLE"),
             CTRP = sum(sub.df$group == "CTRP"),
             GDSC = sum(sub.df$group == "GDSC"),
             'NCI-Sarcoma' = sum(sub.df$group == "NCI-Sarcoma"))
  
})

df.table <- do.call(rbind, df.table)
write.csv(df.table, file = file.path(dir_out, 'STable4_drugInfo.csv'), row.names = FALSE)

# Example binary presence matrix
colnames(df.table)[5] <- 'NCI-Sarcoma'
mat <- df.table %>%
  select(CCLE, CTRP, GDSC, 'NCI-Sarcoma') %>%
  mutate_all(as.integer)

# Function to compute pairwise overlap counts
overlap_matrix <- function(x) {
  m <- as.matrix(x)
  n <- ncol(m)
  result <- matrix(0, nrow = n, ncol = n)
  colnames(result) <- rownames(result) <- colnames(m)
  
  for (i in 1:n) {
    for (j in 1:n) {
      result[i, j] <- sum(m[, i] & m[, j])
    }
  }
  return(result)
}

overlap_mat <- overlap_matrix(mat)

# View full symmetric matrix
print(overlap_mat)

# Show upper triangle only
overlap_upper <- overlap_mat
overlap_upper[lower.tri(overlap_upper)] <- NA

df_heat <- melt(overlap_upper, na.rm = TRUE)

pdf(file.path(dir_out, "Fig1_drug_acrossstudies.pdf"), 
    width = 4, height = 4)

ggplot(df_heat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +  
  geom_text(aes(label = value, color = value), size = 5) +  
  scale_color_gradient(low = "white", high = "black") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),            
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),     
    axis.text.y = element_text(size = 6),     
    axis.ticks = element_blank(),             
    axis.line = element_blank(),             
    panel.grid = element_blank(),             
    panel.background = element_blank(),       
    plot.background = element_blank(),       
    legend.position = "none" )

dev.off()

####################################################################
# Filter to common drugs across three PSets
####################################################################
int <- intersect(intersect(rownames(ctrp_aac), rownames(gdsc_aac)),
                 rownames(nci_aac)) # 69 drugs

ctrp_aac <- ctrp_aac[int, ]
colnames(ctrp_aac) <- paste(colnames(ctrp_aac), "ccle", sep="-")
nci_aac <- nci_aac[int, ]
colnames(nci_aac) <- paste(colnames(nci_aac), "nci", sep="-")
gdsc_aac <- gdsc_aac[int, ]
colnames(gdsc_aac) <- paste(colnames(gdsc_aac), "gdsc", sep="-")

pset_aac <- cbind(ctrp_aac, gdsc_aac, nci_aac)

## load aligned data 
dat <- qread(file.path(dir_aligned, "pgx_rna_celligner_sts.qs"))

dat_aligned <- as.data.frame(Seurat::GetAssayData(dat$comb_obj))
dat_ann <- data.frame(sampleID = dat$comb_obj$sampleID,
                      type = dat$comb_obj$type)

cl_id <- dat_ann[dat_ann$type == "CL", "sampleID"]
pset_aligned <- dat_aligned[ , cl_id]

int <- intersect(colnames(pset_aligned), colnames(pset_aac)) # 62 cell lines in common
pset_aac <- pset_aac[, colnames(pset_aac) %in% int] 
pset_aligned <- pset_aligned[, colnames(pset_aligned) %in% int]

dat <- list(pset_aac = pset_aac, pset_aligned = pset_aligned, gene_ann = gene_ann)
qsave(dat, file= file.path(dir_out, "drug_rna_aligned_sts.qs"))

