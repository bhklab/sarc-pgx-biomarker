# -----------------------------------------------------------
# Soft-Tissue Sarcoma Data Visualization Script
# This script visualizes curated gene expression and drug 
# response data for soft-tissue sarcoma from:
#
#   - GEO (clinical gene expression datasets)
#   - ORCESTRA (pharmacogenomics cell line data)
#
#   Includes:
#     * Heatmap of cell line data distribution
#     * Heatmap of Cell line overlap between data sources 
#     * Table of cell line distribution
#
#   Assumes input data has already been processed
# -----------------------------------------------------------
###################################################################
## Load libraries
###################################################################

library(qs)
library(paletteer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(grid)
library(colorRamp2)
library(reshape2)
library(ComplexHeatmap)

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/procdata'
dir_out <- 'data/results' 
  
##################################################################
# Load curated PGx RNA cell line 
##################################################################
dat <- qread(file.path(dir_in, "PGx_gse_rna_sts.qs"))
dat.annot <- dat$pset_ann

#---- visualizae cell line distribution across PGx RNA data
df <- dat.annot
group <- unique(df$Type)

annot <- lapply(1:length(group), function(k){
  
  sub.df <- df[df$Type == group[k], ]
  
  if(group[k] != 'NCI-Sarcoma'){
    sub.df$sampleID <- substr(sub.df$sampleID, 1, nchar(sub.df$sampleID) - 5)
  }else{
    sub.df$sampleID <- substr(sub.df$sampleID, 1, nchar(sub.df$sampleID) - 4)
      }

  sub.df

  })

annot <- do.call(rbind, annot)

#---- Generate STable 1 & 2 only PGx RNA

cellID <- unique(annot$sampleID) # --> 39 cells
df.table <- lapply(1:length(cellID), function(j){
  
 sub.annot <- annot[annot$sampleID == cellID[j], ] 
 
 data.frame(CellLine = cellID[j],
            OncoTree1 = 'Soft Tissue',
            OncoTree2 = unique(sub.annot$subtype),
            CCLE = sum(sub.annot$Type == "CCLE"),
            GDSC = sum(sub.annot$Type == "GDSC"),
            'NCI-Sarcoma' = sum(sub.annot$Type == "NCI-Sarcoma"))
  
 })

df.table <- do.call(rbind, df.table)
write.csv(df.table, file = file.path(dir_out, 'data', 'cellLineInfo_RNA.csv'), row.names = FALSE)

# Note: manually update the subtypes ---> OncoType level 2

#---- Generate Figure 1 - RNA (heatmap) 

df <- read.csv(file.path(dir_out, 'data', 'STable_PGx Sarcoma_RNA - STable2.csv'))
rownames(df) <- df$CellLine
df <- df[, -c(1,2,3)]
colnames(df)[1] <- 'OncoTree'
colnames(df)[4] <- 'NCI-Sarcoma'

df$OncoTree <- ifelse(df$OncoTree %in% c('Alveolar Soft Part Sarcoma',
                                                'Endometrioid Stromal Sarcoma',
                                                'Epithelioid Sarcoma',
                                                'Malignant Peripheral Nerve Sheath Tumor',
                                                'Spindle Cell Sarcoma', 'Sarcoma'), 'Other', df$OncoTree)

# table(df$OncoTree2)[order(table(df$OncoTree2), decreasing = TRUE)]
custom_order <- c("Rhabdomyosarcoma", 
                  "Synovial Sarcoma", 
                  "Fibrosarcoma", 
                  "Leiomyosarcoma", 
                  "Liposarcoma", 
                  "Chondrosarcoma",
                  "Pleomorphic Sarcoma",
                  "Uterine Corpus Sarcoma",
                  "Other")
df$OncoTree <- factor(df$OncoTree, levels = custom_order)

# Create annotation data frame
anno.df <- data.frame(OncoTree = df$OncoTree)
rownames(anno.df) <- rownames(df)

# Define annotation

colors_vec <- c(
  "#4477AA", # blue
  "#EE6677", # coral
  "#228833", # green
  "#CCBB44", # mustard
  "#66CCEE", # sky blue
  "#AA3377", # purple
  "#DDCC77", # tan
  "#994455", # dark rose
  "#D3D3D3"  # forest green
)

oncotree_colors <- setNames(as.character(colors_vec), levels(df$OncoTree))
  
#--- Create annotation object with colors mapped to levels

ha <- HeatmapAnnotation(
  df = anno.df,
  col = list(OncoTree = oncotree_colors),
  show_annotation_name = TRUE
)

#--- Build matrix with studies as rows, cell lines as columns
mat0 <- as.matrix(df[, c("CCLE","GDSC","NCI-Sarcoma")])
# Ensure discrete mapping for "0"/"1"
storage.mode(mat0) <- "character"
mat <- t(mat0)  # now rows = studies, cols = cell lines

# Order columns by OncoTree
ord <- order(df$OncoTree)
mat <- mat[, ord, drop = FALSE]
anno.ord <- anno.df[ord,, drop = FALSE]
ha <- HeatmapAnnotation(
  df = anno.ord,
  col = list(OncoTree = oncotree_colors),
  show_annotation_name = TRUE
)

# Create binary heatmap
pdf(file.path(dir_out, 'data', "Fig1_cellInfo_RNA.pdf"), width = 12, height = 2.5)

Heatmap(mat,
        name = "Presence",
        col = c("0" = "white", "1" = "#7C7C7CFF"),
        show_column_names = TRUE,
        show_row_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        #right_annotation = ha,
        top_annotation = ha,
        row_title = " ",
        column_title = " ",
        #row_order = order(df$OncoTree),
        column_order = seq_len(ncol(mat)), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y,
                    width = width,
                    height = height,
                    gp = gpar(col = "#242424", fill = NA, lwd = 0.8))
        })

dev.off()

# save ordered and plotted the RNA heatmap (mat is studies x cell lines)
saveRDS(colnames(mat), file.path(dir_out, "data", "cellline_order_Fig1_RNA.rds"))

#---- Generate Figure 1: cell line overlap between RNA data sources 

df <- read.csv(file.path(dir_out, 'data', 'STable_PGx Sarcoma_RNA - STable2.csv'))
rownames(df) <- df$CellLine
df <- df[, -c(2,3)]
colnames(df)[5] <- 'NCI-Sarcoma'

# Example binary presence matrix
mat <- df %>%
  select(CCLE, GDSC, 'NCI-Sarcoma') %>%
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

pdf(file.path(dir_out, 'data', "Fig1_cellInfo_acrossstudies_RNA.pdf"), 
    width = 2, height = 2)

ggplot(df_heat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +  
  geom_text(aes(label = value, color = value), size = 4) +  
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

##################################################################
# Load curated PGx drug response cell line 
##################################################################
dat.drug <- qread(file.path(dir_in, "PGx_sarc_data.qs"))
dat.drug.ctrp <-  dat.drug$aac$CTRP$cell_ann_seq
dat.drug.ctrp <- dat.drug.ctrp[dat.drug.ctrp$tissueid != 'Bone', ]
dat.annot.ctrp <- data.frame(sampleID = paste(dat.drug.ctrp$sampleid, 'ctrp', sep='-'),
                            lineage = 'Soft Tissue',
                            subtype = dat.drug.ctrp$cellosaurus.cellosaurus.disease,
                            subtype_original = 'NA',
                            Primary.Metastasis = 'NA',
                            type = 'CL',
                            Type = 'CTRP')


dat.drug.nci <-  dat.drug$aac$NCI$cell_ann_seq
dat.drug.nci <- dat.drug.nci[dat.drug.nci$dimitrios.soft_vs_bone %in% c("Soft", "Mixed"), ]
dat.annot.nci <- data.frame(sampleID = paste(dat.drug.nci$sampleid, 'nci', sep='-'),
                            lineage = 'Soft Tissue',
                            subtype = dat.drug.nci$cellosaurus.cellosaurus.disease,
                            subtype_original = 'NA',
                            Primary.Metastasis = 'NA',
                            type = 'CL',
                            Type = 'NCI-Sarcoma')

dat.drug.gdsc <-  dat.drug$aac$GDSCv2$cell_ann_seq
dat.drug.gdsc <- dat.drug.gdsc[dat.drug.gdsc$tissueid != 'Bone', ]
dat.annot.gdsc <- data.frame(sampleID = paste(dat.drug.gdsc$sampleid, 'gdsc', sep='-'),
                            lineage = 'Soft Tissue',
                            subtype = dat.drug.gdsc$cellosaurus.cellosaurus.disease,
                            subtype_original = 'NA',
                            Primary.Metastasis = 'NA',
                            type = 'CL',
                            Type = 'GDSC')

#---- visualizae cell line distribution across PGx drug response data

df <- rbind(dat.annot.ctrp, dat.annot.gdsc, dat.annot.nci)
group <- unique(df$Type)

annot <- lapply(1:length(group), function(k){
  
  sub.df <- df[df$Type == group[k], ]
  
  if(group[k] != 'NCI-Sarcoma'){
    sub.df$sampleID <- substr(sub.df$sampleID, 1, nchar(sub.df$sampleID) - 5)
  }else{
    sub.df$sampleID <- substr(sub.df$sampleID, 1, nchar(sub.df$sampleID) - 4)
      }

  sub.df

  })

annot <- do.call(rbind, annot)

#---- Generate STable 1 & 2 only PGx drug-response

cellID <- unique(annot$sampleID)
df.table <- lapply(1:length(cellID), function(j){
  
 sub.annot <- annot[annot$sampleID == cellID[j], ] 
 
 data.frame(CellLine = cellID[j],
            OncoTree1 = 'Soft Tissue',
            OncoTree2 = unique(sub.annot$subtype),
            CTRP = sum(sub.annot$Type == "CTRP"),
            GDSC = sum(sub.annot$Type == "GDSC"),
            'NCI-Sarcoma' = sum(sub.annot$Type == "NCI-Sarcoma"))
  
 })

df.table <- do.call(rbind, df.table)
write.csv(df.table, file = file.path(dir_out, 'data', 'cellLineInfo_drug.csv'), row.names = FALSE)

# Note: manually update the subtypes ---> OncoType level 2

#---- Generate Figure 1 - drug response (heatmap) 

df <- read.csv(file.path(dir_out, 'data', 'STable_PGx Sarcoma_drug - STable2.csv'))
rownames(df) <- df$CellLine
df <- df[, -c(1,2,3)]
colnames(df)[1] <- 'OncoTree'
colnames(df)[4] <- 'NCI-Sarcoma'

df$OncoTree <- ifelse(df$OncoTree %in% c('Alveolar Soft Part Sarcoma',
                                                'Endometrioid Stromal Sarcoma',
                                                'Epithelioid Sarcoma',
                                                'Malignant Peripheral Nerve Sheath Tumor',
                                                'Spindle Cell Sarcoma', 'Sarcoma'), 'Other', df$OncoTree)

# table(df$OncoTree2)[order(table(df$OncoTree2), decreasing = TRUE)]
custom_order <- c("Rhabdomyosarcoma", 
                  "Synovial Sarcoma", 
                  "Fibrosarcoma", 
                  "Leiomyosarcoma", 
                  "Liposarcoma", 
                  "Chondrosarcoma",
                  "Pleomorphic Sarcoma",
                  "Uterine Corpus Sarcoma",
                  "Other")
df$OncoTree <- factor(df$OncoTree, levels = custom_order)

# Create annotation data frame
anno.df <- data.frame(OncoTree = df$OncoTree)
rownames(anno.df) <- rownames(df)

# Define annotation

colors_vec <- c(
  "#4477AA", # blue
  "#EE6677", # coral
  "#228833", # green
  "#CCBB44", # mustard
  "#66CCEE", # sky blue
  "#AA3377", # purple
  "#DDCC77", # tan
  "#994455", # dark rose
  "#D3D3D3"  # forest green
)

oncotree_colors <- setNames(as.character(colors_vec), levels(df$OncoTree))
  
#--- Create annotation object with colors mapped to levels

ha <- HeatmapAnnotation(
  df = anno.df,
  col = list(OncoTree = oncotree_colors),
  show_annotation_name = TRUE
)

#--- Build matrix with studies as rows, cell lines as columns
mat0 <- as.matrix(df[, c("CTRP","GDSC","NCI-Sarcoma")])
# Ensure discrete mapping for "0"/"1"
storage.mode(mat0) <- "character"
mat <- t(mat0)  # now rows = studies, cols = cell lines

#--- Reuse column order from RNA heatmap (if available)
order_path <- file.path(dir_out, "data", "cellline_order_Fig1_RNA.rds")
ref_order <- if (file.exists(order_path)) readRDS(order_path) else NULL

if (!is.null(ref_order)) {
  shared <- intersect(ref_order, colnames(mat))                  
  extra  <- setdiff(colnames(mat), ref_order)                    
  extra_ord <- extra[order(df[extra, "OncoTree"], extra)]
  final_order <- c(shared, extra_ord)
  mat <- mat[, final_order, drop = FALSE]
  anno.df <- anno.df[final_order, , drop = FALSE]
} else {
  final_order <- colnames(mat)[order(df[colnames(mat), "OncoTree"], colnames(mat))]
  mat <- mat[, final_order, drop = FALSE]
  anno.df <- anno.df[final_order, , drop = FALSE]
}

# Order columns by OncoTree
ha <- HeatmapAnnotation(
  df = anno.ord,
  col = list(OncoTree = oncotree_colors),
  show_annotation_name = TRUE
)

# Create binary heatmap
pdf(file.path(dir_out, 'data', "Fig1_cellInfo_Drug.pdf"), width = 12, height = 2.5)

Heatmap(mat,
        name = "Presence",
        col = c("0" = "white", "1" = "#7C7C7CFF"),
        show_column_names = TRUE,
        show_row_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        #right_annotation = ha,
        top_annotation = ha,
        row_title = " ",
        column_title = " ",
        #row_order = order(df$OncoTree),
        column_order = seq_len(ncol(mat)), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y,
                    width = width,
                    height = height,
                    gp = gpar(col = "#242424", fill = NA, lwd = 0.8))
        })

dev.off()

#---- Generate Figure 1: cell line overlap between RNA data sources 

df <- read.csv(file.path(dir_out, 'data', 'STable_PGx Sarcoma_drug - STable2.csv'))
rownames(df) <- df$CellLine
df <- df[, -c(2,3)]
colnames(df)[5] <- 'NCI-Sarcoma'

# Example binary presence matrix
mat <- df %>%
  select(CTRP, GDSC, 'NCI-Sarcoma') %>%
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

pdf(file.path(dir_out, 'data', "Fig1_cellInfo_acrossstudies_drug.pdf"), 
    width = 2, height = 2)

ggplot(df_heat, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +  
  geom_text(aes(label = value, color = value), size = 4) +  
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

##################################################################
# Load curated GEO data
##################################################################
dat <- qread(file.path(dir_in, "PGx_gse_rna_sts.qs"))
dat.annot <- dat$TCGA_ann

dat.annot$Type <- ifelse(dat.annot$Type %in% c('cohort 1-GSE21050', 'cohort 2-GSE21050'), 
                        'GSE21050', dat.annot$Type)

# Summarize counts: Subtype vs Study (Type)
mat_df <- dat.annot %>%
  dplyr::count(subtype, Type) %>%
  pivot_wider(names_from = Type, values_from = n, values_fill = 0)

# Convert to matrix with rownames
mat <- as.data.frame(mat_df)
rownames(mat) <- mat$subtype
mat$subtype <- NULL
mat <- as.matrix(mat)

max_val <- max(mat)
col_fun <- colorRamp2(
  breaks = c(0, max_val / 2, max_val),
  colors = c("#ffffff", "#dddddd", "#888888") # very light to medium blue
)

pdf(file.path(dir_out, 'data', "Fig1_geo.pdf"), width = 4.5, height = 3.5)

Heatmap(mat,
        name = " ",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          # Add text
          grid.text(mat[i, j], x, y, gp = gpar(fontsize = 12))
          # Add border
          grid.rect(x, y, width, height, gp = gpar(col = "#242424", fill = NA, lwd = 0.5))
        },
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))

dev.off()

