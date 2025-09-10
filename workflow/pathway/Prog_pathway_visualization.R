# -----------------------------------------------------------
# STS Pathway–Drug Association Visualization Script
#
# This script visualizes results from pathway enrichment analyses 
# (e.g., HALLMARK gene sets) for drug associations in soft tissue sarcoma (STS) 
# cell lines, using transcriptomic and pharmacogenomic data.
#
# Analyses included:
#   - GSEA (Gene Set Enrichment Analysis) results per drug
#   - ORA (Over-Representation Analysis) results per drug
#
# Key outputs:
#   1. Bar plots summarizing the number of significant pathway–drug associations
#   2. Dot plots for significant pathways by drug
#   3. UpSet plots showing pathway overlap across drugs
#   4. Heatmaps for pathway–drug significance and enrichment scores
# -----------------------------------------------------------
#############################################################
## Load libraries
#############################################################

library(ggplot2)
library(qs)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(VennDiagram)
library(grid)
library(UpSetR)
library(forcats)

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/results/pathway' 
dir_out <- 'data/results/pathway/figure'
dir_moa <- 'data/results/moa'
dir_drug <- 'data/procdata'

top_cutoff <- 30

#################################################################
## load selected drugs
#################################################################
selected_drug <- read.csv(file.path('data/results/drug', 'selected_drugs.csv'))
selected_drug_clin <- selected_drug[selected_drug$sts == 'Yes', 'drug']

# drug → class mapping
load(file.path('data/results/MOA', 'moa_info.rda'))
drug_class <- lapply(1:length(dat_drug_class), function(k){

 data.frame(moa = names(dat_drug_class)[k],
            drug = dat_drug_class[[k]])

})

drug_class <- do.call(rbind, drug_class)

##################################################################################################################
################################################### HALLMARK (GSEA) ##############################################
##################################################################################################################
#################################################################
## HALLMARK pathway significant results (keep direction of NES)
#################################################################
fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
gsea_data <- fgseaRes[fgseaRes$drug %in% selected_drug_clin, ]
gsea_data$pathway <- tools::toTitleCase(tolower(substr(gsea_data$pathway, 10, nchar(gsea_data$pathway))))

thr_gsea <- 0.10

# compute median NES per pathway
pathway_order <- gsea_data %>%
  filter(padj <= thr_gsea) %>%
  group_by(pathway) %>%
  summarise(median_NES = median(NES, na.rm = TRUE)) %>%
  arrange(desc(median_NES)) %>%
  pull(pathway)

# order pathways by median NES (signed) across drugs
gsea_data_filtered <- gsea_data %>%
  filter(padj <= thr_gsea) %>%
  mutate(
    FDR_bin = cut(padj,
                  breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
                  labels = c("< 0.01", "0.01–0.05", "0.05–0.10", "> 0.10")),
    NES_sign = ifelse(NES >= 0, "Upregulated", "Downregulated"),
    Pathway = factor(pathway, levels = pathway_order)
  )


pdf(file=file.path(dir_out, "heatmap_dot_hallmark_gsea.pdf"), 
width = 12, height = 10)

ggplot(gsea_data_filtered, aes(x = drug, y = pathway)) +
  geom_point(aes(size = abs(NES), fill = FDR_bin), shape = 21, color = "black") +
  scale_size_continuous(name = "NES", range = c(2, 8)) +
  scale_fill_manual(
    name = "FDR",
    values = c(
      "< 0.01"     = "#29126dff",
      "0.01–0.05"  = "#c041c9ff",
      "0.05–0.10"  = "#db9273ff",
      "> 0.10"     = "grey90"
    )
  ) +
  facet_grid(. ~ NES_sign) +
  theme_minimal(base_size = 12) +
  theme(
  strip.text = element_text(face = "bold", size = 14),  # Larger facet titles
  axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
  axis.text.y = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  plot.title = element_text(size = 12),
  plot.subtitle = element_text(size = 10),
  panel.spacing.x = unit(1.5, "lines"),
  panel.grid.major = element_line(color = "#fbf8f7ff")
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "",
    y = ""
  )

dev.off()

############################################################
## HALLMARK pathway significant results (abs of NES)
############################################################
fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
gsea_data <- fgseaRes[fgseaRes$drug %in% selected_drug_clin, ]
gsea_data$pathway <- tools::toTitleCase(tolower(substr(gsea_data$pathway, 10, nchar(gsea_data$pathway))))

thr_gsea <- 0.10
top_n_pathways <- 20

# compute median NES per pathway
pathway_order <- gsea_data %>%
  filter(padj <= thr_gsea) %>%
  group_by(pathway) %>%
  summarise(median_abs_NES = median(abs(NES), na.rm = TRUE)) %>%
  arrange(desc(median_abs_NES)) %>%
  slice_head(n = top_n_pathways) %>%
  pull(pathway)

# order pathways by median NES (signed) across drugs
gsea_data_filtered <- gsea_data %>%
  filter(pathway %in% pathway_order, padj <= thr_gsea) %>%
  mutate(
    FDR_bin = cut(padj,
                  breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
                  labels = c("< 0.01", "0.01–0.05", "0.05–0.10", "> 0.10")),
    pathway = factor(pathway, levels = pathway_order)
  ) 

pdf(file=file.path(dir_out, "heatmap_dot_hallmark_gsea_abs_top20.pdf"), 
width = 10, height = 7)

ggplot(gsea_data_filtered, aes(x = drug, y = pathway)) +
  geom_point(aes(size = abs(NES), fill = FDR_bin), shape = 21, color = "black") +
  scale_size_continuous(name = "|NES|", range = c(2, 8)) +
  scale_fill_manual(
    name = "FDR",
    values = c(
      "< 0.01"     = "#29126dff",
      "0.01–0.05"  = "#c041c9ff",
      "0.05–0.10"  = "#db9273ff",
      "> 0.10"     = "grey90"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
  axis.text.y = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10),
  plot.title = element_text(size = 12),
  plot.subtitle = element_text(size = 10),
  panel.spacing.x = unit(1.5, "lines"),
  panel.grid.major = element_line(color = "#fbf8f7ff")
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "",
    y = ""
  ) 

dev.off()

############################################################
## HALLMARK GSEA and ORA along with MOA 
############################################################
fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
foraRes <- qread(file.path(dir_in, 'hallmark_ora_pathway_drug.qs'))

fgseaRes$pathway <- tools::toTitleCase(tolower(substr(fgseaRes$pathway, 10, nchar(fgseaRes$pathway))))
foraRes$pathway <- tools::toTitleCase(tolower(substr(foraRes$pathway, 10, nchar(foraRes$pathway))))

thr_gsea <- 0.05
thr_ora <- 0.15

# GSEA results: keep relevant columns
gsea_data <- fgseaRes %>%
  select(drug, pathway, NES, padj) %>%
  rename(NES_gsea = NES, padj_gsea = padj)

# ORA results: keep relevant columns
ora_data <- foraRes %>%
  select(drug, pathway, padj) %>%
  rename(padj_ora = padj)

# Merge
merged <- gsea_data %>%
  left_join(ora_data, by = c("drug", "pathway")) %>%
  mutate(
    gsea_sig = padj_gsea <= thr_gsea,
    ora_sig  = padj_ora <= thr_ora,
    both_sig = gsea_sig & ora_sig
  )

merged <- merged %>%
  left_join(drug_class, by = "drug")

class_summary <- merged %>%
  group_by(moa, pathway) %>%
  summarise(
    mean_NES = mean(NES_gsea, na.rm = TRUE),
    overlap_count = sum(both_sig, na.rm = TRUE),  
    .groups = "drop"
  ) %>%
  mutate(
    glyph = case_when(
      overlap_count >= 2 ~ "robust",     # 2 drugs in class
      overlap_count == 1 ~ "single",     # more than 2 drug in class
      TRUE ~ "none"
    )
  )

# Filter out pathways with zero hits across all MOAs
filtered_class_summary <- class_summary %>%
  group_by(pathway) %>%
  filter(sum(overlap_count) > 1) %>%
  ungroup()

pathway_order <- filtered_class_summary %>%
  filter(glyph == "robust") %>%
  count(pathway, name = "star_count") %>%
  arrange(desc(star_count)) %>%
  pull(pathway)

missing_pathways <- setdiff(unique(filtered_class_summary$pathway), pathway_order)
pathway_order <- c(pathway_order, sort(missing_pathways))

filtered_class_summary$pathway_order <- factor(filtered_class_summary$pathway,
                                         levels = pathway_order)

# Optional: alphabetically order MOA
filtered_class_summary$moa_order <- factor(filtered_class_summary$moa,
                                     levels = sort(unique(filtered_class_summary$moa)))

pdf(file=file.path(dir_out, "heatmap_hallmark_gsea_ora_moa.pdf"), 
width = 10, height = 7)

ggplot(filtered_class_summary, aes(x = moa_order, y = pathway_order, fill = mean_NES)) +
   geom_tile(color = "#686969ff") +
  # Single-drug overlap → circle
  geom_point(
    data = subset(filtered_class_summary, glyph == "single"),
    aes(x = moa_order, y = pathway_order),
    shape = 19, size = 1.5, color = "black"
  ) +
  # ≥2 drugs overlap → cross
  geom_point(
    data = subset(filtered_class_summary, glyph == "robust"),
    aes(x = moa_order, y = pathway_order),
    shape = 8, size = 1.5, color = "black", fill = "black"
  ) +
  scale_fill_gradient2(
  low = "#093202ff",
  mid = "white",
  high = "#472902ff",
  midpoint = 0,
  limits = c(-2.5, 2.5),   # adjust to your NES range
  oob = scales::squish,
  name = "avg NES"
 ) + 
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "",
    y = ""
  )

dev.off()

############################################################
## UpSet figure 
############################################################
fgseaRes <- qread(file.path(dir_in, "hallmark_gsea_pathway_drug.qs"))

# Format leadingEdge column (into character vectors)
fgseaRes$leadingEdge <- strsplit(fgseaRes$leadingEdge, ",")
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, trimws)

# Clean pathway names
fgseaRes$pathway <- tools::toTitleCase(tolower(substr(fgseaRes$pathway, 10, nchar(fgseaRes$pathway))))

# Keep only selected drugs
fgseaRes <- fgseaRes[fgseaRes$drug %in% selected_drug_clin, ]

# Keep only significant results
fgseaRes <- fgseaRes[fgseaRes$padj < thr_gsea, ]

#id <- "E2F_TARGETS"
#core_drugs <- c(
#  "Temozolomide", "Venetoclax", "Dabrafenib", "Tamoxifen", "Doxorubicin",
#  "Nintedanib", "Panobinostat", "Cabozantinib", "Nilotinib", "Crizotinib")

id <- "G2m_checkpoint"
core_drugs <- c(
  "Temozolomide", "Venetoclax", "Dabrafenib", "Tamoxifen", "Doxorubicin",
  "Nintedanib",  "Cabozantinib", "Nilotinib", "Crizotinib")

df <- fgseaRes[fgseaRes$pathway == id & fgseaRes$drug %in% core_drugs, ]

# Build named list
leading_edge_list <- setNames(df$leadingEdge, df$drug)

# Remove empty sets (if any)
leading_edge_list <- leading_edge_list[sapply(leading_edge_list, length) > 0]

# Rebuild comb_mat
m <- make_comb_mat(leading_edge_list)
m <- m[comb_size(m) > 2]
core_moa <- drug_class[drug_class$drug %in% core_drugs, ]
drug2moa_raw <- core_moa$moa
names(drug2moa_raw) <- core_moa$drug

drug2moa <- drug2moa_raw[set_name(m)]
moa_colors <-   c("Apoptosis regulation" = "#D2C396FF",
                  #"Cell cycle" = "#78847FFF" ,
                  #"Chromatin regulation" = "#ACC2CFFF",
                  "DNA replication" = "#855C75FF",
                  #"EGFR signaling" = "#819574ff", 
                  "ERK MAPK signaling" = "#9C6755FF",
                  #"Genome integrity" = "#877772ff",
                  #"IGF1R signaling" = "#8491BEFF",
                  #"Kinases" = "#8097e8ff",
                  #"Mitosis" = "#B1AF53FF", 
                  "Other" = "#d1d5dfff",
                  #"PI3K MTOR signaling" = "#64894DFF", 
                  "RTK signaling" = "#526A83FF") 

left_annotation <- rowAnnotation(
  MOA = drug2moa,
  col = list(MOA = moa_colors),
  show_annotation_name = FALSE
)

top_annotation <- upset_top_annotation(
  m,
  gp = gpar(fill = "#72874E"),
  height = unit(5, "cm"),
  axis_param = list(labels = FALSE),
  add_numbers = TRUE,
  numbers_gp = gpar(fontsize = 8, col = "black")
)

right_annotation <- upset_right_annotation(
  m,
  gp = gpar(fill = "#023743"),
  width = unit(3, "cm")
)

pdf(file=file.path(dir_out, "upset_hallmark_gsea_g2m.pdf"), 
width = 10, height = 4)

p <- UpSet(
  m,
  top_annotation = top_annotation,
  right_annotation = right_annotation,
  left_annotation = left_annotation,
  set_order = order(set_size(m), decreasing = TRUE),
  comb_order = order(comb_size(m), decreasing = TRUE),
  pt_size = unit(2, "mm"),
  lwd = 1.2
)

draw(p)

dev.off()


############################################################
## Heatmap figure: E2f_targets and/or G2m_checkpoint
############################################################
## pathway result
fgseaRes <- qread(file.path(dir_in, "hallmark_gsea_pathway_drug.qs"))

# Format leadingEdge column (into character vectors)
fgseaRes$leadingEdge <- strsplit(fgseaRes$leadingEdge, ",")
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, trimws)

# Clean pathway names
fgseaRes$pathway <- tools::toTitleCase(tolower(substr(fgseaRes$pathway, 10, nchar(fgseaRes$pathway))))

# Keep only selected drugs
fgseaRes <- fgseaRes[fgseaRes$drug %in% selected_drug_clin, ]

# Keep only significant results
fgseaRes <- fgseaRes[fgseaRes$padj < thr_gsea, ]

id <- "G2m_checkpoint"
core_drugs <- c(
  "Temozolomide", "Venetoclax", "Dabrafenib", "Tamoxifen", "Doxorubicin",
  "Nintedanib", "Panobinostat", "Cabozantinib", "Nilotinib", "Crizotinib"
)
fgsea_sub <- fgseaRes[fgseaRes$drug %in% core_drugs & fgseaRes$pathway == id, ]

# Build matrix of leading-edge genes (TRUE/FALSE)
leading_edge_list <- setNames(fgsea_sub$leadingEdge, fgsea_sub$drug)
all_genes <- unique(unlist(leading_edge_list))

m <- data.frame(Gene = all_genes)
for (drug in names(leading_edge_list)) {
  m[[drug]] <- m$Gene %in% leading_edge_list[[drug]]
}
rownames(m) <- m$Gene
m$Gene <- NULL
m <- as.matrix(m)
m <- m * 1

## integration result
dat <- qread(file= file.path('data/results/drug', "gene_drug_assoc_sts_meta.qs"))
df <- dat[dat$gene_name %in% rownames(m), ]
df <- df[df$drug %in% colnames(m), ] 
df <- df[df$pval < 0.05, ]

m <- m[rownames(m) %in% df$gene_name, ]
min_drug_hits <- ceiling(ncol(m) * 0.5)
m <- m[rowSums(m) >= min_drug_hits, ]

moa_class <- drug_class[drug_class$drug %in% colnames(m), ]
moa_ann <- setNames(moa_class$moa, moa_class$drug)
moa_ann <- moa_ann[colnames(m)]

col = list("MOA" = c("Apoptosis regulation" = "#D2C396FF",
                    # "Cell cycle" = "#78847FFF" ,
                     "Chromatin regulation" = "#ACC2CFFF",
                     "DNA replication" = "#855C75FF",
                    # "EGFR signaling" = "#819574ff", 
                     "ERK MAPK signaling" = "#9C6755FF",
                    # "Genome integrity" = "#877772ff",
                    # "IGF1R signaling" = "#8491BEFF",
                    # "Kinases" = "#8097e8ff",
                    # "Mitosis" = "#B1AF53FF", 
                     "Other" = "#d1d5dfff",
                    # "PI3K MTOR signaling" = "#64894DFF", 
                     "RTK signaling" = "#526A83FF") ) 

#Create the heatmap annotation
ha <- HeatmapAnnotation(
  "MOA" = moa_ann,
   col = col, 
   annotation_name_gp = gpar(fontsize = 7),       
   annotation_legend_param = list(
    title_gp = gpar(fontsize = 7),               
    labels_gp = gpar(fontsize = 6),
    grid_width = unit(3, "mm"),
    grid_height= unit(3, "mm")  
  ))

pdf(file=file.path(dir_out, "heatmap_hallmark_gsea_g2m_v2.pdf"), 
width = 4, height = 7)

Heatmap(
  m,
  name = "Presence",
  col = c("0" = "white", "1" = "#bebbbbff"), 
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation = ha,
  row_title = " ",
  column_title = " ",
  column_order = seq_len(ncol(m)),
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x = x, y = y,
      width = width,
      height = height,
      gp = gpar(col = "#242424", fill = NA, lwd = 0.8)
    )
  }
)

dev.off()

############################################################
## Heatmap figure: G2M_CHECKPOINT
############################################################
## pathway result
fgseaRes <- qread(file.path(dir_in, "hallmark_gsea_pathway_drug.qs"))

# Format leadingEdge column (into character vectors)
fgseaRes$leadingEdge <- strsplit(fgseaRes$leadingEdge, ",")
fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, trimws)

# Clean pathway names
fgseaRes$pathway <- tools::toTitleCase(tolower(substr(fgseaRes$pathway, 10, nchar(fgseaRes$pathway))))

# Keep only selected drugs
fgseaRes <- fgseaRes[fgseaRes$drug %in% selected_drug_clin, ]

# Keep only significant results
fgseaRes <- fgseaRes[fgseaRes$padj < 0.05, ]

id <- "G2m_checkpoint"
core_drugs <- c(
  "Temozolomide", "Venetoclax", "Dabrafenib", "Tamoxifen", "Doxorubicin",
  "Nintedanib", "Panobinostat", "Cabozantinib", "Nilotinib", "Crizotinib"
)
fgsea_sub <- fgseaRes[fgseaRes$drug %in% core_drugs & fgseaRes$pathway == id, ]

# Build matrix of leading-edge genes (TRUE/FALSE)
leading_edge_list <- setNames(fgsea_sub$leadingEdge, fgsea_sub$drug)
all_genes <- unique(unlist(leading_edge_list))

m <- data.frame(Gene = all_genes)
for (drug in names(leading_edge_list)) {
  m[[drug]] <- m$Gene %in% leading_edge_list[[drug]]
}
rownames(m) <- m$Gene
m$Gene <- NULL
m <- as.matrix(m)

## integration result
dat <- qread(file= file.path('data/results/drug', "gene_drug_assoc_sts_meta.qs"))
df <- dat[dat$gene_name %in% rownames(m), ]
df <- df[df$drug %in% colnames(m), ] 
df <- df[df$pval < 0.05, ]
m <- m[rownames(m) %in% df$gene_name, ]

m <- t(m) * 1
min_drug_hits <- ceiling(nrow(m) * 0.5)
m <- m[, colSums(m) >= min_drug_hits]

moa_class <- drug_class[drug_class$drug %in% rownames(m), ]
moa_ann <- setNames(moa_class$moa, moa_class$drug)
moa_ann <- moa_ann[rownames(m)]

col = list("MOA" = c("Apoptosis regulation" = "#D2C396FF",
                    # "Cell cycle" = "#78847FFF" ,
                     "Chromatin regulation" = "#ACC2CFFF",
                     "DNA replication" = "#855C75FF",
                    # "EGFR signaling" = "#819574ff", 
                     "ERK MAPK signaling" = "#9C6755FF",
                    # "Genome integrity" = "#877772ff",
                    # "IGF1R signaling" = "#8491BEFF",
                    # "Kinases" = "#8097e8ff",
                    # "Mitosis" = "#B1AF53FF", 
                     "Other" = "#d1d5dfff",
                    # "PI3K MTOR signaling" = "#64894DFF", 
                     "RTK signaling" = "#526A83FF") ) 

#Create the heatmap annotation
ha <- rowAnnotation(
  "MOA" = moa_ann,
   col = col, 
   annotation_name_gp = gpar(fontsize = 7),       
   annotation_legend_param = list(
    title_gp = gpar(fontsize = 7),               
    labels_gp = gpar(fontsize = 6),
    grid_width = unit(3, "mm"),
    grid_height= unit(3, "mm")  
  ))

pdf(file=file.path(dir_out, "heatmap_hallmark_gsea_g2m.pdf"), 
width = 9, height =3)

Heatmap(
  m,
  name = "Presence",
  col = c("0" = "white", "1" = "#bebbbbff"), 
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  left_annotation = ha,
  row_title = " ",
  column_title = " ",
  column_order = seq_len(ncol(m)),
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(
      x = x, y = y,
      width = width,
      height = height,
      gp = gpar(col = "#242424", fill = NA, lwd = 0.8)
    )
  }
)

dev.off()


