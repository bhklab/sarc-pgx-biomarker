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
#
# Notes:
#   - Significance thresholds:
#       * GSEA: FDR < 0.05
#       * ORA:  FDR < 0.1
#   - Input data loaded from `*.qs` files in `dir_in`
#   - Figures saved to `dir_out`
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

thr_gsea <- 0.05
thr_ora <- 0.10
top_cutoff <- 30

#################################################################
## load selected drugs
#################################################################

selected_drug <- read.csv(file.path('data/results/drug', 'selected_drugs.csv'))
selected_drug_clin <- selected_drug[selected_drug$sts == 'Yes', 'drug']

##################################################################################################################
################################################### HALLMARK (GSEA) ##############################################
##################################################################################################################
#################################################################
## HALLMARK pathway significant results (keep direction of NES)
#################################################################
fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
gsea_data <- fgseaRes[fgseaRes$drug %in% selected_drug_clin, ]
gsea_data$pathway <- substr(gsea_data$pathway, 10, nchar(gsea_data$pathway))

# compute median NES per pathway
pathway_order <- gsea_data %>%
  filter(padj <= 0.10) %>%
  group_by(pathway) %>%
  summarise(median_NES = median(NES, na.rm = TRUE)) %>%
  arrange(desc(median_NES)) %>%
  pull(pathway)

# order pathways by median NES (signed) across drugs
gsea_data_filtered <- gsea_data %>%
  filter(padj <= 0.10) %>%
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
## HALLMARK pathway significant results (abs of NES)
############################################################
fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
gsea_data <- fgseaRes[fgseaRes$drug %in% selected_drug_clin, ]
gsea_data$pathway <- substr(gsea_data$pathway, 10, nchar(gsea_data$pathway))

# compute median NES per pathway
pathway_order <- gsea_data %>%
  filter(padj <= 0.10) %>%
  group_by(pathway) %>%
  summarise(median_abs_NES = median(abs(NES), na.rm = TRUE)) %>%
  arrange(desc(median_abs_NES)) %>%
  pull(pathway)

# order pathways by median NES (signed) across drugs
gsea_data_filtered <- gsea_data %>%
  filter(padj <= 0.10) %>%
  mutate(
    FDR_bin = cut(padj,
                  breaks = c(-Inf, 0.01, 0.05, 0.10, Inf),
                  labels = c("< 0.01", "0.01–0.05", "0.05–0.10", "> 0.10")),
    pathway = factor(pathway, levels = pathway_order)
  )


pdf(file=file.path(dir_out, "heatmap_dot_hallmark_gsea_abs.pdf"), 
width = 10, height = 10)

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

fgseaRes$pathway <- substr(fgseaRes$pathway, 10, nchar(fgseaRes$pathway))
foraRes$pathway <- substr(foraRes$pathway, 10, nchar(foraRes$pathway))

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
    gsea_sig = padj_gsea <= 0.05,
    ora_sig  = padj_ora <= 0.05,
    both_sig = gsea_sig & gsea_sig
  )

# drug → class mapping
load(file.path('data/results/MOA', 'moa_info.rda'))
drug_classes <- lapply(1:length(dat_drug_class), function(k){

 data.frame(moa = names(dat_drug_class)[k],
            drug = dat_drug_class[[k]])

})

drug_classes <- do.call(rbind, drug_classes)

merged <- merged %>%
  left_join(drug_classes, by = "drug")

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

pdf(file=file.path(dir_out, "heatmap_hallmark_gsea_ora_moa.pdf"), 
width = 10, height = 11)

ggplot(class_summary, aes(x = moa, y = pathway, fill = mean_NES)) +
   geom_tile(color = "#686969ff") +
  # Single-drug overlap → circle
  geom_point(
    data = subset(class_summary, glyph == "single"),
    aes(x = moa, y = pathway),
    shape = 19, size = 1.5, color = "black"
  ) +
  # ≥2 drugs overlap → cross
  geom_point(
    data = subset(class_summary, glyph == "robust"),
    aes(x = moa, y = pathway),
    shape = 8, size = 1.5, color = "black", fill = "black"
  ) +
  scale_fill_gradient2(
  low = "#4685A0FF",
  mid = "white",
  high = "#864568FF",
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
## TO BE ADDED











