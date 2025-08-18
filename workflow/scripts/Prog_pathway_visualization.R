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

##################################################################
## Setup directory
##################################################################

dir_in <- 'data/results/pathway' 
dir_out <- 'data/results/pathway/Fig'
dir_moa <- 'data/results/moa'

thr_gsea <- 0.05
thr_ora <- 0.1
top_cutoff <- 30

##################################################################################################################
################################################### HALLMARK (GSEA) ##############################################
##################################################################################################################
############################################################
## Pathway significant results
############################################################
fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
drugs <- unique(fgseaRes$drug)
sig <- sapply(1:length(drugs), function(k){
    
    sub_dat <- fgseaRes[fgseaRes$drug == drugs[k], ]
    nrow(sub_dat[sub_dat$padj < thr_gsea, ])
    
  })
  
df <- data.frame(number_sig = sig, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]
df[df$drug == "Unii-40E3azg1MX", "drug"] <- "BMS-536924" 
  
pdf(file=file.path(dir_out, "bar_hallmark_gsea.pdf"), width = 7, height = 9)
  
  p <- ggplot(df, aes(x = reorder(drug, -number_sig), y = number_sig)) +
    geom_bar(width = 0.4, stat = "identity", fill = "#a6bddb") +
    coord_flip()+
    ylab(paste(paste("HALLMARK", "drug assoication", sep="-"), 
               "\n (FDR < 0.05)", sep="")) +
    xlab("") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      legend.position.inside = c(0.85, 0.85),  
      legend.text = element_text(size = 6, face = "bold"),
      legend.title = element_blank()
    )
  
  print(p)
  
dev.off()

#############################################################################
## Pathway dot plot
#############################################################################

fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
fgseaRes$pathway <- substr(fgseaRes$pathway, 10, nchar(fgseaRes$pathway))
sig <- fgseaRes[fgseaRes$padj < thr_gsea, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])
  
})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]

drugs <- df$drug

for(i in 1:length(drugs)){
  
  df <- sig[sig$drug == drugs[i], ]
  
  pdf(file.path(dir_out, 'GSEA/HALLMARK', paste(drugs[i], "hallmark_FDR_0.05.pdf", sep="_")),
      width = 10, height = 7)
  # Plot significant pathways as points with a gradient color based on p-value and size based on a 'size' variable.
  p <- ggplot(df, aes(y=reorder(pathway, size),x= NES)) + geom_point(aes(color=padj, size=size)) + 
    scale_color_gradientn(colours = c("#084594","#b10026")) + 
    labs(x='normalized enrichment score', y=NULL ) + 
    theme(
      axis.text.x=element_text(size=10),
      axis.title=element_text(size=10, face = "bold"),
      axis.text.y=element_text(size=8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # legend.position="bottom",
      legend.text = element_text(size = 8, face="bold"))
  
  print(p)
  
  dev.off()
  
}

################################################################################
## Upset plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
sig <- fgseaRes[fgseaRes$padj < thr_gsea, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])

})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]

df <- df[df$drug %in% c("Doxorubicin", "Gemcitabine", "Eribulin",
                        "Trabectedin", "Docetaxel", "Dacarbazine",
                        "Pazopanib", "Tazemetostat", "SN-38",
                        "Axitinib", "Linsitinib", "Methotrexate", 
                        "Dabrafenib", "Vorinostat", "Selumetinib"), ]

sig_gene_drug <- lapply(1:nrow(df), function(k){
  
  sub_dat <- fgseaRes[fgseaRes$drug == df$drug[k] & fgseaRes$padj < thr_gsea, ]
  sub_dat$pathway 
  
})

names(sig_gene_drug) <- df$drug

pdf(file=file.path(dir_out, "upset_hallmark_gsea.pdf"),
     width = 12, height = 8)

upset(fromList(sig_gene_drug), nsets = 50 , order.by = "freq", 
      mainbar.y.label = "HALLMARK-drug association (FDR < 0.05)", text.scale =1.2, 
      matrix.color = "#4393c3", main.bar.color = "#053061", sets.bar.color = "#4393c3")

dev.off()

################################################################################
## Heatmap plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))
sig_dat <- fgseaRes[fgseaRes$padj < thr_gsea, ]

drugs <- names(table(sig_dat$drug)[order(table(sig_dat$drug), decreasing  = TRUE)])
sig_dat <- sig_dat[sig_dat$drug %in% drugs, ]
pathway_id <- unique(sig_dat$pathway)
a <- table(sig_dat$pathway, sig_dat$drug)
b <- apply(a, 1, sum)
b <- b[b >= 2]

pathway_heatmap <- lapply(1:length(drugs), function(k){
  
  sub_sig_dat <- sig_dat[sig_dat$drug == drugs[k], ]
  sub_pathway_name <- sub_sig_dat$pathway
  
  cor_val <- sapply(1:length(pathway_id), function(i){
    
    if(pathway_id[i] %in% sub_pathway_name){
      
      val <- sub_sig_dat[sub_sig_dat$pathway == pathway_id[i], "NES"]
      
    }else{ val <- NA}
    
    val
  })
  
  cor_val
  
})

pathway_heatmap <- do.call(cbind, pathway_heatmap)
rownames(pathway_heatmap) <- substr(pathway_id, 10, nchar(pathway_id))
colnames(pathway_heatmap) <- drugs

pdf(file=file.path(dir_out, "heatmap_hallmark_gsea.pdf"),
     width = 16, height = 9)

p1 <- Heatmap(pathway_heatmap, 
              col = colorRamp2(c(-4, 0, 4), c("#006837", "#e0e0e0", "#8e0152")),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_names = TRUE, 
              show_column_names = TRUE,
              na_col = "#e0e0e0",
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(legend_direction = "horizontal",
                                          title = "enrichment score",
                                          title_position = "topcenter"),
              rect_gp = gpar(col = "grey", lwd = 0.5),
              width = ncol(pathway_heatmap)*unit(5, "mm")
              )

p <-  draw(p1, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")  
p

dev.off()

##################################################################################################################
################################################### HALLMARK (ORA) ##############################################
##################################################################################################################
############################################################
## Pathway significant results
############################################################
fgseaRes <- qread(file.path(dir_in, 'hallmark_ora_pathway_drug.qs'))
drugs <- unique(fgseaRes$drug)
sig <- sapply(1:length(drugs), function(k){
    
    sub_dat <- fgseaRes[fgseaRes$drug == drugs[k], ]
    nrow(sub_dat[sub_dat$p.adjust < thr_ora, ])
    
  })
  
df <- data.frame(number_sig = sig, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]
df[df$drug == "Unii-40E3azg1MX", "drug"] <- "BMS-536924" 
  
pdf(file=file.path(dir_out, "bar_hallmark_ora.pdf"), width = 7, height = 7)
  
  p <- ggplot(df, aes(x = reorder(drug, -number_sig), y = number_sig)) +
    geom_bar(width = 0.4, stat = "identity", fill = "#a6bddb") +
    coord_flip()+
    ylab(paste(paste("HALLMARK", "drug assoication", sep="-"), 
               "\n (FDR < 0.1)", sep="")) +
    xlab("") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      legend.position.inside = c(0.85, 0.85),  # ✅ new argument
      legend.text = element_text(size = 6, face = "bold"),
      legend.title = element_blank()
    )
  
  print(p)
  
dev.off()

#############################################################################
## Pathway dot plot
#############################################################################

fgseaRes <- qread(file.path(dir_in, 'hallmark_ora_pathway_drug.qs'))
fgseaRes <- fgseaRes[order(fgseaRes$p.adjust), ]
fgseaRes$ID <- substr(fgseaRes$ID, 10, nchar(fgseaRes$ID))
sig <- fgseaRes[fgseaRes$p.adjust < thr_ora, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])
  
})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]

drugs <- df$drug

for(i in 1:length(drugs)){
  
  df <- sig[sig$drug == drugs[i], ]
  
  pdf(file.path(dir_out, 'ORA/HALLMARK', paste(drugs[i], "hallmark_FDR_0.1.pdf", sep="_")),
      width = 8, height = 5)
  # Plot significant pathways as points with a gradient color based on p-value and size based on a 'size' variable.
  p <- ggplot(df, aes(y=reorder(ID, Count),x= FoldEnrichment)) + geom_point(aes(color=p.adjust, size=Count)) + 
    scale_color_gradientn(colours = c("#084594","#b10026")) + 
    labs(x='fold Enrichment', y=NULL ) + 
    theme(
      axis.text.x=element_text(size=10),
      axis.title=element_text(size=10, face = "bold"),
      axis.text.y=element_text(size=8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # legend.position="bottom",
      legend.text = element_text(size = 8, face="bold"))
  
  print(p)
  
  dev.off()
  
}

################################################################################
## Upset plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'hallmark_ora_pathway_drug.qs'))
sig <- fgseaRes[fgseaRes$p.adjust < thr_ora, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])

})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$drug %in% c("Doxorubicin", "Gemcitabine", "Eribulin",
                        "Trabectedin", "Docetaxel", "Dacarbazine",
                        "Pazopanib", "Tazemetostat", "SN-38",
                        "Axitinib", "Linsitinib", "Methotrexate", 
                        "Dabrafenib", "Vorinostat", "Selumetinib"), ]

sig_gene_drug <- lapply(1:nrow(df), function(k){
  
  sub_dat <- fgseaRes[fgseaRes$drug == df$drug[k] & fgseaRes$p.adjust < thr_ora, ]
  sub_dat$ID
  
})

names(sig_gene_drug) <- df$drug

pdf(file=file.path(dir_out, "upset_hallmark_ora.pdf"),
     width = 12, height = 8)

upset(fromList(sig_gene_drug), nsets = 50 , order.by = "freq", 
      mainbar.y.label = "HALLMARK-drug association (FDR < 0.1)", text.scale =1.2, 
      matrix.color = "#4393c3", main.bar.color = "#053061", sets.bar.color = "#4393c3")

dev.off()

################################################################################
## Heatmap plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'hallmark_ora_pathway_drug.qs'))
sig_dat <- fgseaRes[fgseaRes$p.adjust < thr_ora, ]

drugs <- names(table(sig_dat$drug)[order(table(sig_dat$drug), decreasing  = TRUE)])
sig_dat <- sig_dat[sig_dat$drug %in% drugs, ]
pathway_id <- unique(sig_dat$ID)
tab <- table(sig_dat$ID, sig_dat$drug)           
bin_mat <- (tab > 0) * 1L                             
keep_rows <- rowSums(bin_mat) >= 2
bin_mat <- bin_mat[keep_rows, , drop = FALSE]

# order columns (drugs) by total significant pathways; rows likewise
col_order <- order(colSums(bin_mat), decreasing = TRUE)
row_order <- order(rowSums(bin_mat), decreasing = TRUE)
bin_mat <- bin_mat[row_order, col_order, drop = FALSE]
rownames(bin_mat) <- substr(rownames(bin_mat) , 10, nchar(rownames(bin_mat) ))

pdf(file=file.path(dir_out, "heatmap_hallmark_ora.pdf"),
     width = 16, height = 9)

p <- Heatmap(
  bin_mat,
  name = "Significance",
  col = c(`0` = "white", `1` = "#2C7FB8"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  na_col = "white",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "grey", lwd = 0.5),
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("NS", "FDR < 0.1"),
    legend_direction = "horizontal",
    title_position = "topcenter",
    title = NULL
  ),
  width = ncol(bin_mat) * unit(5, "mm")
)

draw(p, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")

dev.off()

##################################################################################################################
################################################## HALLMARK GSEA MOA #############################################
##################################################################################################################
## load class of drugs information
load(file.path(dir_moa, 'moa_info.rda'))
moa_names <- names(dat_drug_class)

## load GSEA results
fgseaRes <- qread(file.path(dir_in, 'hallmark_gsea_pathway_drug.qs'))

for (j in 1:length(moa_names)){

sig_dat <- fgseaRes[fgseaRes$padj < thr_gsea, ]
drugs <- dat_drug_class[[j]]
sig_dat <- sig_dat[sig_dat$drug %in% drugs, ]
pathway_id <- unique(sig_dat$pathway)
a <- table(sig_dat$pathway, sig_dat$drug)
b <- apply(a, 1, sum)
b <- b[b >= 2]

pathway_heatmap <- lapply(1:length(drugs), function(k){
  
  sub_sig_dat <- sig_dat[sig_dat$drug == drugs[k], ]
  sub_pathway_name <- sub_sig_dat$pathway
  
  cor_val <- sapply(1:length(pathway_id), function(i){
    
    if(pathway_id[i] %in% sub_pathway_name){
      
      val <- sub_sig_dat[sub_sig_dat$pathway == pathway_id[i], "NES"]
      
    }else{ val <- NA}
    
    val
  })
  
  cor_val
  
})

pathway_heatmap <- do.call(cbind, pathway_heatmap)
rownames(pathway_heatmap) <- substr(pathway_id, 10, nchar(pathway_id))
colnames(pathway_heatmap) <- drugs

pdf(file=file.path(dir_out, 'heatmap/HALLMARK', paste(moa_names[j], "heatmap_hallmark_gsea.pdf", sep='_')),
     width = 8, height = 8)

p1 <- Heatmap(pathway_heatmap, 
              col = colorRamp2(c(-4, 0, 4), c("#006837", "#e0e0e0", "#8e0152")),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_names = TRUE, 
              show_column_names = TRUE,
              na_col = "#e0e0e0",
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(legend_direction = "horizontal",
                                          title = "enrichment score",
                                          title_position = "topcenter"),
              rect_gp = gpar(col = "grey", lwd = 0.5),
              width = ncol(pathway_heatmap)*unit(5, "mm")
              )

p <-  draw(p1, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")  
p

dev.off()

}





















gsea <- qread(file.path(dir_in, "hallmark_gsea_pathway_drug.qs"))
ora  <- qread(file.path(dir_in, "hallmark_ora_pathway_drug.qs"))

# Harmonize pathway names like you do elsewhere
norm_name <- function(x) substr(x, 10, nchar(x))

sig_gsea <- gsea[gsea$padj < thr_gsea, ]
sig_ora <- ora[ora$padj < thr_ora, ]
drugs <- intersect(sig_gsea$drug, sig_ora$drug)

for(i in 1:length(drugs)){
 
 print(i)
gsea_d <- sig_gsea[sig_gsea$drug == drugs[i], ]
ora_d <- sig_ora[sig_ora$drug == drugs[i], ]

# Sets
set_gsea <- unique(gsea_d$pathway)
set_ora  <- unique(ora_d$pathway)

pdf(file=file.path(dir_out, 'venn/HALLMARK', paste(drugs[i], "ora_gsea.pdf", sep="_")),
     width = 5.5, height = 5.5)

venn.plot <- draw.pairwise.venn(
  area1 = length(set_gsea),
  area2 = length(set_ora),
  cross.area = length(intersect(set_gsea, set_ora)),
  category = c("GSEA", "ORA"),
  fill = c("#98adbdff", "#e5b3caff"),
  alpha = c(0.5, 0.5),
  cex = 1.2, cat.cex = 1.2
)
grid.draw(venn.plot)

dev.off()

}

##################################################################################################################
################################################### GO:BP (GSEA) ##############################################
##################################################################################################################
############################################################
## Pathway significant results
############################################################
fgseaRes <- qread(file.path(dir_in, 'go_gsea_pathway_drug.qs'))
drugs <- unique(fgseaRes$drug)
sig <- sapply(1:length(drugs), function(k){
    
    sub_dat <- fgseaRes[fgseaRes$drug == drugs[k], ]
    nrow(sub_dat[sub_dat$padj < thr_gsea, ])
    
  })
  
df <- data.frame(number_sig = sig, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]
df[df$drug == "Unii-40E3azg1MX", "drug"] <- "BMS-536924" 
  
pdf(file=file.path(dir_out, "bar_go_gsea.pdf"), width = 7, height = 9)
  
  p <- ggplot(df, aes(x = reorder(drug, -number_sig), y = number_sig)) +
    geom_bar(width = 0.4, stat = "identity", fill = "#a6bddb") +
    coord_flip()+
    ylab(paste(paste("GO Biological Process", "drug assoication", sep="-"), 
               "\n (FDR < 0.01)", sep="")) +
    xlab("") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      legend.position.inside = c(0.85, 0.85),  
      legend.text = element_text(size = 6, face = "bold"),
      legend.title = element_blank()
    )
  
  print(p)
  
dev.off()

#############################################################################
## Pathway dot plot
#############################################################################

fgseaRes <- qread(file.path(dir_in, 'go_gsea_pathway_drug.qs'))
fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
fgseaRes$pathway <- substr(fgseaRes$pathway, 6, nchar(fgseaRes$pathway))
sig <- fgseaRes[fgseaRes$padj < thr_gsea, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])
  
})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]

drugs <- df$drug

for(i in 1:length(drugs)){
  
  df <- sig[sig$drug == drugs[i], ]

  if(nrow(df) >= top_cutoff){
    df <- df[1:top_cutoff, ]
  }

  pdf(file.path(dir_out, 'GSEA/GO', paste(drugs[i], "go_FDR_0.01.pdf", sep="_")),
      width = 10, height = 7)
  # Plot significant pathways as points with a gradient color based on p-value and size based on a 'size' variable.
  p <- ggplot(df, aes(y=reorder(pathway, size),x= NES)) + geom_point(aes(color=padj, size=size)) + 
    scale_color_gradientn(colours = c("#084594","#b10026")) + 
    labs(x='normalized enrichment score', y=NULL ) + 
    theme(
      axis.text.x=element_text(size=10),
      axis.title=element_text(size=10, face = "bold"),
      axis.text.y=element_text(size=8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # legend.position="bottom",
      legend.text = element_text(size = 8, face="bold"))
  
  print(p)
  
  dev.off()
  
}

################################################################################
## Upset plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'go_gsea_pathway_drug.qs'))
sig <- fgseaRes[fgseaRes$padj < thr_gsea, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])

})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]

df <- df[df$drug %in% c("Doxorubicin", "Gemcitabine", "Eribulin",
                        "Trabectedin", "Docetaxel", "Dacarbazine",
                        "Pazopanib", "Tazemetostat", "SN-38",
                        "Axitinib", "Linsitinib", "Methotrexate", 
                        "Dabrafenib", "Vorinostat", "Selumetinib", "Bms-754807"), ]

sig_gene_drug <- lapply(1:nrow(df), function(k){
  
  sub_dat <- fgseaRes[fgseaRes$drug == df$drug[k] & fgseaRes$padj < thr_gsea, ]
  sub_dat$pathway 
  
})

names(sig_gene_drug) <- df$drug

pdf(file=file.path(dir_out, "upset_go_gsea.pdf"),
     width = 12, height = 8)

upset(fromList(sig_gene_drug), nsets = 50 , order.by = "freq", 
      mainbar.y.label = "GO:BP-drug association (FDR < 0.01)", text.scale =1.2, 
      matrix.color = "#4393c3", main.bar.color = "#053061", sets.bar.color = "#4393c3")

dev.off()

################################################################################
## Heatmap plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'go_gsea_pathway_drug.qs'))
sig_dat <- fgseaRes[fgseaRes$padj < thr_gsea, ]

drugs <- names(table(sig_dat$drug)[order(table(sig_dat$drug), decreasing  = TRUE)])
sig_dat <- sig_dat[sig_dat$drug %in% drugs, ]
pathway_id <- unique(sig_dat$pathway)
a <- table(sig_dat$pathway, sig_dat$drug)
b <- apply(a, 1, sum)
b <- b[b >= 2]

pathway_heatmap <- lapply(1:length(drugs), function(k){
  
  sub_sig_dat <- sig_dat[sig_dat$drug == drugs[k], ]
  sub_pathway_name <- sub_sig_dat$pathway
  
  cor_val <- sapply(1:length(pathway_id), function(i){
    
    if(pathway_id[i] %in% sub_pathway_name){
      
      val <- sub_sig_dat[sub_sig_dat$pathway == pathway_id[i], "NES"]
      
    }else{ val <- NA}
    
    val
  })
  
  cor_val
  
})

pathway_heatmap <- do.call(cbind, pathway_heatmap)
rownames(pathway_heatmap) <- substr(pathway_id, 6, nchar(pathway_id))
colnames(pathway_heatmap) <- drugs

# fraction of missingness allowed
max_missing_frac <- 0.20    

# number of drugs (columns)
n_drugs <- ncol(pathway_heatmap)

# maximum number of NA values allowed per pathway (20% of drugs)
max_missing <- floor(n_drugs * max_missing_frac)

# keep rows (pathways) with <= max_missing NA values
keep <- rowSums(!is.na(pathway_heatmap)) >= max_missing 
pathway_heatmap_filt <- pathway_heatmap[keep, , drop = FALSE]

# drop columns (drugs) with all NA values
pathway_heatmap_filt <- pathway_heatmap_filt[, colSums(is.na(pathway_heatmap_filt)) != nrow(pathway_heatmap_filt), drop = FALSE]

pdf(file=file.path(dir_out, "heatmap_go_gsea.pdf"),
     width = 16, height = 12)

p1 <- Heatmap(pathway_heatmap_filt, 
              col = colorRamp2(c(-3, 0, 3), c("#006837", "#e0e0e0", "#8e0152")),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_names = TRUE, 
              show_column_names = TRUE,
              na_col = "#e0e0e0",
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(legend_direction = "horizontal",
                                          title = "enrichment score",
                                          title_position = "topcenter"),
              rect_gp = gpar(col = "grey", lwd = 0.5),
              width = ncol(pathway_heatmap)*unit(5, "mm")
              )

p <-  draw(p1, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")  
p

dev.off()

##################################################################################################################
################################################### GO:BP (ORA) ##############################################
##################################################################################################################
############################################################
## Pathway significant results
############################################################
fgseaRes <- qread(file.path(dir_in, 'go_ora_pathway_drug.qs'))
drugs <- unique(fgseaRes$drug)
sig <- sapply(1:length(drugs), function(k){
    
    sub_dat <- fgseaRes[fgseaRes$drug == drugs[k], ]
    nrow(sub_dat[sub_dat$p.adjust < thr_ora, ])
    
  })
  
df <- data.frame(number_sig = sig, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]
df[df$drug == "Unii-40E3azg1MX", "drug"] <- "BMS-536924" 
  
pdf(file=file.path(dir_out, "bar_go_ora.pdf"), width = 7, height = 7)
  
  p <- ggplot(df, aes(x = reorder(drug, -number_sig), y = number_sig)) +
    geom_bar(width = 0.4, stat = "identity", fill = "#a6bddb") +
    coord_flip()+
    ylab(paste(paste("GO Biological Process", "drug assoication", sep="-"), 
               "\n (FDR < 0.1)", sep="")) +
    xlab("") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      legend.position.inside = c(0.85, 0.85),  # ✅ new argument
      legend.text = element_text(size = 6, face = "bold"),
      legend.title = element_blank()
    )
  
  print(p)
  
dev.off()

#############################################################################
## Pathway dot plot
#############################################################################

fgseaRes <- qread(file.path(dir_in, 'go_ora_pathway_drug.qs'))
fgseaRes <- fgseaRes[order(fgseaRes$p.adjust), ]
fgseaRes$ID <- substr(fgseaRes$ID, 6, nchar(fgseaRes$ID))
sig <- fgseaRes[fgseaRes$p.adjust < thr_ora, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])
  
})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig != 0, ]

drugs <- df$drug

for(i in 1:length(drugs)){
  
  df <- sig[sig$drug == drugs[i], ]
  
    if(nrow(df) >= top_cutoff){
    df <- df[1:top_cutoff, ]
  }

  pdf(file.path(dir_out, 'ORA/GO', paste(drugs[i], "go_FDR_0.1.pdf", sep="_")),
      width = 8, height = 5)
  # Plot significant pathways as points with a gradient color based on p-value and size based on a 'size' variable.
  p <- ggplot(df, aes(y=reorder(ID, Count),x= FoldEnrichment)) + geom_point(aes(color=p.adjust, size=Count)) + 
    scale_color_gradientn(colours = c("#084594","#b10026")) + 
    labs(x='fold Enrichment', y=NULL ) + 
    theme(
      axis.text.x=element_text(size=10),
      axis.title=element_text(size=10, face = "bold"),
      axis.text.y=element_text(size=8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.line = element_line(colour = "black"),
      # legend.position="bottom",
      legend.text = element_text(size = 8, face="bold"))
  
  print(p)
  
  dev.off()
  
}

################################################################################
## Upset plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'go_ora_pathway_drug.qs'))
sig <- fgseaRes[fgseaRes$p.adjust < thr_ora, ]
drugs <- unique(sig$drug)

df <- sapply(1:length(drugs), function(k){
  
  nrow(sig[sig$drug == drugs[k], ])

})

df <- data.frame(number_sig = df, drug= drugs) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$drug %in% c("Doxorubicin", "Gemcitabine", "Eribulin",
                        "Trabectedin", "Docetaxel", "Dacarbazine",
                        "Pazopanib", "Tazemetostat", "SN-38",
                        "Axitinib", "Linsitinib", "Methotrexate", 
                        "Dabrafenib", "Vorinostat", "Selumetinib", "Bms-754807"), ]

sig_gene_drug <- lapply(1:nrow(df), function(k){
  
  sub_dat <- fgseaRes[fgseaRes$drug == df$drug[k] & fgseaRes$p.adjust < thr_ora, ]
  sub_dat$ID
  
})

names(sig_gene_drug) <- df$drug

pdf(file=file.path(dir_out, "upset_go_ora.pdf"),
     width = 12, height = 8)

upset(fromList(sig_gene_drug), nsets = 50 , order.by = "freq", 
      mainbar.y.label = "GO:BP-drug association (FDR < 0.1)", text.scale =1.2, 
      matrix.color = "#4393c3", main.bar.color = "#053061", sets.bar.color = "#4393c3")

dev.off()

################################################################################
## Heatmap plot
################################################################################

fgseaRes <- qread(file.path(dir_in, 'go_ora_pathway_drug.qs'))
sig_dat <- fgseaRes[fgseaRes$p.adjust < thr_ora, ]

drugs <- names(table(sig_dat$drug)[order(table(sig_dat$drug), decreasing  = TRUE)])
sig_dat <- sig_dat[sig_dat$drug %in% drugs, ]
pathway_id <- unique(sig_dat$ID)
tab <- table(sig_dat$ID, sig_dat$drug)           
pathway_heatmap <- (tab > 0) * 1L                             
keep_rows <- rowSums(pathway_heatmap) >= 2
pathway_heatmap <- pathway_heatmap[keep_rows, , drop = FALSE]

# order columns (drugs) by total significant pathways; rows likewise
col_order <- order(colSums(pathway_heatmap), decreasing = TRUE)
row_order <- order(rowSums(pathway_heatmap), decreasing = TRUE)
pathway_heatmap <- pathway_heatmap[row_order, col_order, drop = FALSE]
rownames(pathway_heatmap) <- substr(rownames(pathway_heatmap) , 6, nchar(rownames(pathway_heatmap) ))

# drop columns (drugs) with all zero values
pathway_heatmap <- pathway_heatmap[, colSums(pathway_heatmap) != 0, drop = FALSE]


pdf(file=file.path(dir_out, "heatmap_go_ora.pdf"),
     width = 16, height = 16)

p <- Heatmap(
  pathway_heatmap,
  name = "Significance",
  col = c(`0` = "white", `1` = "#2C7FB8"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  na_col = "white",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "grey", lwd = 0.5),
  heatmap_legend_param = list(
    at = c(0, 1),
    labels = c("NS", "FDR < 0.1"),
    legend_direction = "horizontal",
    title_position = "topcenter",
    title = NULL
  ),
  width = ncol(pathway_heatmap_filt) * unit(5, "mm")
)

draw(p, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")

dev.off()

##################################################################################################################
################################################## GO GSEA MOA #############################################
##################################################################################################################
## load class of drugs information
load(file.path(dir_moa, 'moa_info.rda'))
moa_names <- names(dat_drug_class)

## load GSEA results
fgseaRes <- qread(file.path(dir_in, 'go_gsea_pathway_drug.qs'))

for (j in 1:length(moa_names)){

sig_dat <- fgseaRes[fgseaRes$padj < thr_gsea, ]
drugs <- dat_drug_class[[j]]
sig_dat <- sig_dat[sig_dat$drug %in% drugs, ]
pathway_id <- unique(sig_dat$pathway)
a <- table(sig_dat$pathway, sig_dat$drug)
b <- apply(a, 1, sum)
b <- b[b >= 2]

pathway_heatmap <- lapply(1:length(drugs), function(k){
  
  sub_sig_dat <- sig_dat[sig_dat$drug == drugs[k], ]
  sub_pathway_name <- sub_sig_dat$pathway
  
  cor_val <- sapply(1:length(pathway_id), function(i){
    
    if(pathway_id[i] %in% sub_pathway_name){
      
      val <- sub_sig_dat[sub_sig_dat$pathway == pathway_id[i], "NES"]
      
    }else{ val <- NA}
    
    val
  })
  
  cor_val
  
})

pathway_heatmap <- do.call(cbind, pathway_heatmap)
rownames(pathway_heatmap) <- substr(pathway_id, 6, nchar(pathway_id))
colnames(pathway_heatmap) <- drugs

# fraction of missingness allowed
#max_missing_frac <- 0.20    

# number of drugs (columns)
#n_drugs <- ncol(pathway_heatmap)

# maximum number of NA values allowed per pathway (20% of drugs)
#max_missing <- floor(n_drugs * max_missing_frac)

# keep rows (pathways) with <= max_missing NA values
#keep <- rowSums(!is.na(pathway_heatmap)) >= max_missing 
#pathway_heatmap_filt <- pathway_heatmap[keep, , drop = FALSE]

# drop columns (drugs) with all NA values
pathway_heatmap <- pathway_heatmap[, colSums(is.na(pathway_heatmap)) != nrow(pathway_heatmap), drop = FALSE]

if(nrow(pathway_heatmap) >= 150){

pathway_heatmap <- pathway_heatmap[1:150, ]

}

pdf(file=file.path(dir_out, 'heatmap/GO', paste(moa_names[j], "heatmap_go_gsea.pdf", sep='_')),
     width = 12, height = 14)

p1 <- Heatmap(pathway_heatmap, 
              col = colorRamp2(c(-4, 0, 4), c("#006837", "#e0e0e0", "#8e0152")),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              show_row_names = TRUE, 
              show_column_names = TRUE,
              na_col = "#e0e0e0",
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 9),
              heatmap_legend_param = list(legend_direction = "horizontal",
                                          title = "enrichment score",
                                          title_position = "topcenter"),
              rect_gp = gpar(col = "grey", lwd = 0.5),
              width = ncol(pathway_heatmap)*unit(5, "mm")
              )

p <-  draw(p1, show_heatmap_legend = TRUE, heatmap_legend_side = "bottom")  
p

dev.off()

}







