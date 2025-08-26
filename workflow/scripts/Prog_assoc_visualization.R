# -----------------------------------------------------------
# STS Gene–Drug Association Visualization Script
#
# This script visualizes results from a meta-analysis of gene–drug
# associations in soft-tissue sarcoma (STS) cell lines, using aligned
# transcriptomic and pharmacogenomic data from CCLE/CTRP, GDSCv2, and NCI-Sarcoma.
#
# Key outputs:
#   - Barplot and upset plot summarizing significant gene–drug associations
#   - Volcano plots for selected drugs 
#   - Heatmap and boxplot illustrating Linsitinib response and expression
#
# -----------------------------------------------------------
#############################################################
## Load libraries
#############################################################

library(qs)
library(UpSetR)
library(ggrepel)
library(ggplot2)
library(VennDiagram)
library(ComplexHeatmap)
library(paletteer)
library(dplyr)

##################################################################
## Setup directory
##################################################################
dir_in <- 'data/results/drug' 
dir_out <- 'data/results/drug' 

##################################################################
## Load data
##################################################################
dat.meta <- qread(file.path(dir_in, "gene_drug_assoc_sts_meta.qs"))
dat.ccle <- qread(file.path(dir_in, "gene_drug_assoc_sts_ccle_ctrp.qs"))
dat.gdsc <- qread(file.path(dir_in, "gene_drug_assoc_sts_gdsc.qs"))
dat.nci <- qread(file.path(dir_in, "gene_drug_assoc_sts_nci.qs"))

r.cutoff <- 0.30

sig.meta <- dat.meta[dat.meta$padj < 0.05 & abs(dat.meta$r) >= r.cutoff, ]
sig.ccle <- dat.ccle[dat.ccle$padj < 0.05 & abs(dat.ccle$estimate) >= r.cutoff, ]
sig.gdsc <- dat.gdsc[dat.gdsc$padj < 0.05 & abs(dat.gdsc$estimate) >= r.cutoff, ]
sig.nci <- dat.nci[dat.nci$padj < 0.05 & abs(dat.nci$estimate) >= r.cutoff, ]

# seelcted drugs
selected_drugs <- read.csv(file = file.path(dir_in, 'selected_drugs.csv'))

###################################################################
## Bar plot of each studies vs meta
###################################################################
before_union <- TRUE

if (before_union) {
  # union of genes across studies per drug
  before_tbl <- dplyr::bind_rows(
      dplyr::select(sig.ccle, Ensembl_ID, drug),
      dplyr::select(sig.gdsc, Ensembl_ID, drug),
      dplyr::select(sig.nci,  Ensembl_ID, drug)
    ) |>
    dplyr::distinct(Ensembl_ID, drug) |>
    dplyr::count(drug, name = "n_before")
} else {
  # sum of counts across studies per drug
  cnt_ccle <- sig.ccle |> dplyr::count(drug, name = "n_ccle")
  cnt_gdsc <- sig.gdsc |> dplyr::count(drug, name = "n_gdsc")
  cnt_nci  <- sig.nci  |> dplyr::count(drug, name = "n_nci")

  before_tbl <- cnt_ccle |>
    dplyr::full_join(cnt_gdsc, by = "drug") |>
    dplyr::full_join(cnt_nci,  by = "drug") |>
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("n_"), ~ tidyr::replace_na(.x, 0L)),
      n_before = n_ccle + n_gdsc + n_nci
    ) |>
    dplyr::select(drug, n_before)
}

after_tbl <- sig.meta |>
  dplyr::distinct(Ensembl_ID, drug) |>
  dplyr::count(drug, name = "n_after")

# Align drugs and fill zeros where missing
dat <- dplyr::full_join(before_tbl, after_tbl, by = "drug") |>
  dplyr::mutate(
    n_before = tidyr::replace_na(n_before, 0L),
    n_after  = tidyr::replace_na(n_after,  0L)
  )

# Long format for plotting
dat_long <- dat |>
  tidyr::pivot_longer(
    cols = c(n_before, n_after),
    names_to = "Stage",
    values_to = "n_sig"
  ) |>
  dplyr::mutate(
    Stage = dplyr::recode(Stage,
                          n_before = "Before integration",
                          n_after  = "After integration")
  )

# Order drugs by "After" (then total) for nicer plotting
drug_order <- dat |>
  dplyr::mutate(total = n_before + n_after) |>
  dplyr::arrange(dplyr::desc(n_after), dplyr::desc(total)) |>
  dplyr::pull(drug)

dat_long <- dat_long |>
  dplyr::mutate(
    drug  = factor(drug, levels = drug_order),
    Stage = factor(Stage, levels = c("Before integration", "After integration"))
  )

dat_long <- dat_long[dat_long$drug %in% selected_drugs[selected_drugs$sts == 'Yes', "drug"], ]
####################### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ############################

pdf(file= file.path(dir_out, "bar_sig_assoc_meta_studies_clinicalapproved.pdf"),
     width = 8, height = 5)

p <- ggplot2::ggplot(dat_long,
              ggplot2::aes(x = drug, y = n_sig, fill = Stage)) +
  ggplot2::geom_bar(stat = "identity", width = 0.5) +
  ggplot2::scale_fill_manual(
    values = c("After integration" = "#d5ced7ff",   
               "Before integration"  = "#401840ff") )  +
  ggplot2::labs(
    x = " ",
    y = "drug response–associated genes (FDR < 0.05, |r| ≥ 0.30)",
    title = " "
  ) +
  ggplot2::theme(axis.text.x=element_text(size=8, angle = 45, hjust = 1),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position=c(0.5,0.5),
        legend.text = element_text(size = 7.5, face="bold"),
        legend.title = element_blank())

p

dev.off()

###################################################################
## violin plots studies vs meta
###################################################################

ccle_long <- dat.ccle %>%
  dplyr::transmute(gene_name, drug, r = estimate, padj, study = "CCLE/CTRP")
gdsc_long <- dat.gdsc %>%
  dplyr::transmute(gene_name, drug, r = estimate, padj, study = "GDSC")
nci_long  <- dat.nci %>%
  dplyr::transmute(gene_name, drug, r = estimate, padj, study = "NCI-Sarcoma")
meta_long <- dat.meta %>%
  dplyr::transmute(gene_name, drug, r = r, padj, study = "Meta")

eff_long <- dplyr::bind_rows(ccle_long, gdsc_long, nci_long, meta_long)
eff_long <- eff_long [eff_long $padj < 0.05 & abs(eff_long$r) >= 0.30, ]

eff_long <- eff_long %>%
  dplyr::mutate(study = factor(study, levels = c("CCLE/CTRP","GDSC","NCI-Sarcoma","Meta")))

pdf(file= file.path(dir_out, "violin_sig_assoc_meta_estimates.pdf"),
     width = 5, height = 4)

p <- ggplot2::ggplot(eff_long, ggplot2::aes(x = study, y = r, fill = study)) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.85) +
  ggplot2::geom_boxplot(width = 0.05, outlier.shape = NA, color = "black") +
    ggplot2::scale_fill_manual(
    values=c('CCLE/CTRP'="#1b7837", 
             'GDSC' = "#542788", 
            'NCI-Sarcoma' = "#bf812d",
            "Meta" = "#d5ced7ff")
 )  +
  ggplot2::coord_cartesian(ylim = c(-1, 1)) +        # for correlations
  ggplot2::labs(
    x = NULL,
    y = "gene–drug correlation (FDR < 0.05, |r| ≥ 0.30)",
    title = ""
  ) +
    ggplot2::theme(axis.text.x=element_text(size=8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text(size=10),
        axis.text.y=element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position='none',
        legend.text = element_text(size = 7.5, face="bold"),
        legend.title = element_blank())

p

dev.off()

###################################################################
## Bar plot of late integration
###################################################################
sig.meta <- sig.meta[!is.na(sig.meta$r), ]
drug <- unique(sig.meta$drug)
sig <- sapply(1:length(drug), function(k){
  
  sub.meta <- sig.meta[sig.meta$drug == drug[k], ]
  nrow(sub.meta[sub.meta$padj < 0.05 & abs(sub.meta$r) >= r.cutoff, ]) 
  
})

df <- data.frame(number_sig = sig, drug= drug) 
df <- df[df$number_sig != 0 , ]
df <- df[order(df$number_sig, decreasing = TRUE), ]
df$drug <- ifelse( df$number_sig <= 10, "Other", df$drug )

df[df$drug == "Unii-40E3azg1MX", "drug"] <- "BMS-536924" 
df <- df[df$drug %in% selected_drugs[selected_drugs$sts == 'Yes', "drug"], ]

pdf(file= file.path(dir_out, "bar_sig_assoc_meta_30PercCutoffCor_clinicalapproved.pdf"),
     width = 6, height = 6)

p <- ggplot(df, aes(x = reorder(drug, -number_sig), y = number_sig)) +
  geom_bar(width = 0.4, stat = "identity") +
  scale_fill_manual(values = c("#96A5A5FF")) +
  coord_flip()+
  ylab(paste("drug response–associated genes (FDR < 0.05, |r| ≥ 0.30)")) +
  xlab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_text(size=6,  face="bold"),
        axis.title.x=element_text(size=8,face="bold"),
        axis.title.y=element_text(size=8,face="bold"),
        axis.text.y=element_text(size=6, face = "bold"),
        strip.text = element_text(size=8, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position=c(0.85,0.85),
        legend.text = element_text(size = 7.5, face="bold"),
        legend.title = element_blank())

p

dev.off()

################################################################################
## Upset plot 
################################################################################
drug <- unique(sig.meta$drug)
sig <- sapply(1:length(drug), function(k){
  
  sub.meta <- sig.meta[sig.meta$drug == drug[k], ]
  nrow(sub.meta[sub.meta$padj < 0.05 & abs(sub.meta$r) >= r.cutoff, ])
  
})

df <- data.frame(number_sig = sig, drug= drug) 
df <- df[df$number_sig != 0, ]
df <- df[order(df$number_sig, decreasing = TRUE), ]
df$drug <- ifelse( df$number_sig <= 10, "Other", df$drug )
df[df$drug == "Unii-40E3azg1MX", "drug"] <- "BMS-536924" 
df <- df[df$drug %in% selected_drugs[selected_drugs$sts == 'Yes', "drug"], ]

df <- df[1:10, ] # select top 20 drugs

sig_gene_drug <- lapply(1:nrow(df), function(k){
  
  sub_dat <- sig.meta[sig.meta$drug == df$drug[k], ]
  sub_dat <- sub_dat[sub_dat$padj < 0.05 & abs(sub_dat$r) >= r.cutoff, ]
  sub_dat$gene_name
  
})

names(sig_gene_drug) <- df$drug

pdf(file= file.path(dir_out, "upset_sig_assoc_meta_top10_drugs_30PercCutoffCor_clinicalapproved.pdf"),
     width = 12, height = 7)

upset(fromList(sig_gene_drug), nsets = 100 ,  order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE),
      mainbar.y.label = "drug response–associated genes (FDR < 0.05, |r| ≥ 0.30)", text.scale =1.2, 
      matrix.color = "#023743FF", main.bar.color = "#72874EFF", sets.bar.color = "#023743FF")

dev.off()

################################################################################
## Volcano plots STS 
################################################################################

drug <- unique(sig.meta$drug)
sig <- sapply(1:length(drug), function(k){
  
 nrow(sig.meta[sig.meta$drug == drug[k], ])
  
})

df <- data.frame(number_sig = sig, drug= drug) 
df <- df[order(df$number_sig, decreasing = TRUE), ]
df <- df[df$number_sig >= 10, ]

drug <- unique(df$drug)

for(i in 1:length(drug)){
  
  res <- dat.meta[dat.meta$drug == drug[i], ]
  
  res$diffexpressed <- "NO"
  res$diffexpressed[res$r > r.cutoff & res$padj < 0.05] <- "FDR < 0.05, r > 0.30"
  res$diffexpressed[res$r < (-r.cutoff) & res$padj < 0.05] <- "FDR < 0.05, r < -0.30"  
  
  mycolors <- c( "#EF8A62","#67A9CF", "#999999" )
  names(mycolors) <- c("FDR < 0.05, r > 0.30", 
                       "FDR < 0.05, r < -0.30", 
                       "NO")  
 
  res$delabel <- NA
  res <- res[order(res$padj, decreasing = FALSE), ]
  id_pos <- res[res$r > r.cutoff & res$padj < 0.05, "gene_name"][1:25]
  id_neg <- res[res$r < (-r.cutoff) & res$padj < 0.05, "gene_name"][1:25]
  id <- c(id_pos, id_neg)
  id <- id[!is.na(id)]
  
  for(j in 1:length(id)){
    k <- which(res$gene_name == id[j])
    res$delabel[k] <- res[k, ]$gene_name
  }

 pdf(file=paste(paste(file.path(dir_out, "volcano/cutoff_0.30"), 
          paste("volcano", drug[i], sep="_"), sep="/"), "pdf", sep="."), 
     width = 6, height = 6)
  
  p <- ggplot(data=res, aes(x=r, y= -log10(pval), 
                            col=diffexpressed)) + 
    geom_point(size = 1.7) + 
   # geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", color = "darkgray", size = 0.4) +
   # geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "darkgray", size = 0.4) +
    ylab("-log10 P value") +  
    xlab("pooled correlation estimate") +
    scale_colour_manual(values = mycolors) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=12,face="bold"),
          axis.text.y=element_text(size=10, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position="bottom",
          legend.text = element_text(size = 10, face="bold"),
          legend.title = element_blank()) +
    geom_text_repel(aes(label= delabel),
                    size = 2.5,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both", 
                    seed = 2356,
                    fontface= "bold",
                    max.overlaps = 50)
  
  print(p)
  
  dev.off()
  
}

################################################################################
## Heatmap: top 100 genes nad or significant ones
################################################################################
drug <- unique(sig.meta$drug)
sig <- sapply(1:length(drug), function(k){
  
  sub.meta <- sig.meta[sig.meta$drug == drug[k], ]
  nrow(sub.meta[sub.meta$padj < 0.05 & abs(sub.meta$r) >= r.cutoff, ])
  
})

df <- data.frame(number_sig = sig, drug= drug) 
df <- df[df$number_sig != 0, ]
df <- df[order(df$number_sig, decreasing = TRUE), ]
df$drug <- ifelse( df$number_sig <= 10, "Other", df$drug )

drugs <- df[df$number_sig >= 100, 'drug']

sig.meta <- sig.meta[!is.na(sig.meta$r), ]
sig <- sig.meta[sig.meta$padj < 0.05 & abs(sig.meta$r) >= r.cutoff, ]
sig <- sig[order(sig$padj), ]

for(k in 1:length(drugs)){

print(k)
gene_name <- sig[sig$drug == drugs[k], "gene_name"][1:100]

# load drug response data
dat <- qread(file.path(dir_in, "drug_rna_aligned_sts.qs"))   
pset_aac <- dat$pset_aac # it needs to be changed

pset_aac <- pset_aac[rownames(pset_aac) == drugs[k], ]

# load aligned expression data
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

df <- pset_aligned[rownames(pset_aligned) %in% gene_name, ]
expr <- t(scale(t(df)))
pset_aac <- pset_aac[, !is.na(pset_aac), ]
int <- intersect(colnames(pset_aac), colnames(expr))
pset_aac <- as.numeric(pset_aac[, int])
expr <- expr[, int]

group <- colnames(expr)
group[grep("-ccle", colnames(expr))] <- "CCLE"
group[grep("-gdsc", colnames(expr))] <- "GDSC"
group[grep("-nci", colnames(expr))] <- "NCI"

col = list( "PSets" = c( "CCLE" = "#1b7837",
                         "GDSC" = "#542788",
                         "NCI" = "#bf812d"),
            "AAC" = colorRamp2(c(0, median(pset_aac), max(pset_aac)), 
                               c("#354823FF", "white",  "#D4613EFF" )) ) 

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  "PSets" = group, 
  "AAC" = pset_aac,
  col = col,
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(0.5, "cm")
)

# Combine the heatmap and the annotation
pdf(file=file.path(dir_out, "heatmap/cutoff_0.30" , 
         paste('heatmap', paste(drugs[k], '.pdf', sep=""), sep="_")),
          width =6, height = 8)

ht_data <- Heatmap(expr, name = "expr",
                   top_annotation = ha,
                   show_row_names = TRUE,
                   #column_split = split,
                   row_names_gp = gpar(fontsize = 6),
                   #column_names_gp = gpar(fontsize = 8),
                   show_column_names = FALSE,
                   cluster_columns = TRUE,
                   cluster_rows = TRUE,
                   column_title = NULL,
                   show_column_dend = FALSE,
                   show_row_dend = TRUE,
                   colorRamp2(c(min(expr), median(expr), max(expr)), 
                              c("#053061", "white", "#67001f")))

draw(ht_data,
     column_title_gp = gpar(fontsize = 8, fontface = "bold"),
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side="bottom")

dev.off()


}

################################################################################
## Boxplot for AAC distributions 
################################################################################
# consider specific drug 'Linsitinib'
dat <- qread(file.path(dir_in, "drug_rna_aligned_sts.qs"))   
pset_aac <- dat$pset_aac # it needs to be changed
pset_aligned <- dat$pset_aligned
gene_ann <- dat$gene_ann

# remove drugs with more than 70% missing
missing_aac <- sapply(1:nrow(pset_aac), function(k){
  length(which(is.na(pset_aac[k,])))
})

missing_df <- data.frame(perc = round((missing_aac/ncol(pset_aac)) * 100), drug= rownames(pset_aac)) 
missing_df <- missing_df[missing_df$perc < 70, ] # no missing with cut-off 70%

pset_aac <- pset_aac[rownames(pset_aac) %in% missing_df$drug, ] 
durgs <- missing_df$drug

for(k in 1:length(drugs)){
  
df <- pset_aac[rownames(pset_aac) == drugs[k], ]
group <- names(df)
group[grep("-ccle", names(df))] <- "CCLE/CTRP"
group[grep("-gdsc", names(df))] <- "GDSC"
group[grep("-nci", names(df))] <- "NCI"

df <- data.frame(cell = names(df),
                 PSets = group,
                 AAC = as.numeric(df))
df <- df[!is.na(df$AAC), ]


pdf(file= file.path(dir_out, 'boxplot', 
          paste('boxplot', paste(drugs[k], '.pdf', sep=""), sep="_")), 
          width = 3, height = 3)

p <- ggplot(df, aes(x=reorder(PSets, -AAC), y=AAC, fill = PSets)) + 
  geom_boxplot(width = 0.5) +
  scale_fill_manual(values = c("#1b7837", "#542788", "#bf812d")) +
  coord_flip() + 
  #ylim(c(-0.5, 1)) +
  xlab("") +
  ylab(paste(drugs[k], "drug response", sep=" ")) +
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size=9),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.text = element_text(size = 6, face="bold"),
        legend.title = element_blank()) 

print(p)

dev.off()

}
