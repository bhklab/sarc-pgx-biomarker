# -------------------------------------------------------------------------
# TCGA Validation and Survival Association Analysis Script
#
# This script performs comprehensive validation of gene–drug associations
# in soft-tissue sarcoma (STS) using TCGA transcriptomic and clinical data.
# The analysis integrates drug sensitivity (from CCLE/CTRP, GDSCv2, NCI-Sarcoma)
# with survival correlations and pathway-level enrichment in STS patients.
#
# Key outputs:
#   - Scatter plots: Cell line vs. TCGA survival correlations (gene-level)
#   - Kaplan–Meier (KM) plots: Pan-STS and LMS subtype for key genes
#   - Boxplots: Gene expression across histological subtypes
#   - Volcano plot: TCGA gene-level survival associations
#   - Pathway volcano plot: Hallmark GSEA enrichment with drug support
#
# Used in: Figures 6D and Supplementary Figures (gene survival and pathway)
# --------------------------------------------------------------------------
########################################################################
## Load library
########################################################################

library(qs)
library(survival)
library(survminer)
library(data.table)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(scales)
library(SummarizedExperiment)

########################################################################
## Set up directroy
########################################################################

dir_in_cl <- 'data/results/drug'
dir_in_cl_pathway <- 'data/results/pathway'
dir_in_tcga <- 'data/results/validation/TCGA'
dir_in_tcga_dat <- 'data/procdata/validation'
dir_out <- 'data/results/validation/TCGA'
dir_moa <- 'data/results/MOA'

thr_ora <- 0.15
thr_gsea <- 0.05

#########################################################################
## load TCGA data
#########################################################################

tcga_dat <- qread(file.path(dir_in_tcga_dat, 'curated_TCGA_se.qs'))
clin <- data.frame(colData(tcga_dat))
expr <- assay(tcga_dat)

genes <- c("DLX2", "STMN1", "SMO", 'BCL2', 'GLI1', "MAPK","RIPK4", "MDM2", 
           "RPTOR", "PRKAA2", "AGL" , "RB1", "ATRX", "MYOCD",
           "MCM2", "PTEN", 'TP53', 'MED12', 'HMGA2', 'RAD51B',
           'TOP2A', 'ASPM', 'BUB1B', 'CEP55', 'PRC1',
           'TXP2', 'ANLN', 'MELK', 'CDC20',
           'KPNA2', 'CENPF', 'NUF2', 'BRCA2',
           'ALK', 'FGFR3', 'FGFR4', 'FLT3',
           'PAX3', 'PAX7', 'RET', 'PGR', 'RIC3', 'DPP6', 'ACKR1', 'TTK', 'KIF4A')
########################################################################
## Load TCGA and Cell Line results (gene level)
########################################################################
selected_drug <- read.csv(file.path(dir_in_cl, 'selected_drugs.csv'))
load(file.path(dir_in_tcga, 'meta_cox_tcga_histo.rda'))
res.meta <- qread(file.path(dir_in_cl, 'gene_drug_assoc_sts_meta.qs'))

drug <- selected_drug[selected_drug$sts == "Yes", "drug"]
sig.meta <- res.meta[res.meta$padj < 0.05 & abs(res.meta$r) >= 0.3 & res.meta$drug %in% drug, ]
sig.cox <- meta.cox[meta.cox$padj < 0.1, ]

int <- lapply(1:length(drug), function(k){

df.meta <- sig.meta[sig.meta$drug == drug[k], ]
data.frame(drug = drug[k],
           n = length(intersect(sig.cox$Gene, df.meta$gene_name)) )

})

int <- do.call(rbind, int) 
int <- int[order(int$n, decreasing=TRUE), ]
selected_drug <- int[int$n >= 5, 'drug']


#######################################################
## scatter plot
#######################################################
cor_res <- numeric()

for(i in 1:length(selected_drug)){

# Filter cell line data for Axitinib + FDR + r cutoff
drug_genes <- res.meta %>%
  filter(drug == selected_drug[i], padj < 0.05, abs(r) >= 0.3) 

cox_genes <- meta.cox %>%
  filter(padj < 0.1) 

merged <- merge(
  drug_genes,
  cox_genes,
  by.x = "gene_name",
  by.y = "Gene"
)

merged$tcga_sig <- cut(
  merged$padj.y,
  breaks = c(-Inf, 0.05, 0.1, Inf),
  labels = c("FDR < 0.05", "0.05 < FDR < 0.1", "Not Sig"),
  right = FALSE
)

merged$tcga_sig <- factor(merged$tcga_sig, levels = c("FDR < 0.05", "0.05 < FDR < 0.1", "Not Sig"))

fit <- cor.test(merged$r, merged$Coef, method = 'spearman')
cor_res[i] <- fit$estimate

pdf(file=file.path(dir_out, 'Fig/scatter', paste(selected_drug[i], "scatter.pdf", sep="_")), 
width = 4.5, height = 4.5)

p <- ggplot(merged, aes(x = r, y = Coef)) +
  geom_point(aes(color = tcga_sig), size = 1.5, alpha = 0.85) +
 # geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
 # geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
 # geom_smooth(method = "lm", se = FALSE, color = "steelblue", linetype = "dashed") +
  geom_text_repel(
    aes(label = ifelse(tcga_sig != "Not Sig", gene_name, "")),
    size = 2.5, max.overlaps = 15
  ) +
  scale_color_manual(
    name = "TCGA Significance",
    values = c(
      "FDR < 0.05" = "red",
      "0.05 < FDR < 0.1" = "orange",
      "Not Sig" = "gray80"
    )
  ) +
  labs(
    title = " ",
    x = "Pearson correlation (r)",
    y = "Cox Coefficient (logHR)"
  ) +
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="right",
        legend.text = element_text(size = 6, face="bold"),
        legend.title = element_blank()) 

print(p)

dev.off()


}

#################################################################
## KM plot (pan-STS) 
#################################################################
sig.meta <- res.meta[res.meta$padj < 0.05 & abs(res.meta$r) >= 0.3 & res.meta$drug %in% drug, ]
sig.cox <- meta.cox[meta.cox$padj < 0.05, ]

int <- intersect(sig.cox$Gene, sig.meta$gene_name)
intersect(sig.cox$Gene, genes)
# "DLX2"  "STMN1" "HMGA2"
intersect(sig.meta$gene_name, genes)
# "ATRX" "GLI1" "SMO"  "ASPM" "DLX2" "PTEN"
intersect(intersect(sig.cox$Gene, sig.meta$gene_name), genes)
# [1] "DLX2"

gene <- c(intersect(intersect(sig.cox$Gene, sig.meta$gene_name), genes), int)  

for(i in 1:length(gene)){

print(i)
clin_subset <- clin
expr_gene <- expr[rownames(expr) %in% gene[i], ]
expr_gene <- expr_gene[names(expr_gene) %in% rownames(clin_subset)]
clin_subset <- clin_subset[names(expr_gene), ]  
clin_subset$gene_expr <- as.numeric(expr_gene)

# Create quartiles
qtiles <- quantile(clin_subset$gene_expr, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
qtiles_unique <- unique(qtiles)

if (length(qtiles_unique) == 5) {
  # All quantiles are unique: assign quartiles safely
  clin_subset$quartile <- cut(
    clin_subset$gene_expr,
    breaks = qtiles,
    include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4")
  )


clin_sub <- clin_subset[clin_subset$quartile %in% c("Q1", "Q4"), ]
clin_sub$quartile <- droplevels(clin_sub$quartile)
clin_sub$quartile <- factor(clin_sub$quartile, levels = c("Q1", "Q4"))

# Survival object
surv_obj <- Surv(time = clin_sub$os_time, event = clin_sub$os_event)

# Fit survival curves and cox model
fit <- survfit(surv_obj ~ quartile, data = clin_sub)

# Generate KM plot
km_plot <- ggsurvplot(
  fit,
  data = clin_sub,
  risk.table = FALSE,
  pval = TRUE,
  pval.method = TRUE,
  pval.size = 2.5,
  conf.int = FALSE,
  palette = c("#0072B2", "#D55E00"),
  legend.title = gene[i],
  legend.labs = c("Low", "High"),
  xlab = "Time (months)",
  ylab = "Overall Survival",
  title = '',
  font.title = c(10, "bold"),
  font.x = c(10, "plain"),
  font.y = c(10, "plain"),
  font.legend = c(8),
  font.tickslab = 8,
  risk.table.height = 0.25,
  size = 0.5,         
  linetype = 1,        
  ggtheme = theme_classic()
)


pdf(file=file.path(dir_out, 'Fig/KM/pan', paste(gene[i], "km_tcga.pdf", sep="_")), 
width = 2.5, height =3)

# Print plot
print(km_plot)

dev.off()

   } else { 
  print("Not enough unique quantile cut points. Quartile classification skipped.")
  }

}

############################################
# Violin figures 
############################################
clin$histological <- ifelse(clin$histological == 'UPS-MFH', 'UPS', clin$histological)
clin$histological <- ifelse(clin$histological == 'synovial-mpnst', 'SS-MPNST', clin$histological)

for(j in 1:length(gene)){

df <- clin[, c('patientid', 'histological')]
expr_gene <- expr[rownames(expr) %in% gene[j], ]
expr_gene <- expr_gene[names(expr_gene) %in% rownames(df)]
df <- df[names(expr_gene), ]  
df$gene_expr <- as.numeric(expr_gene)

fit <- kruskal.test(df$gene_expr ~ df$histological)
p_value <- fit$p.value
p_label <- ifelse(p_value < 0.001,
                  paste0("p = ", formatC(p_value, format = "e", digits = 2)),
                  paste0("p = ", signif(p_value, 3)))

top <- levels(reorder(df$histological, -df$gene_expr))[5]
rng <- range(df$gene_expr, na.rm = TRUE)
y_annot <- rng[2] + 0.02 * diff(rng) 

pdf(file= file.path(dir_out, 'Fig/boxplot', paste(gene[j], "subtype_tcga.pdf", sep="-")), 
     width = 3, height = 3)

p <- ggplot(df, aes(x=reorder(histological, -gene_expr), y=gene_expr)) + 
  geom_boxplot(width = 0.5, fill= "#D8D8D8FF", color = "#181830FF") +
  coord_flip() + 
  #ylim(c(-0.5, 1)) +
  xlab("") +
  ylab(" ") +
  annotate("text", x = top, y =y_annot, label = p_label, hjust = 1, vjust = 1, size = 3) + 
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
        legend.text = element_text(size = 3, face="bold"),
        legend.title = element_blank()) 

print(p)

dev.off()

}

#############################################################
## KM plot (LMS subtype)
#############################################################
# subtype association results
sig.meta <- res.meta[res.meta$padj < 0.05 & abs(res.meta$r) >= 0.3 & res.meta$drug %in% drug, ]

dat <- read.csv(file.path(dir_in_tcga, 'cox_tcga_histo.csv'))
dat <- dat[dat$padj < 0.05, ]
dat <- dat[order(dat$padj, decreasing=FALSE), ]

int <- intersect(sig.meta$gene_name, dat$gene_name)
intersect(dat$gene_name, genes)
#[1] "ATRX"  "RB1"   "FLT3"  "HMGA2" "MYOCD"
intersect(int, genes)
#[1] "ATRX" 

gene <- c(intersect(int, genes), int)  # 230 

# dat <- dat[dat$gene_name %in% gene, ]
# Note: mainly focus on LMS subtype

# Extract expression vector and match to clinical data

for(i in 1:length(gene)){

print(i)

clin_subset <- clin[clin$histological == 'LMS', ]
expr_gene <- expr[rownames(expr) %in% gene[i], ]
expr_gene <- expr_gene[names(expr_gene) %in% rownames(clin_subset)]
clin_subset <- clin_subset[names(expr_gene), ]  
clin_subset$gene_expr <- as.numeric(expr_gene)

# Create quartiles
qtiles <- quantile(clin_subset$gene_expr, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
qtiles_unique <- unique(qtiles)

if (length(qtiles_unique) == 5) {
  # All quantiles are unique: assign quartiles safely
  clin_subset$quartile <- cut(
    clin_subset$gene_expr,
    breaks = qtiles,
    include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4")
  )


clin_sub <- clin_subset[clin_subset$quartile %in% c("Q1", "Q4"), ]
clin_sub$quartile <- droplevels(clin_sub$quartile)
clin_sub$quartile <- factor(clin_sub$quartile, levels = c("Q1", "Q4"))

# Survival object
surv_obj <- Surv(time = clin_sub$os_time, event = clin_sub$os_event)

# Fit survival curves and cox model
fit <- survfit(surv_obj ~ quartile, data = clin_sub)

# Generate KM plot
km_plot <- ggsurvplot(
  fit,
  data = clin_sub,
  risk.table = FALSE,
  pval = TRUE,
  pval.method = TRUE,
  pval.size = 2.5,
  conf.int = FALSE,
  palette = c("#0072B2", "#D55E00"),
  legend.title = gene[i],
  legend.labs = c("Low", "High"),
  xlab = "Time (months)",
  ylab = "Overall Survival",
  title = '',
  font.title = c(10, "bold"),
  font.x = c(10, "plain"),
  font.y = c(10, "plain"),
  font.legend = c(8),
  font.tickslab = 8,
  risk.table.height = 0.25,
  size = 0.5,         
  linetype = 1,        
  ggtheme = theme_classic()
)


pdf(file=file.path(dir_out, 'Fig/KM/histo/LMS', paste(gene[i], "km_tcga.pdf", sep="_")), 
width = 2.5, height =3)

# Print plot
print(km_plot)

dev.off()

   } else { 
  print("Not enough unique quantile cut points. Quartile classification skipped.")
  }

}



#############################################################
## KM plot (DDLPS subtype)
#############################################################
# subtype association results
sig.meta <- res.meta[res.meta$padj < 0.05 & abs(res.meta$r) >= 0.3 & res.meta$drug %in% drug, ]

dat <- read.csv(file.path(dir_in_tcga, 'cox_tcga_histo.csv'))
dat <- dat[dat$padj < 0.1, ]
dat <- dat[order(dat$padj, decreasing=FALSE), ]

int <- intersect(sig.meta$gene_name, dat$gene_name)
intersect(dat$gene_name, genes)
#[1] "ATRX"   "RB1"    "FLT3"   "HMGA2"  "MYOCD"  "PRKAA2" "STMN1"  "RIC3"  
#[9] "DLX2"
intersect(int, genes)
#[1] "ATRX" "DLX2"

gene <- c(intersect(int, genes), int)  # 511 

# dat <- dat[dat$gene_name %in% gene, ]
# Note: mainly focus on DDLPS subtype

# Extract expression vector and match to clinical data

for(i in 1:length(gene)){

print(i)

clin_subset <- clin[clin$histological == 'DDLPS', ]
expr_gene <- expr[rownames(expr) %in% gene[i], ]
expr_gene <- expr_gene[names(expr_gene) %in% rownames(clin_subset)]
clin_subset <- clin_subset[names(expr_gene), ]  
clin_subset$gene_expr <- as.numeric(expr_gene)

# Create quartiles
qtiles <- quantile(clin_subset$gene_expr, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
qtiles_unique <- unique(qtiles)

if (length(qtiles_unique) == 5) {
  # All quantiles are unique: assign quartiles safely
  clin_subset$quartile <- cut(
    clin_subset$gene_expr,
    breaks = qtiles,
    include.lowest = TRUE,
    labels = c("Q1", "Q2", "Q3", "Q4")
  )


clin_sub <- clin_subset[clin_subset$quartile %in% c("Q1", "Q4"), ]
clin_sub$quartile <- droplevels(clin_sub$quartile)
clin_sub$quartile <- factor(clin_sub$quartile, levels = c("Q1", "Q4"))

# Survival object
surv_obj <- Surv(time = clin_sub$os_time, event = clin_sub$os_event)

# Fit survival curves and cox model
fit <- survfit(surv_obj ~ quartile, data = clin_sub)

# Generate KM plot
km_plot <- ggsurvplot(
  fit,
  data = clin_sub,
  risk.table = FALSE,
  pval = TRUE,
  pval.method = TRUE,
  pval.size = 2.5,
  conf.int = FALSE,
  palette = c("#0072B2", "#D55E00"),
  legend.title = gene[i],
  legend.labs = c("Low", "High"),
  xlab = "Time (months)",
  ylab = "Overall Survival",
  title = '',
  font.title = c(10, "bold"),
  font.x = c(10, "plain"),
  font.y = c(10, "plain"),
  font.legend = c(8),
  font.tickslab = 8,
  risk.table.height = 0.25,
  size = 0.5,         
  linetype = 1,        
  ggtheme = theme_classic()
)


pdf(file=file.path(dir_out, 'Fig/KM/histo/DDLPS', paste(gene[i], "km_tcga.pdf", sep="_")), 
width = 2.5, height =3)

# Print plot
print(km_plot)

dev.off()

   } else { 
  print("Not enough unique quantile cut points. Quartile classification skipped.")
  }

}

########################################################################
## Volcano plot (gene level)
########################################################################
sig.meta <- res.meta[res.meta$padj < 0.05 & abs(res.meta$r) >= 0.3 & res.meta$drug %in% drug, ]
volcano.meta.cox <- meta.cox
volcano.meta.cox$gene <- volcano.meta.cox$Gene %in% sig.meta$gene_name
volcano.meta.cox$FDR <- cut(
  volcano.meta.cox$padj,
  breaks = c(-Inf, 0.05, 0.1, Inf),
  labels = c("FDR < 0.05", "0.05 < FDR < 0.1", "NS"),
  right = TRUE  # include right endpoint
)

volcano.meta.cox <- volcano.meta.cox %>%
  arrange(padj) %>%
  mutate(label = ifelse(row_number() <= 45, Gene, ""))

fdr_colors <- c("FDR < 0.05" = "#377540ff", 
                "0.05 < FDR < 0.1" = "#c18569ff", 
                "NS" = "grey80")


pdf(file=file.path(dir_out, 'Fig', "volcano_tcga.pdf"), 
width = 3.5, height =4)

# Plot
p <- ggplot(volcano.meta.cox, aes(x = Coef, y = -log10(Pval))) +
  geom_point(aes(color = FDR, shape = gene), alpha = 0.7) +
  scale_color_manual(values = fdr_colors) +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16), labels = c("No", "Yes")) +
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "log(HR)",
    y = "-log10(p-value)",
    color = "Sig FDR",
    shape = "Gene"
  ) +
   geom_text_repel(
    aes(label = label),
    size = 1.8,
    max.overlaps = 15,
    min.segment.length = 0
  ) +
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 6, face="bold"),
        legend.title = element_blank()) 

print(p)

dev.off()

########################################################################
## Load TCGA and Cell Line results (pathway level)
########################################################################
selected_drug <- read.csv(file.path(dir_in_cl, 'selected_drugs.csv'))
res.tcga <- qread(file.path(dir_in_tcga, 'hallmark_gsea_pathway.qs'))
res.cl <- qread(file.path(dir_in_cl_pathway, 'hallmark_gsea_pathway_drug.qs'))

moa <- read.csv(file.path(dir_moa, 'moa_dat.csv'))
drug <- selected_drug[selected_drug$sts == "Yes", "drug"]
sig.tcga <- res.tcga[res.tcga$padj < thr_gsea, ]
sig.tcga <- sig.tcga[order(sig.tcga$padj, decreasing= FALSE), ]

sig.cl <- res.cl[res.cl$padj < thr_gsea & res.cl$drug %in% drug, ]
sig.cl <- sig.cl[order(sig.cl$padj, decreasing= FALSE), ]

int <- intersect(sig.tcga$pathway, sig.cl$pathway)
sig.tcga <- sig.tcga[sig.tcga$pathway %in% int, ]
sig.tcga <- as.data.frame(sig.tcga)
sig.cl <- sig.cl[sig.cl$pathway %in% int, ]

sig.cl$moa <- sapply(1:nrow(sig.cl), function(k){
moa[moa$treatmentid == sig.cl$drug[k], 'TARGET_PATHWAY_updated']
})

## Visualize 
res <- lapply(1:length(int), function(k){

 df <- sig.cl[sig.cl$pathway == int[k], ]
 data.frame(pathway = int[k],
            n_drug = nrow(df),
            NES = sig.tcga[sig.tcga$pathway == int[k], 'NES'],
            pval = sig.tcga[sig.tcga$pathway == int[k], 'pval'],
            padj = sig.tcga[sig.tcga$pathway == int[k], 'padj'],
            drug = paste(df$drug, collapse = "|"))

})

res <- do.call(rbind, res)
res <- res[order(res$n_drug, decreasing= TRUE), ]

res_tcga <- res %>%
filter(pathway %in% int) %>%
  mutate(
    NES = NES,                            
    FDR   = pmax(padj, 1e-300),
    mlp   = -log10(FDR),     
    direction = ifelse(NES > 0, "Protective", "Risk"),
    label = tools::toTitleCase(tolower(substr(res$pathway, 10, nchar(res$pathway))))
  ) 


# y-threshold for FDR < 0.05
thr_y <- -log10(0.05)

pdf(file= file.path(dir_out, 'Fig', paste("volcano_hallmark_tcga.pdf", sep="-")), 
     width = 5, height = 5)

p <- ggplot(res_tcga, aes(x = NES, y = mlp)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.6, color = "grey55") +
  geom_vline(xintercept = 0,            linetype = 3, linewidth = 0.5, color = "grey65") +
  geom_point(aes(color = direction, size = n_drug), alpha = 0.9) +
  scale_x_continuous("NES") +
  scale_y_continuous(expression(-log[10]("FDR")),
                     breaks = c(2, 4, 6, 8, 10)) +
  scale_color_manual(values = c("Risk" = "#C44E52", "Protective" = "#4C72B0")) +
  scale_size_continuous(name = "drug", range = c(2.2, 7), breaks = pretty_breaks(4)) +
  labs(title = "") +
  theme(axis.text.x=element_text(size=8),
        axis.title=element_text(size=8),
        axis.text.y=element_text(size=8),
        strip.text = element_text(size=8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="right",
        legend.text = element_text(size = 5)) 

  lab_df <- res_tcga %>% slice_min(FDR, n = 23)
  p + ggrepel::geom_text_repel(data = lab_df, aes(label = sub("^HALLMARK_", "", pathway)),
                             size = 2, max.overlaps = Inf, box.padding = 0.3,
                             segment.size = 0.3, show.legend = FALSE)
print(p)

dev.off()












