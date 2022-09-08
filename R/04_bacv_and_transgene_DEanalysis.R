#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Spodoptera frugiperda Sf9 cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of figure 3 of the manuscript

##Load necessary packages
library("DESeq2")
library("tidyverse")
library("dplyr")
library("pheatmap")
library("ggplot2")
library("xlsx")
library("RColorBrewer")
library("ggplotify")
library("ggpubr")

### Define folders
result_dir <- "./DESeq_analysis/deseq_results"
intermdir <- "./DESeq_analysis/intermediate_results"
figsdir <- "./DESeq_analysis/figures"

### Set plotting parameters
windowsFonts("Helvetica" = windowsFont("Helvetica"))
theme_script <- function(base_size = 10, base_family = "Helvetica")
{
  library(grid)
  (theme_bw(base_size = base_size, base_family = base_family) + 
      theme(text = element_text(color = "black"), 
            axis.title = element_text(face ="bold", size = rel(1)), 
            plot.title = element_text(size= rel(1), color = "black"), 
            legend.title = element_text (face = "bold"), 
            legend.background = element_rect(fill= "transparent"), 
            legend.key.size = unit(0.8, "lines"), 
            panel.border = element_rect(color ="black", size =1), 
            panel.grid = element_blank()))
}

##Load DESeq results infected 48 vs 24 hpi for baculovirus genes
load(paste(intermdir, "/DESeq_data_bacv.rData", sep=""))

## Vst normalize data and make euclidean distance matrix and PCA
vst <- varianceStabilizingTransformation(DESeq_data_bacv, blind = T)

#Make a sample distance matrix (not included in the manuscript)
sampleDists <- dist(t(assay(vst)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vst$time
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

bacv_eucld <- pheatmap(sampleDistMatrix,
                       clustering_distance_rows = sampleDists,
                       clustering_distance_cols = sampleDists,
                       col = colors)

### Perform principle component analysis
pcaData_bacv<- plotPCA(vst, intgroup = c( "condition"), returnData = TRUE)
pcaData_bacv
percentVar_bacv <- round(100 * attr(pcaData_bacv, "percentVar"))

### Plot PCA
pca_plot_bacv <- ggplot(pcaData_bacv, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =6, alpha=0.8, shape = c(16,16,16,16, 15,15,15,15))+
  xlab(paste0("PC1: ", percentVar_bacv[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_bacv[2], "% variance")) +
  coord_fixed() +
  theme_script()+
  scale_y_continuous(limits=c(-15,15))+
  scale_x_continuous(limits=c(-10,10))+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  #scale_color_okabe_ito(order = c(2, 5, 1, 6, 3), labels=c("inf_24hpi"="Infected - 24 hpi", "inf_48hpi"="Infected - 48 hpi", "non_24hpi" = "Non-Infected - 24 hpi", "non_48hpi" = "Non-Infected - 48 hpi", "pre_0hpi" = "Pre-Culture - 0 hpi"))+
  scale_color_manual(labels=c("24"="Infected - 24 hpi", "48"="Infected - 48 hpi"), values = c(  "#F04949", "#B30C0C"))+
  theme(legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.6)),
        axis.text = element_text(size=rel(1.6)),
        axis.title = element_text(size=rel(1.8)),
        aspect.ratio = 0.7)+
  guides(color=guide_legend(override.aes=list(shape=c(16, 15))))

pca_plot_bacv

### Save PCA plot
ggsave(filename= "/pca_plot_bacv.tiff", plot= pca_plot_bacv, path = figsdir, device="tiff",  width=7, height=5, dpi = 600)

############################################################################################
###########################Make volcano plot for DE viral genes##########################
############################################################################################

#Paste results of DESeq2 for bacv genes
res_bacv <- results(DESeq_data_bacv, independentFiltering = T, lfcThreshold = 0)

# Subset and order differentially expressed genes with adjusted p-value < 0.05 and |FC| > 1.5
de_res_bacv <- subset(res_bacv,
                      abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >= 50)
de_res_bacv <- de_res_bacv[order(de_res_bacv$log2FoldChange, decreasing = T),]

#Get annotation file
annotation_bacv <- read.table("./DESeq_analysis/reference_files/GCF_000838485.1_ViralProj14023_feature_table.txt", header=T, sep="\t")

#Annotate results file
de_res_bacv_df <- as.data.frame(de_res_bacv)
de_res_bacv_annotated <- dplyr::left_join(rownames_to_column(de_res_bacv_df), annotation_bacv %>% filter(feature =="CDS"), by = c("rowname" = "locus_tag"))

res_bacv_df <- as.data.frame(res_bacv)
res_bacv_annotated <- dplyr::left_join(rownames_to_column(res_bacv_df), annotation_bacv %>% filter(feature =="CDS"), by = c("rowname" = "locus_tag"))

# save the results to Excel
fn <- paste(result_dir, "/de_res_bacv_annotated.xlsx",sep="")
if (file.exists(fn)) {file.remove(fn)}
write.xlsx(de_res_bacv_annotated,
          file = fn,
          sheetName = "DE genes baculovirus and transgenes",
          row.names=F)

### save the results to Excel
fn <- paste(result_dir, "/res_bacv_annotated.xlsx",sep="")
if (file.exists(fn)) {file.remove(fn)}
write.xlsx(res_bacv_annotated,
           file = fn,
           sheetName = "All genes baculovirus and transgenes",
           row.names=F)

### Remove transgenes from analysis
transgenes <- c("rep78", "cap", "TN", "CMV", "ITR_L", "eGFP")
genesToRemove <- which(!rownames(res_bacv_df) %in%  transgenes)
res_bacv_final <- res_bacv_df[genesToRemove,]

### Manipulate data for volcano plot
de_res_bacv_volcano <- res_bacv_final[res_bacv_final$baseMean >=50,]
de_res_bacv_volcano$threshold_DE <- de_res_bacv_volcano$padj < 0.05 & abs(de_res_bacv_volcano$log2FoldChange) >= 0.5849625

DE_eventbacv <- cbind(de_res_bacv_volcano$log2FoldChange, de_res_bacv_volcano$padj, as.logical(de_res_bacv_volcano$threshold_DE))

DE_eventbacv <- data.frame(geneid = row.names(de_res_bacv_volcano),
                           log2FoldChange = de_res_bacv_volcano$log2FoldChange,
                           padj = de_res_bacv_volcano$padj,
                           threshold = de_res_bacv_volcano$threshold_DE)

DE_eventbacv <- DE_eventbacv %>% mutate(pointcolor = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "#29BF12",
  log2FoldChange < 0 & threshold == 1 ~ "#009ACD",
  threshold == 0 ~ "gray"
))
DE_eventbacv <- DE_eventbacv %>% mutate(pointclass = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "Upregulated",
  log2FoldChange < 0 & threshold == 1 ~ "Downregulated",
  threshold == 0 ~ "Not Significant"
))

DE_eventbacv <- DE_eventbacv %>% mutate(pointsize = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "0.5",
  log2FoldChange < 0 & threshold == 1 ~ "0.5",
  threshold == 0 ~ "0.2"
))

### In order for genes with padj = 0 to appear we assume a padj of 0.0000001 for those
DE_eventbacv <- DE_eventbacv %>% mutate(padj_plot = case_when(
  padj == 0 ~ paste0(padj+0.0000001),
  padj != 0 ~ paste0(padj)
))

### Plot volcano plot
volcanobacv <- ggplot(data=DE_eventbacv) +
  geom_point(aes(x=log2FoldChange, y=-log10(as.numeric(padj_plot)), color=pointcolor, fill=pointclass), alpha=0.4, size=2) +
  scale_colour_manual(
    values = c("#DB2763","#29BF12", "gray"),
    labels = c("Downregulated", "Upregulated", "Not significant"))+
  theme_script() +
  theme(legend.position = "none", legend.title = element_blank()) +
  xlab(bquote(bold(~Log[2] ~ "fold change"))) + ylab(bquote(bold(~-log[10] ~ "BH adjusted p-value")))+
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "gray"
  ) +
  geom_vline(
    xintercept = c(-0.5849625, 0.5849625),
    linetype = "dashed",
    color = "light gray"
  ) +
  theme(legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(1.6)),
        axis.text = element_text(size=rel(1.6)),
        axis.title = element_text(size=rel(1.8)),
        plot.title = element_text(size = rel(1.6), hjust = 0.5, face = "bold"),
        legend.spacing.x = unit(0.5, 'cm'),
        aspect.ratio = 0.7) + 
  guides(fill="none", colour = guide_legend(override.aes = list(size = 3)))+
  scale_x_continuous(limits = c(-4,4))

volcanobacv

### Save plot
ggsave(filename= "/bacv_volcano_annotated.tiff", plot= volcanobacv, path = figsdir, device="tiff", height=5, width=7, dpi=600)

sum(DE_eventbacv$pointclass == "Downregulated", na.rm=T)
sum(DE_eventbacv$pointclass == "Upregulated", na.rm=T)

##############################################################################################
###########################################Figure 6c###########################################
##############################################################################################
load(paste(intermdir, "/transgenes_norm_count.rData", sep=""))
labels_new <- c("Inf.24hpi", "Inf.48hpi")
names(labels_new) <- c("24", "48")

axis_text_y1 <- bquote(atop(bold('Normalized AAV transgene'),
                            bold('counts'~'('~10^{3}~')')))

transgenes <- ggplot(transgene_norm_count, aes(x=gene, y=count/1000, color = condition)) + 
  geom_boxplot() + 
  theme_script()+
  facet_grid(~condition, labeller =labeller(condition = labels_new))+
  xlab("Transgenes") +
  scale_y_continuous(limits = c(0, 300), name = axis_text_y1)+
  scale_x_discrete(labels = c("eGFP"="gfp", "rep" = "rep", "cap" = "cap"))+
  theme(axis.title = element_text(size = rel(1.8)),
        legend.position = "none",
        aspect.ratio = 1.8,
        legend.text = element_text(size = rel(1.8)),
        axis.text.x = element_text(size = rel(1.7), face = "italic"),
        axis.text.y = element_text(size = rel(1.7)),
        strip.text.x = element_text(size = rel(2), face = "bold"),
        strip.background = element_blank())

transgenes

ggsave(filename= "/transgene_expression.tiff", plot= transgenes, path = figsdir, device="tiff", height=7, width=7, dpi = 600)

rep_all <- transgene_norm_count %>% filter(gene == "rep")
rep_FC <- mean(rep_all$count[5:8])/mean(rep_all$count[1:4])
rep_FC


### Combine figures together
combined_fig3 <- ggarrange("",
                      pca_plot_bacv,
                       "",
                      volcanobacv,
                       "", 
                      transgenes,
                      "",
                       nrow = 7, heights = c(0.01, 1, 0.2, 1, 0.2, 1, 0.01))

combined_fig3

ggsave(filename= "./fig3.tiff", plot= combined_fig3, path = figsdir, device="tiff", width=8, height=12, dpi = 600)
