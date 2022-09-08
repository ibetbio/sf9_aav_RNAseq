#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Spodoptera frugiperda Sf9 cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of figure 3 of the manuscript

##Load necessary packages
library("DESeq2")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("gplots")
library("ggplotify")
library("ggpubr")
library("RColorBrewer")
library("pheatmap")

#Define folders
result_dir <- "./DESeq_analysis/deseq_results"
intermdir <- "./DESeq_analysis/intermediate_results"
figsdir <- "./DESeq_analysis/figures"

# Set plotting parameters
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

######################################################################################
################# Make overall QC plots - fig 2b and 2c ##############################
######################################################################################
# Load necessary data into R script
load(paste(intermdir, "/DESeq_data_overall_IC.rData", sep=""))

# Run vst-normalization on Insect cell data
vst_overall <- vst(DESeq_data_overall_IC, blind = T)
sampleDists_overall <- dist(t(assay(vst_overall)))

# Make euclidean sample distance matrix
sampleDistMatrix <- as.matrix(sampleDists_overall)
rownames(sampleDistMatrix) <- c("inf-24hpi_R01", "inf-24hpi_R02", "inf-24hpi_R03", "inf-24hpi_R04",
                                "inf-48hpi_R01", "inf-48hpi_R02", "inf-48hpi_R03", "inf-48hpi_R04",
                                "non-24hpi_R01", "non-24hpi_R02", "non-24hpi_R03", "non-24hpi_R04",
                                "non-48hpi_R01", "non-48hpi_R02", "non-48hpi_R03", "non-48hpi_R04",
                                "inoc-0hpi_R01", "inoc-0hpi_R02", "inoc-0hpi_R03", "inoc-0hpi_R04")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)

eucl_distance_plot_overall <- pheatmap(sampleDistMatrix,
                                       clustering_distance_rows = sampleDists_overall,
                                       clustering_distance_cols = sampleDists_overall,
                                       col = colors, treeheight_col = 0, fontsize = 12)
eucl_distance_plot_overall

ggsave(filename= "/eucl_distance_plot_overall.tiff", plot= eucl_distance_plot_overall, path = figsdir, device="tiff", width=7, height=5)

pcaData_overall <- plotPCA(vst_overall, intgroup = c( "condition"), returnData = TRUE)
pcaData_overall
percentVar_overall <- round(100 * attr(pcaData_overall, "percentVar"))

pca_plot_overall <- ggplot(pcaData_overall, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =6, alpha=0.8, shape = c(16,16,16,16, 15,15,15,15, 16,16,16,16, 15,15,15,15, 18,18,18,18))+
  xlab(paste0("PC1: ", percentVar_overall[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_overall[2], "% variance")) +
  coord_fixed() +
  theme_script()+
  scale_y_continuous(limits=c(-15,15))+
  scale_x_continuous(limits=c(-25,75))+
  theme(legend.position = "bottom",
        legend.title = element_blank())+
  #scale_color_okabe_ito(order = c(2, 5, 1, 6, 3), labels=c("inf_24hpi"="Infected - 24 hpi", "inf_48hpi"="Infected - 48 hpi", "non_24hpi" = "Non-Infected - 24 hpi", "non_48hpi" = "Non-Infected - 48 hpi", "pre_0hpi" = "Pre-Culture - 0 hpi"))+
  scale_color_manual(labels=c("inf_24hpi"="Infected - 24 hpi", "inf_48hpi"="Infected - 48 hpi", "non_24hpi" = "Non-infected - 24 hpi", "non_48hpi" = "Non-infected - 48 hpi", "seed_0hpi" = "Pre-infection - 0 hpi"), values = c(  "#F04949", "#B30C0C", "#84D5DE", "#1868AA", "#4DAF4A"))+
  theme(legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=rel(2)),
        axis.text = element_text(size=rel(2)),
        axis.title = element_text(size=rel(2.2)),
        aspect.ratio = 0.5)+
  guides(color=guide_legend(nrow=2, byrow=TRUE, override.aes=list(shape=c(16, 15, 16, 15, 18))))

pca_plot_overall

ggsave(filename= "/pca_plot_overall.tiff", plot= pca_plot_overall, path = figsdir, device="tiff",  width=7, height=5, dpi = 600)

########## Make Volcano plots########
load(paste(intermdir, "/res_infvsnon24.rData", sep=""))

res_infvsnon24 <- res_infvsnon24[res_infvsnon24$baseMean >=50,]
res_infvsnon24$threshold_DE <- res_infvsnon24$padj < 0.05 & abs(res_infvsnon24$log2FoldChange) >= 0.5849625

DE_event24 <- data.frame(geneid = rownames(res_infvsnon24),
                         log2FoldChange = res_infvsnon24$log2FoldChange,
                         padj = res_infvsnon24$padj,
                         threshold = res_infvsnon24$threshold_DE)

DE_event24 <- DE_event24 %>% mutate(pointcolor = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "#DB2763",
  log2FoldChange < 0 & threshold == 1 ~ "#29BF12",
  threshold == 0 ~ "gray"
))
DE_event24 <- DE_event24 %>% mutate(pointclass = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "Upregulated",
  log2FoldChange < 0 & threshold == 1 ~ "Downregulated",
  threshold == 0 ~ "Not Significant"
))

DE_event24 <- DE_event24 %>% mutate(pointsize = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "0.5",
  log2FoldChange < 0 & threshold == 1 ~ "0.5",
  threshold == 0 ~ "0.2"
))


DE_event24 <- na.omit(DE_event24)

volcano_infvsnon24 <- ggplot(data=DE_event24) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), color=pointcolor, fill=pointclass), alpha=0.4, size=1) +
  scale_colour_manual(
    values = c("#DB2763","#29BF12", "gray"),
    labels = c("Significantly downregulated" ," Significantly upregulated", "Not significant"))+theme_script() +
  theme(legend.position = "none", legend.title = element_blank()) +
  xlab(bquote(bold(~log[2] ~ "fold change"))) + ylab(bquote(bold(~-log[10] ~ "BH adjusted p-value")))+
  ggtitle("Early infection")+
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
        legend.text = element_text(size=rel(2)),
        axis.text = element_text(size=rel(2)),
        axis.title = element_text(size=rel(2.2)),
        plot.title = element_text(size = rel(2.2), hjust = 0.5, face = "bold"),
        legend.spacing.x = unit(0.5, 'cm')) +  
  guides(fill="none", colour = guide_legend(override.aes = list(size = 3)))+
  scale_x_continuous(limits = c(-6,6))

volcano_infvsnon24
sum(DE_event24$pointclass == "Downregulated", na.rm=T)
sum(DE_event24$pointclass == "Upregulated", na.rm=T)

ggsave(filename= "/volcano_infvsnon24.tiff", plot= volcano_infvsnon24, path = figsdir, device="tiff", height=5, width=7, dpi = 600)

####Volcano plot 48 vs 24
load(paste(intermdir, "/res_inf48vs24.rData", sep=""))

res_inf48vs24 <- res_inf48vs24[res_inf48vs24$baseMean >=50,]
res_inf48vs24$threshold_DE <- res_inf48vs24$padj < 0.05 & abs(res_inf48vs24$log2FoldChange) >= 0.5849625

DE_event48 <- data.frame(geneid = rownames(res_inf48vs24),
                         log2FoldChange = res_inf48vs24$log2FoldChange,
                         padj = res_inf48vs24$padj,
                         threshold = res_inf48vs24$threshold_DE)

DE_event48 <- DE_event48 %>% mutate(pointcolor = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "#DB2763",
  log2FoldChange < 0 & threshold == 1 ~ "#29BF12",
  threshold == 0 ~ "gray"
))
DE_event48 <- DE_event48 %>% mutate(pointclass = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "Upregulated",
  log2FoldChange < 0 & threshold == 1 ~ "Downregulated",
  threshold == 0 ~ "Not Significant"
))

DE_event48 <- DE_event48 %>% mutate(pointsize = case_when(
  log2FoldChange > 0 & threshold == 1 ~ "0.5",
  log2FoldChange < 0 & threshold == 1 ~ "0.5",
  threshold == 0 ~ "0.2"
))


DE_event48 <- na.omit(DE_event48)

volcano_inf_48vs24 <- ggplot(data=DE_event48) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), color=pointcolor, fill=pointclass), alpha=0.4, size=0.7) +
  scale_colour_manual(
    values = c("#DB2763","#29BF12", "gray"),
    labels = c("Significantly downregulated" ," Significantly upregulated", "Not significant"))+
  theme_script() +
  theme(legend.position = "none", legend.title = element_blank()) +
  xlab(bquote(bold(~log[2] ~ "fold change"))) + ylab(bquote(bold(~-log[10] ~ "BH adjusted p-value")))+
  ggtitle(("Late infection"))+
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
        legend.text = element_text(size=rel(2)),
        axis.text = element_text(size=rel(2)),
        axis.title = element_text(size=rel(2.2)),
        plot.title = element_text(size = rel(2.2), hjust = 0.5, face = "bold"),
        legend.spacing.x = unit(0.5, 'cm')) + 
  guides(fill="none", colour = guide_legend(override.aes = list(size = 3)))+
  scale_x_continuous(limits = c(-10,10))

volcano_inf_48vs24
sum(DE_event48$pointclass == "Downregulated", na.rm=T)
sum(DE_event48$pointclass == "Upregulated", na.rm=T)

ggsave(filename= "/volcano_inf_48vs24.tiff", plot= volcano_inf_48vs24, path = figsdir, device="tiff", height=5, width=7, dpi = 600)

# Combined plot for figure 2
combined4 <- ggarrange("",
                       ggarrange("", 
                                 pca_plot_overall,
                                 "",
                                 ncol = 3,  widths = c(0.01, 1.5, 0.01)), 
                       "",
                       ggarrange("", 
                                 volcano_infvsnon24,
                                 "", 
                                 volcano_inf_48vs24,
                                 "", ncol = 5, widths = c(0.5, 1, 0.01, 1, 0.5), hjust = 1, vjust=-0.5, align = "hv", common.legend = T, legend = "bottom"), 
                       
                       "", 
                       nrow = 5, heights = c(0.01, 1.1, 0.1, 1.3, 0.01))

combined4

ggsave(filename= "./comb4.tiff", plot= combined4, path = figsdir, device="tiff", width=16, height=12, dpi = 600)

###End of script
