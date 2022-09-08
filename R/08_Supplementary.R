#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Sf9 insect cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of the supplementary figure 3 of the manuscript.

## Load necessary packages
library("dplyr")
library("ggplot2")
library("tidyverse")
library("gplots")
library("ggplotify")
library("ggpubr")
library("RColorBrewer")
library("pheatmap")
library( "genefilter")
library("DESeq2")
library("ggplotify")

result_dir <- "./DESeq_analysis/deseq_results"
intermdir <- "./DESeq_analysis/intermediate_results"
figsdir <- "./supplements"

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


remove <- read.delim("DESeq_analysis/reference_files/genes_to_remove_deseq.txt", header = F)

##############################################################################################
############DESeq analysis of bacv genes WITHOUT cellular genes in the design#################
#################################### FigS3 ###################################################
##############################################################################################

###Load data
load(paste(intermdir, "/DESeq_data_bacv.rData", sep=""))
#Remove transgenes
vst_bacv <- varianceStabilizingTransformation(DESeq_data_bacv, blind = T)

transgenes <- c("rep78", "cap", "TN", "CMV", "ITR_L", "eGFP")
genesToRemove <- which(!rownames(vst_bacv) %in%  transgenes)
vst_bacv_new <- vst_bacv[genesToRemove,]

#Find top variable genes of all of viral transcriptome
topVarGenes <- head(order(rowVars(assay(vst_bacv_new)), decreasing=TRUE), 156)

# Make heatmap of viral genes
mat <- assay(vst_bacv_new)[topVarGenes, ]
mat_zscore <- t(scale(t(mat)))
mat_zscore_new <- na.omit(mat_zscore)
Sample <- c("Infected - 24 hpi", "Infected - 24 hpi", "Infected - 24 hpi", "Infected - 24 hpi", "Infected - 48 hpi", "Infected - 48 hpi", "Infected - 48 hpi", "Infected - 48 hpi")
df <- as.data.frame(Sample)
rownames(df) <- colnames(mat_zscore_new)
mat_zscore_plot <- mat_zscore_new
mycolors <- list(Sample = c("Infected - 24 hpi" ="#F04949", "Infected - 48 hpi" = "#B30C0C"))

z_score_heatmap_top100 <- as.ggplot(pheatmap(mat_zscore_plot,clustering_distance_rows = "correlation", scale = "row", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), cluster_cols = F,  cluster_rows = T, show_rownames = F, show_colnames = F, annotation_col = df, annotation_names_col = F, fontsize_row = 5, border_col = NA, annotation_colors =  mycolors, legend=T, annotation_legend = F))

ggsave(filename= "/z_score_top100_bacv.tiff", plot= z_score_heatmap_top100, path = figsdir, device="tiff", height = 11,width = 10, units = "cm", dpi = 600)

####################################################################################################
#######################DESeq analysis of bacv genes with cellular genes in the design######################
####################################################################################################
count_dir <- "./DESeq_analysis/raw_data/counts/"

# Define File names
count_file_names <- grep("count", list.files(count_dir), value=T)

cond_bacv_all <- c("24", "24", "24", "24",
                   "48", "48", "48", "48")


si_bacv_all <- data.frame(sampleName = count_file_names[c(1:8)],
                          fileName =  count_file_names[c(1:8)],
                          condition = cond_bacv_all)

DESeq_data_bacv_all <- DESeqDataSetFromHTSeqCount(
  sampleTable = si_bacv_all,
  directory = count_dir,
  design = ~ condition)

colData(DESeq_data_bacv_all)$condition <- factor(colData(DESeq_data_bacv_all)$condition,
                                                 levels = c("48", "24"))

DESeq_data_bacv_all$condition <- relevel(DESeq_data_bacv_all$condition, "24")

DESeq_data_bacv_all <- DESeq(DESeq_data_bacv_all)

vst_bacv_all <- varianceStabilizingTransformation(DESeq_data_bacv_all, blind = T)

remove <- read.delim("DESeq_analysis/reference_files/genes_to_remove_deseq.txt", header =F)

genesToRemove <- which(rownames(vst_bacv_all) %in%  remove$V1)
vst_bacv_all_new <- vst_bacv_all[genesToRemove,]

genesToRemove <- which(!rownames(vst_bacv_all_new) %in%  transgenes)
vst_bacv_final_all <- vst_bacv_all_new[genesToRemove,]
topVarGenes <- head(order(rowVars(assay(vst_bacv_final_all)), decreasing=TRUE), 156)

mat_all <- assay(vst_bacv_final_all)[ topVarGenes, ]
mat_all_zscore <- t(scale(t(mat_all)))
mat_all_zscore_new <- na.omit(mat_all_zscore)
Sample <- c("Infected - 24 hpi", "Infected - 24 hpi", "Infected - 24 hpi", "Infected - 24 hpi", "Infected - 48 hpi", "Infected - 48 hpi", "Infected - 48 hpi", "Infected - 48 hpi")
df <- as.data.frame(Sample)
rownames(df) <- colnames(mat_all_zscore_new)
mat_all_zscore_plot <- mat_all_zscore_new

z_score_heatmap_top_alldesign <- as.ggplot(pheatmap(mat_all_zscore_plot, scale = "row",clustering_distance_rows = "correlation", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255), cluster_cols = F,  cluster_rows = T, show_rownames = F, show_colnames = F, annotation_col = df, annotation_names_col = F, fontsize_row = 5, border_col = NA, annotation_colors =  mycolors, legend=T, annotation_legend = F))

ggsave(filename= "/z_score_top_bacv_design_all.tiff", plot= z_score_heatmap_top_alldesign, path = figsdir, device="tiff", height = 11,width = 10, units = "cm", dpi = 600)

#### End of script
