#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Spodoptera frugiperda Sf9 cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of figS5 and figSx of the manuscript.

### Load necessary packages
library("DESeq2")
library("tidyverse")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("ggplotify")

#Load data
load("./functional_enrichment_analysis/intermediate_data/combined_enrichment_results.rData")
load("./functional_enrichment_analysis/intermediate_data/combined_gsea_results.rData")
load("./functional_enrichment_analysis/intermediate_data/protein_folding.rData")
load("./functional_enrichment_analysis/intermediate_data/cell_division.rData")
load("./functional_enrichment_analysis/intermediate_data/mitotic_cc.rData")
load("./functional_enrichment_analysis/intermediate_data/virus_response_cc.rData")

#Define figure directory
figsdir <- "./functional_enrichment_analysis/figures"

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

#Define labels for plot
labels_new <- c("Inf.24hpi vs\n Non-Inf.24hpi", "Inf.48hpi vs\n Inf.24hpi")
names(labels_new) <- c("24hpi", "48hpi")

ora_data <- combined_enrichment_results %>% filter(V2=="biological_process")

#Make ORA plot
barplot_ora_plot <- ggplot(data = ora_data, aes(x = ratio, y = fct_reorder(V3, ratio), fill = p.adjust)) + 
  geom_col()+
  scale_fill_distiller(palette ="Blues", direction = -1)+
  ylab("Enriched terms")+
  xlab("Enrichment Ratio")+
  facet_wrap(~day, labeller =labeller(day = labels_new))+
  labs(fill = "Adjusted p-value")+
  theme_script()+
  theme(axis.title = element_text(size = rel(1.6)),
        legend.position = "bottom",
        aspect.ratio = 2,
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.4)),
        axis.text = element_text(size = rel(1.3)),
        strip.text.x = element_text(size = rel(1.9), face = "bold"),
        strip.background = element_blank(),
        legend.key.width = unit(1, "cm"))

barplot_ora_plot

ggsave(filename= "/barplot_ORA_new.tiff", plot= barplot_ora_plot, path = figsdir, device="tiff", height=7, width=9, dpi = 600)

####Make Heatmaps using pheatmap
#Load apoptosis reference genes and data
apoptosis_az_paper <- read.delim("C:/Users/Nikolaus Virgolini/OneDrive - iBET/iBET/PhD_NV/Experiments/Transcriptomics/200921_BulkSeq#2/Sequencing/afternewdata/test3_210721/apoptosis_genes_azidrazin_paper_noNA.txt", header=T)
load("./DESeq_analysis/intermediate_results/DESeq_data_overall_IC.rData")

#Perform vst normalization
vst_overall <- vst(DESeq_data_overall_IC, blind = T)

# Prepare data and plot for Apoptosis
mat_ap  <- assay(vst_overall)[apoptosis_az_paper$locus,]
mat_ap_scaled <- t(scale(t(mat_ap)))
mat_ap_new <- na.omit(mat_ap_scaled)
Sample <- c("inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_48hpi", "non_48hpi", "non_48hpi", "non_48hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi")
df <- as.data.frame(Sample)
rownames(df) <- colnames(mat_ap_new)
mat_ap_plot <- mat_ap_new[,c(17,18,19,20,9:16,1:8)]
mycolors <- list(Sample = c("inf_24hpi" = "#F04949", "inf_48hpi" = "#B30C0C",
                            "pre_0hpi"="#4DAF4A", "non_24hpi" ="#84D5DE", "non_48hpi" = "#1868AA"))

z_score_heatmap_azidrazin <- as.ggplot(pheatmap(mat_ap_plot, cluster_cols = F,col = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(255), cluster_rows = T, show_rownames = F, show_colnames = F, annotation_col = df, annotation_names_col = F, fontsize_row = 5, border_col = NA, annotation_colors =  mycolors, treeheight_row = 0, legend=T, annotation_legend = F, gaps_col = c(4,8,12,16)))
heatmap_azidrazin_final <- z_score_heatmap_azidrazin + ggtitle("Cell apoptosis associated genes") + theme(plot.title = element_text(hjust=0.5, size = 16, face = "bold"))

heatmap_azidrazin_final

ggsave(filename= "/heatmap_azidrazin_final.tiff", plot= heatmap_azidrazin_final, path = figsdir, device="tiff", dpi = 600, height = 5, width = 5)

# Prepare data and plot for protein folding
mat_protein_folding  <- assay(vst_overall)[protein_folding,]
mat_protein_folding_scaled <- t(scale(t(mat_protein_folding)))
mat_protein_folding_new <- na.omit(mat_protein_folding_scaled)
Sample <- c("inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_48hpi", "non_48hpi", "non_48hpi", "non_48hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi")
df <- as.data.frame(Sample)
rownames(df) <- colnames(mat_protein_folding_new)
mat_protein_folding_plot <- mat_protein_folding_new[,c(17,18,19,20,9:16,1:8)]

z_score_heatmap_protein_folding <- as.ggplot(pheatmap(mat_protein_folding_plot,col = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(255), cluster_cols = F,  cluster_rows = T, show_rownames = F, show_colnames = F, annotation_col = df, annotation_names_col = F, fontsize_row = 5, border_col = NA, annotation_colors =  mycolors, treeheight_row = 0, legend=T, annotation_legend = F, gaps_col = c(4,8,12,16)))

heatmap_protein_folding_final <- z_score_heatmap_protein_folding + ggtitle("Protein folding (GO:0006457)") + theme(plot.title = element_text(hjust=0.5, size = 16, face = "bold"))

ggsave(filename= "/z_score_heatmap_protein_folding.tiff", plot= heatmap_protein_folding_final, path = figsdir, device="tiff", height = 5, width = 5)

# Prepare data and plot for cell division
mat_cell_division  <- assay(vst_overall)[cell_division,]
mat_cell_division_scaled <- t(scale(t(mat_cell_division)))
mat_cell_division_new <- na.omit(mat_cell_division_scaled)
Sample <- c("inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_48hpi", "non_48hpi", "non_48hpi", "non_48hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi")
df <- as.data.frame(Sample)
rownames(df) <- colnames(mat_cell_division_new)
mat_cell_division_plot <- mat_cell_division_new[,c(17,18,19,20,9:16,1:8)]

z_score_heatmap_cell_division <- as.ggplot(pheatmap(mat_cell_division_plot,col = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(255), cluster_cols = F,  cluster_rows = T, show_rownames = F, show_colnames = F, annotation_col = df, annotation_names_col = F, fontsize_row = 5, border_col = NA, annotation_colors =  mycolors, treeheight_row = 0, legend=T, annotation_legend = F, gaps_col = c(4,8,12,16)))

heatmap_cell_division_final <- z_score_heatmap_cell_division + ggtitle("Cell division (GO:0051301)") + theme(plot.title = element_text(hjust=0.5, size = 16, face = "bold"))

ggsave(filename= "/z_score_heatmap_cell_division.tiff", plot= heatmap_cell_division_final, path = figsdir, device="tiff", height = 5, width = 5)

# Prepare data and plot for mitotic cell cycle
mat_mitotic_cc  <- assay(vst_overall)[mitotic_cc,]
mat_mitotic_cc_scaled <- t(scale(t(mat_mitotic_cc)))
mat_mitotic_cc_new <- na.omit(mat_mitotic_cc_scaled)
Sample <- c("inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_48hpi", "non_48hpi", "non_48hpi", "non_48hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi")
df <- as.data.frame(Sample)
rownames(df) <- colnames(mat_mitotic_cc_new)
mat_mitotic_cc_plot <- mat_mitotic_cc_new[,c(17,18,19,20,9:16,1:8)]


z_score_heatmap_mitotic_cc <- as.ggplot(pheatmap(mat_mitotic_cc_plot, cluster_cols = F,col = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(255),   cluster_rows = T, show_rownames = F, show_colnames = F, annotation_col = df, annotation_names_col = F, fontsize_row = 5, border_col = NA, annotation_colors =  mycolors, treeheight_row = 0, legend=T, annotation_legend = F, annotation_legend_position = "bottom", gaps_col = c(4,8,12,16)))
heatmap_mitotic_cc_final <- z_score_heatmap_mitotic_cc + ggtitle("Mitotic cell cycle (GO:0051301)") + theme(plot.title = element_text(hjust=0.5, size = 16, face = "bold"))

ggsave(filename= "/real_z_score_heatmap_mitotic_cc.tiff", plot= heatmap_mitotic_cc_final, path = figsdir, device="tiff", height = 5, width = 5)


# Prepare data and plot for viral response
mat_viral_response_cc  <- assay(vst_overall)[virus_response_cc,]
mat_viral_response_cc_scaled <- t(scale(t(mat_viral_response_cc)))
mat_viral_response_cc_new <- na.omit(mat_viral_response_cc_scaled)
Sample <- c("inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_24hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "inf_48hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_24hpi", "non_48hpi", "non_48hpi", "non_48hpi", "non_48hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi", "pre_0hpi")
df <- as.data.frame(Sample)
rownames(df) <- colnames(mat_viral_response_cc_new)
mat_viral_response_cc_plot <- mat_viral_response_cc_new[,c(17,18,19,20,9:16,1:8)]


z_score_heatmap_viral_response_cc <- as.ggplot(pheatmap(mat_viral_response_cc_plot, cluster_cols = F,col = colorRampPalette( rev(brewer.pal(9, "RdYlBu")) )(255),   cluster_rows = T, show_rownames = F, show_colnames = F, annotation_col = df, annotation_names_col = F, fontsize_row = 5, border_col = NA, annotation_colors =  mycolors, treeheight_row = 0, legend=T, annotation_legend = F, annotation_legend_position = "bottom", gaps_col = c(4,8,12,16)))
heatmap_viral_response_cc_final <- z_score_heatmap_viral_response_cc + ggtitle("Defense response to virus (GO:0051607)") + theme(plot.title = element_text(hjust=0.5, size = 16, face = "bold"))

ggsave(filename= "/real_z_score_heatmap_viral_response_cc.tiff", plot= heatmap_viral_response_cc_final, path = figsdir, device="tiff", height = 5, width = 5)

#############End of script

