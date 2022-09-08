#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Spodoptera frugiperda Sf9 cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of the read distribution calculations and figure 2 of the manuscript.

### Load necessary packages
library(tidyverse)
library(ggpubr)
library(ggplot)

### Manipulate data for read distribution plot
#Get raw data
rrna_data <- read.table("./bioinformatic_statistics/raw_data/rrna_mapping_summary_header", sep="\t",header=T)
mapping_data <- read.table("./bioinformatic_statistics/raw_data/star_mapping_summary_header", sep="\t",header=T)
initial_data <- read.csv("./bioinformatic_statistics/raw_data/sequencing_statistics.csv", sep=",",header=T)

#Define replicates and timepoints. 
##Note: replicate R02 from Non-Infected is first because it has the largest data size and needed to be run first. Therefore it is mixed up here for this analysis.
condition <-c("Non-infected", rep("Infected", 8), rep("Non-infected",7), rep("Pre-infection",4))
replicate <- c("2", 
               "1", "2", "3", "4", 
               "1", "2", "3", "4", 
               "1", "3", "4",
               "1", "2", "3", "4", 
               "1", "2", "3", "4")

timepoint <- c("24 hpi", rep("24 hpi",4),rep("48 hpi",4),rep("24 hpi",3),rep("48 hpi",4),rep("0 hpi",4))               

#Define initial read number
initial_reads <- initial_data$reads
#Define trimmed reads
trimming <- initial_reads-rrna_data$input_reads
trimming_percent <- trimming/initial_reads*100
#Define rRNA removed reads
rrna_removed <- rrna_data$input_reads - mapping_data$input_reads
rrna_percent <- rrna_removed/initial_reads*100
#Define mapped reads
mapped <- mapping_data$uniquely_mapped
mapped_percent <- mapped/initial_reads*100

#Define unmapped reads
unmapped <- mapping_data$input_reads-mapping_data$uniquely_mapped
unmapped_percent <- unmapped/initial_reads*100

#Combine reads for plotting
read_distribution <- data.frame(trimming_percent, rrna_percent, unmapped_percent, mapped_percent)               

read_distribution_plot <- read_distribution  %>% gather(key="option", value="value")

read_distribution_plot_final <- data.frame(read_distribution_plot, condition, replicate, timepoint, sampleNames = mapping_data$sampleNames) 

replicates <- c("1", "2", "3", "4", 
                "1", "2", "3", "4", 
                "1", "2","3", "4",
                "1", "2", "3", "4", 
                "1", "2", "3", "4")

#Save data for plotting in R script specific for figure 2
save(read_distribution_plot_final, file = "./bioinformatic_statistics/intermediate_data/read_distribution_plot_final.rData")

#################Percentage of viral mapping calculation - viral content
# Get data and define row and column names
overall_inf24_1 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP01-R01.count", header=F)
colnames(overall_inf24_1) <- c("gene_id", "counts")
overall_inf24_2 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP01-R02.count", header=F)
colnames(overall_inf24_2) <- c("gene_id", "counts")
overall_inf24_3 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP01-R03.count", header=F)
colnames(overall_inf24_3) <- c("gene_id", "counts")
overall_inf24_4 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP01-R04.count", header=F)
colnames(overall_inf24_4) <- c("gene_id", "counts")

overall_inf48_1 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP02-R01.count", header=F)
colnames(overall_inf48_1) <- c("gene_id", "counts")
overall_inf48_2 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP02-R02.count", header=F)
colnames(overall_inf48_2) <- c("gene_id", "counts")
overall_inf48_3 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP02-R03.count", header=F)
colnames(overall_inf48_3) <- c("gene_id", "counts")
overall_inf48_4 <- read.delim("./DESeq_analysis/raw_data/counts/Sf-inf-TP02-R04.count", header=F)
colnames(overall_inf48_4) <- c("gene_id", "counts")

#Summarize read counts for viral genes
virus_inf24_1 <- overall_inf24_1 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE")
virus_inf24_2 <- overall_inf24_2 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE") 
virus_inf24_3 <- overall_inf24_3 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE")
virus_inf24_4 <- overall_inf24_4 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE")
virus_inf48_1 <- overall_inf48_1 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE")
virus_inf48_2 <- overall_inf48_2 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE")
virus_inf48_3 <- overall_inf48_3 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE")
virus_inf48_4 <- overall_inf48_4 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% filter (Matches =="TRUE")

virus_summary <- data.frame(virus_inf24_1 = virus_inf24_1$counts,
                                 virus_inf24_2 = virus_inf24_2$counts,
                                 virus_inf24_3 = virus_inf24_3$counts,
                                 virus_inf24_4 = virus_inf24_4$counts,
                                 virus_inf48_1 = virus_inf48_1$counts,
                                 virus_inf48_2 = virus_inf48_2$counts,
                                 virus_inf48_3 = virus_inf48_3$counts,
                                 virus_inf48_4 = virus_inf48_4$counts)

total_virus <- colSums(virus_summary, na.rm=TRUE)

virus_summary_24 <- data.frame(virus_24inf_1 = virus_inf24_1$counts,
                               virus_24inf_2 = virus_inf24_2$counts,
                               virus_24inf_3 = virus_inf24_3$counts,
                               virus_24inf_4 = virus_inf24_4$counts)

virus_summary_48 <- data.frame(virus_48inf_1 = virus_inf48_1$counts,
                               virus_48inf_2 = virus_inf48_2$counts,
                               virus_48inf_3 = virus_inf48_3$counts,
                               virus_48inf_4 = virus_inf48_4$counts)

#Define transgenes 
transgenes <- list("eGFP", "rep78", "cap", "TN", "TN7R", "TN7L", "ITR_R", "ITR_L", "CMV")

#Get number of transgenes counts per sample and summarize
transgenes_inf24_1 <- overall_inf24_1 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")
transgenes_inf24_2 <- overall_inf24_2 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")
transgenes_inf24_3 <- overall_inf24_3 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")
transgenes_inf24_4 <- overall_inf24_4 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")
transgenes_inf48_1 <- overall_inf48_1 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")
transgenes_inf48_2 <- overall_inf48_2 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")
transgenes_inf48_3 <- overall_inf48_3 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")
transgenes_inf48_4 <- overall_inf48_4 %>%  mutate(Matches=grepl("ACNV", gene_id)) %>% mutate(Matches2=grepl("LOC", gene_id)) %>% mutate(Matches3=grepl("AOB78", gene_id)) %>% mutate(Matches4=grepl("Trna", gene_id)) %>% mutate(Matches5=grepl("__", gene_id)) %>% filter (Matches !="TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE" & Matches4 != "TRUE" & Matches5 != "TRUE")

transgenes_summary <- data.frame(transgenes_inf24_1 = transgenes_inf24_1$counts,
                                 transgenes_inf24_2 = transgenes_inf24_2$counts,
                                 transgenes_inf24_3 = transgenes_inf24_3$counts,
                                 transgenes_inf24_4 = transgenes_inf24_4$counts,
                                 transgenes_inf48_1 = transgenes_inf48_1$counts,
                                 transgenes_inf48_2 = transgenes_inf48_2$counts,
                                 transgenes_inf48_3 = transgenes_inf48_3$counts,
                                 transgenes_inf48_4 = transgenes_inf48_4$counts)

total_transgenes <- colSums(transgenes_summary, na.rm=TRUE)

#Get number of insect cell counts per sample and summarize
insect_inf24_1 <- overall_inf24_1 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove24.1 <- which(insect_inf24_1$gene_id %in%  transgenes)
insect_inf24_1 <- insect_inf24_1[-genesToRemove24.1,]

insect_inf24_2 <- overall_inf24_2 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove24.2 <- which(insect_inf24_2$gene_id %in%  transgenes)
insect_inf24_2 <- insect_inf24_2[-genesToRemove24.2,]

insect_inf24_3 <- overall_inf24_3 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove24.3 <- which(insect_inf24_3$gene_id %in%  transgenes)
insect_inf24_3 <- insect_inf24_3[-genesToRemove24.3,]

insect_inf24_4 <- overall_inf24_4 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove24.4 <- which(insect_inf24_4$gene_id %in%  transgenes)
insect_inf24_4 <- insect_inf24_4[-genesToRemove24.4,]

insect_inf48_1 <- overall_inf48_1 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove48.1 <- which(insect_inf48_1$gene_id %in%  transgenes)
insect_inf48_1 <- insect_inf48_1[-genesToRemove48.1,]

insect_inf48_2 <- overall_inf48_2 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove48.2 <- which(insect_inf48_2$gene_id %in%  transgenes)
insect_inf48_2 <- insect_inf48_2[-genesToRemove48.2,]

insect_inf48_3 <- overall_inf48_3 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove48.3 <- which(insect_inf48_3$gene_id %in%  transgenes)
insect_inf48_3 <- insect_inf48_3[-genesToRemove48.3,]

insect_inf48_4 <- overall_inf48_4 %>%  mutate(Matches=grepl("__", gene_id)) %>%  mutate(Matches2=grepl("ACNV", gene_id)) %>%  mutate(Matches3=grepl("ACNV", gene_id)) %>% filter (Matches != "TRUE" & Matches2 != "TRUE" & Matches3 != "TRUE")
genesToRemove48.4 <- which(insect_inf48_4$gene_id %in%  transgenes)
insect_inf48_4 <- insect_inf48_4[-genesToRemove48.4,]

insect_summary <- data.frame(insect_inf24_1 = insect_inf24_1$counts,
                             insect_inf24_2 = insect_inf24_2$counts,
                             insect_inf24_3 = insect_inf24_3$counts,
                             insect_inf24_4 = insect_inf24_4$counts,
                             insect_inf48_1 = insect_inf48_1$counts,
                             insect_inf48_2 = insect_inf48_2$counts,
                             insect_inf48_3 = insect_inf48_3$counts,
                             insect_inf48_4 = insect_inf48_4$counts)

# Define identifier, timepoints and replicates
identifier <- rep(c("inf24_1", "inf24_2", "inf24_3", "inf24_4", "inf48_1", "inf48_2", "inf48_3", "inf48_4"),3)
timepoint_virus <- rep(c(rep("24 hpi",4), rep("48 hpi",4)),3)
replicate_virus <- rep(c("R01", "R02", "R03", "R04"),3)

overall_sf9 <- colSums(insect_summary)
virus_overall<- colSums(virus_summary)
total_transgenes <- colSums(transgenes_summary)
#Summarize data for plotting
summary <- cbind(overall_sf9, virus_overall, total_transgenes)
summary_df <- as.data.frame(summary)
summary_combined <- gather(summary_df)

data_plot <- data.frame(summary_combined, identifier, timepoint_virus, replicate_virus)

#### Plot Figure 2
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

read_distribution_plot_final$condition <- factor(read_distribution_plot_final$condition, levels = c("Pre-infection", "Non-infected", "Infected"))

read_distribution_barplot <- read_distribution_plot_final %>% 
  mutate(option = fct_relevel(option, "trimming_percent", "rrna_percent", "unmapped_percent", "mapped_percent")) %>%
  ggplot(aes(x=sampleNames, y=value, fill=option))+
  geom_bar(position = "fill", stat="identity", width=0.9)+
  theme_script() +
  ylab("% Sequenced reads")+
  facet_grid( ~ condition + timepoint, scales="free_x", space="free_x", switch="x")+
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        strip.background = element_blank())+
  scale_fill_manual(labels = c("trimming_percent"="Quality trimmed", "rrna_percent" ="rRNA", "unmapped_percent"="Unmapped", "mapped_percent"="Uniquely mapped"),values=c("red","#009ACD","#FFA500", "#2C932C"))+ #29BF12
  scale_x_discrete(labels=rep(c("1","2","3","4"),5))+
  scale_y_continuous(labels=scales::percent_format(accuracy=1))+
  theme(strip.placement = "outside",
        axis.text = element_text(size = rel(1.6), colour = "black"),
        axis.title.y = element_text(size = rel(1.8)),
        axis.title.x = element_blank(),
        legend.text = element_text(size = rel(1.6)),
        strip.text.x = element_text(size = rel(1.8)),
        axis.ticks.x = element_blank())

read_distribution_barplot

ggsave(filename= "./bioinformatic_statistics/figures/read_distribution.tiff", plot= read_distribution_barplot, device="tiff", width=8, height=5, dpi = 600)


mapping_distribution <- data_plot %>% 
  mutate(key = fct_relevel(key, "total_transgenes", "virus_overall", "overall_sf9")) %>%
  ggplot(aes(x=identifier, y=value, fill=key))+
  geom_bar(position = "fill", stat="identity", width=0.9)+
  theme_script() +
  ylab("% Mapped reads")+
  facet_grid(. ~ timepoint_virus, scales="free_x", space="free_x", switch="x")+
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        strip.background = element_blank())+
  scale_fill_manual(labels = c("total_transgenes"="AAV transgenes", "overall_sf9"="Host cell genes", "virus_overall"="Baculovirus genes"),values=c("#015501","#147714", "#4CAE4C"))+
  scale_x_discrete(labels=rep(c("1","2","3","4"),2))+
  scale_y_continuous(labels=scales::percent_format(accuracy=1))+
  theme(strip.placement = "outside",
        axis.text = element_text(size = rel(1.6), colour = "black"),
        axis.title.y = element_text(size = rel(1.8)),
        axis.title.x = element_blank(),
        legend.text = element_text(size = rel(1.6)),
        strip.text.x = element_text(size = rel(1.8)),
        axis.ticks.x = element_blank())+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))

mapping_distribution
ggsave(filename= "./bioinformatic_statistics/figures/mapping_distribution.tiff", plot= mapping_distribution, device="tiff", width=5, height=5, dpi = 600)


combined <- ggarrange("", 
                       read_distribution_barplot, 
                       "", 
                       mapping_distribution, ncol = 4, nrow = 1, widths= c(0.01,1.8,0.3,1.2)) 

combined
ggsave(filename= "./bioinformatic_statistics/figures/combined.tiff", plot= combined, device="tiff", width=14, height=6, dpi = 600)
