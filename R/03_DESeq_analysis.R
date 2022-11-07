#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Spodoptera frugiperda Sf9 cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of the DESeq2 analysis needed for all scripts from 4 onwards.

# Load necessary packages
library("DESeq2")
library("dplyr")
library("tidyverse")
library("xlsx")

# Define input and output directory
count_dir <- "./DESeq_analysis/raw_data/counts/"
result_dir <- "./DESeq_analysis/deseq_results"
intermdir <- "./DESeq_analysis/intermediate_results"
remove <- read.delim("DESeq_analysis/reference_files/genes_to_remove_deseq.txt", header =F)

# Define File names
count_file_names <- grep("count", list.files(count_dir), value=T)

##########################################################################################################
######### Run DESeq2 for all samples for sample-distance, PCA analysis and heatmaps ######################
##########################################################################################################
# Define conditions
replicate <-   c("R01","R02","R03","R04",
                 "R01","R02","R03","R04",
                 "R01","R02","R03","R04",
                 "R01","R02","R03","R04",
                 "R01","R02","R03","R04")

condition <- c("inf_24hpi","inf_24hpi","inf_24hpi","inf_24hpi",
               "inf_48hpi","inf_48hpi","inf_48hpi","inf_48hpi",
               "non_24hpi","non_24hpi","non_24hpi","non_24hpi",
               "non_48hpi","non_48hpi","non_48hpi","non_48hpi",
               "seed_0hpi","seed_0hpi","seed_0hpi","seed_0hpi")

sample_information <- data.frame(sampleName = count_file_names,
                                 fileName =  count_file_names,
                                 condition = condition,
                                 replicate = replicate)

DESeq_data_overall <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_information,
  directory = count_dir,
  design = ~ replicate + condition)

# Remove viral and transgene counts
genesToRemove <- which(!rownames(DESeq_data_overall) %in%  remove$V1)
DESeq_data_overall_IC <- DESeq_data_overall[genesToRemove,]

# Run DESeq
DESeq_data_overall_IC <- DESeq(DESeq_data_overall_IC)

## Save DESeq_datasets for all genes and insect only genes for further analysis in other scripts
save(DESeq_data_overall_IC, file = paste(intermdir, "/DESeq_data_overall_IC.rData", sep=""))

##########################################################################################################
######################### Run DESeq analysis of infected vs non-infected samples at 24 hpi###################
##########################################################################################################
# Define conditions
cond_24infvsnon <-   c("inf", "inf", "inf", "inf",
                       "non", "non", "non", "non")


si_24infvsnon <- data.frame(sampleName = count_file_names[c(1:4, 9:12)],
                            fileName =  count_file_names[c(1:4, 9:12)],
                            condition = cond_24infvsnon)

DESeq_data_24infvsnon <- DESeqDataSetFromHTSeqCount(
  sampleTable = si_24infvsnon,
  directory = count_dir,
  design = ~ condition)

# Set 24 hpi Non-Infected as reference
colData(DESeq_data_24infvsnon)$condition <- factor(colData(DESeq_data_24infvsnon)$condition,
                                                   levels = c("inf", "non"))

DESeq_data_24infvsnon$condition <- relevel(DESeq_data_24infvsnon$condition, "non")

# Remove viral and transgene counts
genesToRemove <- which(!rownames(DESeq_data_24infvsnon) %in%  remove$V1)
DESeq_data_24infvsnon <- DESeq_data_24infvsnon[genesToRemove,]

# Run DESeq
DESeq_data_24infvsnon <- DESeq(DESeq_data_24infvsnon)

# Extract differentially expressed genes
res_infvsnon24 <- results(DESeq_data_24infvsnon, independentFiltering = T)

deres_infvsnon24 <- subset(res_infvsnon24,
                            abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >= 50)
deres_infvsnon24 <- deres_infvsnon24[order(deres_infvsnon24$log2FoldChange, decreasing = T),]

#Save a .csv file with results for all genes
fn <- paste(result_dir, "/res_infvsnon24.csv",sep="")
if (file.exists(fn)) {file.remove(fn)}
write.csv(res_infvsnon24,
          file = fn)

#Save data for downstream analysis in other scripts
save(DESeq_data_24infvsnon, file = paste(intermdir, "/DESeq_data_24infvsnon.rData", sep=""))
save(res_infvsnon24, file = paste(intermdir, "/res_infvsnon24.rData", sep=""))
save(deres_infvsnon24, file = paste(intermdir, "/deres_infvsnon24.rData", sep=""))

#Annotate significant differentially expressed gene list with RefSeq names and Trinotate output
annotation <- read.csv("./DESeq_analysis/reference_files/gene_names.csv", header = T)

trinotate_annotation <- read.table("./DESeq_analysis/reference_files/annot_feature_map.txt", header=F, sep="\t")

deres_infvsnon24_df <- as.data.frame(deres_infvsnon24)
deres_infvsnon24_annotated <- dplyr::left_join(rownames_to_column(deres_infvsnon24_df), annotation, by = c("rowname" = "symbol"))
deres_infvsnon24_final <- dplyr::left_join(deres_infvsnon24_annotated, trinotate_annotation, by = c("rowname" = "V1"))

# Save the results as an Excel file
fn_24inf <- paste(result_dir, "/tableS3.xlsx",sep="")
if (file.exists(fn_24inf)) {file.remove(fn_24inf)}

write.xlsx(deres_infvsnon24_final,
          file = fn_24inf,
          sheetName = "DE genes Sf9 - 24 hpi",
          row.names=F)

####################################################################################################
#######################Run DESeq analysis of infected samples 48 hpi vs 24 hpi######################
####################################################################################################
# Define conditions
cond_inf48vs24 <- c("24", "24", "24", "24",
                    "48", "48", "48", "48")


si_inf48vs24 <- data.frame(sampleName = count_file_names[c(1:8)],
                           fileName =  count_file_names[c(1:8)],
                           condition = cond_inf48vs24)

DESeq_data_inf48vs24 <- DESeqDataSetFromHTSeqCount(
  sampleTable = si_inf48vs24,
  directory = count_dir,
  design = ~ condition)

# Set 24 hpi Infected as reference
colData(DESeq_data_inf48vs24)$condition <- factor(colData(DESeq_data_inf48vs24)$condition,
                                                  levels = c("48", "24"))

DESeq_data_inf48vs24$condition <- relevel(DESeq_data_inf48vs24$condition, "24")

# Remove viral and transgene counts
genesToRemove <- which(!rownames(DESeq_data_inf48vs24) %in%  remove$V1)
DESeq_data_inf48vs24 <- DESeq_data_inf48vs24[genesToRemove,]

# Run DESeq
DESeq_data_inf48vs24 <- DESeq(DESeq_data_inf48vs24)

# Extract differentialy expressed genes
res_inf48vs24 <- results(DESeq_data_inf48vs24, independentFiltering = T)

deres_inf48vs24 <- subset(res_inf48vs24,
                           abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >=50)
deres_inf48vs24 <- deres_inf48vs24[order(deres_inf48vs24$log2FoldChange, decreasing = T),]

#Save a .csv file with results for all genes
fn <- paste(result_dir, "/res_inf48vs24.csv",sep="")
if (file.exists(fn)) {file.remove(fn)}
write.csv(res_inf48vs24,
          file = fn)

# Save data for other analysis and plotting
save(DESeq_data_inf48vs24, file = paste(intermdir, "/DESeq_data_inf48vs24.rData", sep=""))
save(res_inf48vs24, file = paste(intermdir, "/res_inf48vs24.rData", sep=""))
save(deres_inf48vs24, file = paste(intermdir, "/deres_inf48vs24.rData", sep=""))

#Annotate significant differentially expressed gene list with RefSeq names and Trinotate output
deres_inf48vs24_df <- as.data.frame(deres_inf48vs24)
deres_inf48vs24_annotated <- dplyr::left_join(rownames_to_column(deres_inf48vs24_df), annotation, by = c("rowname" = "symbol"))
deres_inf48vs24_final <- dplyr::left_join(deres_inf48vs24_annotated, trinotate_annotation, by = c("rowname" = "V1"))

# Save the results to Excel
fn_48inf <- paste(result_dir, "/tableS4.xlsx",sep="")
if (file.exists(fn_48inf)) {file.remove(fn_48inf)}

write.xlsx(deres_inf48vs24_final,
           file = fn_48inf,
           sheetName = "DE genes Sf9 - 48 hpi",
           row.names=F)

###################################################################################################
#######################Run DESeq analysis of infected samples 48 hpi vs 24 hpi######################
#######################################for bacv genes###############################################
####################################################################################################
#Define conditions
cond_bacv <- c("24", "24", "24", "24",
               "48", "48", "48", "48")


si_bacv <- data.frame(sampleName = count_file_names[c(1:8)],
                      fileName =  count_file_names[c(1:8)],
                      condition = cond_bacv)

DESeq_data_bacv <- DESeqDataSetFromHTSeqCount(
  sampleTable = si_bacv,
  directory = count_dir,
  design = ~ condition)

# Set 24 hpi - Infected as reference
colData(DESeq_data_bacv)$condition <- factor(colData(DESeq_data_bacv)$condition,
                                             levels = c("48", "24"))

DESeq_data_bacv$condition <- relevel(DESeq_data_bacv$condition, "24")

# Remove insect genes
genesToKeep <- which(rownames(DESeq_data_bacv) %in%  remove$V1)
DESeq_data_bacv <- DESeq_data_bacv[genesToKeep,]

# Run DESeq
DESeq_data_bacv <- DESeq(DESeq_data_bacv)


# Extract differentially expressed genes
res_bacv <- results(DESeq_data_bacv, independentFiltering = T)

deres_bacv <- subset(res_bacv,
                     abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >=50)
deres_bacv <- deres_bacv[order(deres_bacv$log2FoldChange, decreasing = T),]

# Save overall results as .csv
fn <- paste(result_dir, "/res_bacv.csv",sep="")
if (file.exists(fn)) {file.remove(fn)}
write.csv(res_bacv,
          file = fn)

# Save intermediate data for further analysis
save(DESeq_data_bacv, file = paste(intermdir, "/DESeq_data_bacv.rData", sep=""))
save(res_bacv, file = paste(intermdir, "/res_bacv.rData", sep=""))
save(deres_bacv, file = paste(intermdir, "/deres_bacv.rData", sep=""))

#Annotate significant differential expressed gene list
bacv_annotation <- read.table("./DESeq_analysis/reference_files/GCF_000838485.1_ViralProj14023_feature_table.txt", header=T, sep="\t")

deres_bacv_df <- as.data.frame(deres_bacv)
deres_bacv_final <- dplyr::left_join(rownames_to_column(deres_bacv_df), bacv_annotation, by = c("rowname" = "locus_tag"))

# Save the results to Excel
fn_48inf <- paste(result_dir, "/tableS5.xlsx",sep="")
if (file.exists(fn_48inf)) {file.remove(fn_48inf)}

write.xlsx(deres_bacv_final,
           file = fn_48inf,
           sheetName = "DE genes bacv and transgenes - 48 hpi",
           row.names=F)

####################################################################################################
#######################Extract transgene sequences important for rAAV###############################
####################################################################################################
# Extract normalized counts for all transgenes
eGFP <- data.frame(plotCounts(DESeq_data_bacv, gene=c("eGFP"), intgroup="condition", returnData = T), gene = rep("eGFP",8))
rep_gene <- data.frame(plotCounts(DESeq_data_bacv, gene=c("rep78"), intgroup="condition", returnData = T), gene = rep("rep",8))
cap_gene <- data.frame(plotCounts(DESeq_data_bacv, gene=c("cap"), intgroup="condition", returnData = T), gene = rep("cap",8))

# Combine and save data for plotting
transgene_norm_count <- rbind(eGFP, rep_gene, cap_gene)
save(transgene_norm_count, file = paste(intermdir, "/transgenes_norm_count.rData", sep=""))

####################################################################################################
#######################DESeq analysis for Non-Infected samples over time######################
####################################################################################################

cond_non <- c("24", "24", "24", "24",
              "48", "48", "48", "48",
              "0", "0", "0", "0")


si_non <- data.frame(sampleName = count_file_names[c(9:20)],
                     fileName =  count_file_names[c(9:20)],
                     condition = cond_non)

DESeq_data_non <- DESeqDataSetFromHTSeqCount(
  sampleTable = si_non,
  directory = count_dir,
  design = ~ condition)

colData(DESeq_data_non)$condition <- factor(colData(DESeq_data_non)$condition,
                                            levels = c("48", "24", "0"))

DESeq_data_non$condition <- relevel(DESeq_data_non$condition, "0")

genesToRemove <- which(!rownames(DESeq_data_non) %in%  remove$V1)
DESeq_data_non <- DESeq_data_non[genesToRemove,]

DESeq_data_non <- DESeq(DESeq_data_non)

res_non24 <- results(DESeq_data_non, independentFiltering = T, lfcThreshold = 0, contrast = c("condition", "24", "0"))
res_non48 <- results(DESeq_data_non, independentFiltering = T, lfcThreshold = 0, contrast = c("condition", "48", "24"))


deres_non24 <- subset(res_non24,
                       abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >= 50)
deres_non24 <- deres_non24[order(deres_non24$log2FoldChange, decreasing = T),]

deres_non48 <- subset(res_non48,
                       abs(log2FoldChange) >= 0.5849625 & padj < 0.05 & baseMean >= 50)
deres_non48 <- deres_non48[order(deres_non48$log2FoldChange, decreasing = T),]


save(DESeq_data_non, file = paste(result_dir, "/DESeq_data_non.rData", sep=""))
save(res_non24, file = paste(result_dir, "/res_non24.rData", sep=""))
save(res_non48, file = paste(result_dir, "/res_non48.rData", sep=""))
save(deres_non24, file = paste(result_dir, "/deres_non24.rData", sep=""))
save(deres_non48, file = paste(result_dir, "/deres_non48.rData", sep=""))

### End of Script