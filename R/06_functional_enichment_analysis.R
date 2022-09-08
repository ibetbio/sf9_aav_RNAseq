#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Spodoptera frugiperda Sf9 cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of the DESeq2 analysis needed for all scripts from 4 onwards.

### Load necessary packages
library("clusterProfiler")
library("DESeq2")
library("apeglm")
library("tidyverse")
library("xlsx")

intermdir <- c("./functional_enrichment_analysis/intermediate_data")
result_dir <- c("./functional_enrichment_analysis/result_dir")
load("./DESeq_analysis/intermediate_results/deres_infvsnon24.rData")
load("./DESeq_analysis/intermediate_results/deres_inf48vs24.rData")

#Read Drosophila melanogaster GO-slim reference
go_slim <- read.delim("functional_enrichment_analysis/raw_data/go_slim_drosophila.tsv", header=F)

######This code is only needed for the first time, when the trinotate output file is shortened to Drosophila GO slim reference###########
#########################################################################################################################################
#reference <- read.delim("functional_enrichment_analysis/raw_data/gos.txt", header=F)
#gos_keep <- which(reference$V1 %in% go_slim$V1)
#reference_new <- reference[gos_keep,]
#write.csv(reference_new, file="./functional_enrichment_analysis/raw_data/gos_slim.csv", quote=F)
#########################################################################################################################################
# Manipulate reference file for ClusterProfiler input
reference <- read.csv("./functional_enrichment_analysis/raw_data/gos_slim.csv")
reference <- reference[,-1]
reference_new <- as.data.frame(t(reference))
colnames(reference_new) <- reference_new[1,]
reference_new <- reference_new[-1,]
reference_final <- gather(reference_new, key="gene", value="geneSet")

### Overrepresentation analysis of DGE's of early infection stages
enrichment_24 <- clusterProfiler::enricher(gene=row.names(deres_infvsnon24),
                                                    pvalueCutoff=0.05,
                                                    pAdjustMethod="BH",
                                                    minGSSize=10,
                                                    maxGSSize=500,
                                                    qvalueCutoff = 0.2,
                                                    TERM2GENE = reference_final)

# Summarise result to create a results table
results_24infvsnon <- as.data.frame(enrichment_24@result)
results_24infvsnon_sign <- results_24infvsnon %>% filter(p.adjust<=0.05)
ratio_24 <- strsplit(results_24infvsnon_sign$GeneRatio, split="/")
ratio_24 <- unlist(ratio_24)
calc_ratio_24 <- c()
x <- 1
for(i in seq(from=1,to=c(length(ratio_24)-1), by=2)){
  calc_ratio_24[x] <- as.numeric(ratio_24[i])/as.numeric(ratio_24[i+1])
  x <- x+1}
enrichment_24infvsnon_results <- data.frame(results_24infvsnon_sign, ratio = calc_ratio_24, day = c(rep("24hpi", length(results_24infvsnon_sign$geneID))))
enrichment_24infvsnon_annotated <- dplyr::left_join(enrichment_24infvsnon_results, go_slim, by=c("ID" = "V1"))

fn_24inf <- paste(result_dir, "/tableS6.xlsx",sep="")
if (file.exists(fn_24inf)) {file.remove(fn_24inf)}

# save the results to Excel
write.xlsx(enrichment_24infvsnon_annotated,
           file = fn_24inf,
           sheetName = "ORA - Early Infection",
           row.names=F)


### Overrepresentation analysis of DGE's of late infection stages
enrichment_48vs24 <- clusterProfiler::enricher(gene=row.names(deres_inf48vs24),
                                                     pvalueCutoff=0.05,
                                                     pAdjustMethod="BH",
                                                     minGSSize=10,
                                                     maxGSSize=500,
                                                     qvalueCutoff = 0.2,
                                                     TERM2GENE = reference_final)

# Summarise result to create a results table
results_48vs24 <- as.data.frame(enrichment_48vs24@result)
results_48vs24_sign <- results_48vs24 %>% filter(p.adjust<=0.05)
ratio_48vs24<- strsplit(results_48vs24_sign$GeneRatio, split="/")
ratio_48vs24 <- unlist(ratio_48vs24)
calc_ratio_48vs24 <- c()
x <- 1
for(i in seq(from=1,to=c(length(ratio_48vs24)-1), by=2)){
  calc_ratio_48vs24[x] <- as.numeric(ratio_48vs24[i])/as.numeric(ratio_48vs24[i+1])
  x <- x+1}
enrichment_48vs24_results <- data.frame(results_48vs24_sign, ratio = calc_ratio_48vs24, day = c(rep("48hpi", length(results_48vs24_sign$geneID))))
enrichment_48vs24_annotated <- dplyr::left_join(enrichment_48vs24_results, go_slim, by=c("ID" = "V1"))

fn_48inf <- paste(result_dir, "/tableS7.xlsx",sep="")
if (file.exists(fn_48inf)) {file.remove(fn_48inf)}

# save the results to Excel
write.xlsx(enrichment_48vs24_annotated,
           file = fn_48inf,
           sheetName = "ORA - Late Infection",
           row.names=F)

#Combine data for later plotting
combined_enrichment_results <- rbind(enrichment_24infvsnon_annotated, enrichment_48vs24_annotated)

##############################################################################
##########################Save data for other scripts#########################
##############################################################################
save(combined_enrichment_results, file="./functional_enrichment_analysis/intermediate_data/combined_enrichment_results.rData")


protein_folding <- as.data.frame(reference_new[,"GO:0006457"])
colnames(protein_folding) <- "genes"
protein_folding <- protein_folding %>% filter(genes != "na")
protein_folding <- c(protein_folding$genes)
save(protein_folding, file="./functional_enrichment_analysis/intermediate_data/protein_folding.rData")

cell_division <- as.data.frame(reference_new[,"GO:0051301"])
colnames(cell_division) <- "genes"
cell_division <- cell_division %>% filter(genes != "na")
cell_division <- c(cell_division$genes)
save(cell_division, file="./functional_enrichment_analysis/intermediate_data/cell_division.rData")

mitotic_cc <- as.data.frame(reference_new[,"GO:0000278"])
colnames(mitotic_cc) <- "genes"
mitotic_cc <- mitotic_cc %>% filter(genes != "na")
mitotic_cc <- c(mitotic_cc$genes)
save(mitotic_cc, file="./functional_enrichment_analysis/intermediate_data/mitotic_cc.rData")

virus_response_cc <- as.data.frame(reference_new[,"GO:0051607"])
colnames(virus_response_cc) <- "genes"
virus_response_cc <- virus_response_cc %>% filter(genes != "na")
virus_response_cc <- c(virus_response_cc$genes)
save(virus_response_cc, file="./functional_enrichment_analysis/intermediate_data/virus_response_cc.rData")

#### End of script

