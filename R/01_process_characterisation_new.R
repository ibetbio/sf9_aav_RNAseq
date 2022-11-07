#!/usr/bin/env Rscript
#===============================================================================
#' Author: Nikolaus Virgolini
#' Date: 2022/May
#' Transcriptome analysis of Sf9 insect cells during production of recombinant AAV 
#===============================================================================
#This script allows the reproduction of the AAV production characterization shown in Figure 1 of the manuscript.

## Load necessary packages
library("dplyr")
library ("tidyr")
library("ggplot2")
library("DT")
library("ggpubr")
library("forcats")
library("RColorBrewer")
library("ggbreak")
library("scales")

### Set theme for plotting
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

### Create a directory for saving figures
figsdir <- "./process_characterisation/figures/"
if (!dir.exists(figsdir)) {
  dir.create(figsdir)
}

### Create a directory for saving intermediate data
intermediatadir <- "./process_characterisation/intermediate_data/"
if (!dir.exists(intermediatadir)) {
  dir.create(intermediatadir)
}

### Get growth and viability data
data <- read.csv("process_characterisation/raw_data/Raw_data_R.csv",header=TRUE,fill=TRUE, sep=",")

### Define sample names, timepoints and replicates
replicate_names <- c("R01","R02","R03","R04")
timepoints <- c("TP00", "TP01", "TP02", "TP03", "TP04", "TP05", "TP06")
conditions <- c("pre", "non", "inf")

### Manipulate data for easier plotting
data_plot_inf <- data %>% filter(condition != "non") %>% mutate(condition_2 = "inf")
data_plot_non <- data %>% filter(condition != "inf") %>% mutate(condition_2 = "non")
data_plot <- rbind(data_plot_inf, data_plot_non)

################################################################################################
####################################Plot Figure 1b##############################################
################################################################################################
# Define x-, y1-, and y2-axis for plot 1b
axis_text_y1 <- bquote(atop(bold('Viable cell concentration'),
                             '['~10^{6}~'cell.mL'^{-1}~']'))
axis_text_y2 <- bquote(bold('Cell viability') ~ '[%]')
axis_text_x <- bquote(bold('Time after infection') ~ '[hours]')

plotting <- data_plot[-c(29:40),]
plotting$condition <- factor(plotting$condition, levels = c("pre", "non", "inf"))
plotting$condition2 <- factor(plotting$condition2, levels = c("pre", "non", "inf"))

growth.data <- ggplot(plotting, aes(x=time, y=VCD, color=condition)) + 
  geom_point(aes(fill = condition),alpha=0.5, size=6) +
  scale_shape_manual(labels=c("pre" = " Pre-infection","inf" = "Infected", "non" = "Non-infected"), values = c(21,21,21),guide = "none")+
  geom_point(aes(x=time, y=viability /10*1.6, color=condition, shape=condition, fill = "white"),alpha=0.5, size=6) +
  geom_smooth(data_plot, method="lm", formula=y ~ poly(x, 5), se=FALSE, mapping=aes(x=time, y=VCD, color=condition_2),linetype="dotted", size=0.6) +
  geom_smooth(data_plot, method="lm", formula=y ~ poly(x, 5), se=FALSE, mapping=aes(x=time, y=viability /10 * 1.6, color=condition_2),linetype="dotted",size=0.6)+ 
  geom_vline(xintercept=c(0), linetype="dashed", size=0.5, color="grey") +
  scale_y_continuous(name=axis_text_y1, breaks=seq(0,16, by=4), sec.axis=sec_axis(~.*10/1.6, name=axis_text_y2))+
  scale_x_continuous(name = axis_text_x, breaks=seq(-48,100, by=24))+
  scale_color_manual(values=c("#4daf4a", "#377eb8","#e41a1c", "white"), labels=c("pre" = "Pre-infection", "inf" = "Infected", "non" = "Non-infected")) +
  scale_fill_manual(values=c("#4daf4a", "#377eb8","#e41a1c", "white"), guide = "none") +
  theme_script()+
  geom_text(aes(x=0, y=4.1), label="\U2193", colour="black", size=10)+
  geom_text(aes(x=24, y=6.1), label="\U2193", colour="black", size=10)+
  geom_text(aes(x=48, y=10.1), label="\U2193", colour="black", size=10)+
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.text = element_text(size=rel(1.8)),
    axis.title = element_text(size = rel(2.2)),
    axis.text=element_text(size=rel(1.8))
  ) +
  guides(color = guide_legend(override.aes=list(shape = c(21,21,21), fill = c("#4daf4a", "#377eb8","#e41a1c"), color = c("#4daf4a", "#377eb8","#e41a1c"))))

growth.data 

### Save figure 1b in figures directory
ggsave(filename= "/figure_1b.tiff", plot= growth.data, path = figsdir, device="tiff", width=21.3, height=12.7, unit = "cm", dpi = 600)


################################################################################################
################################# Metabolic flux determination #################################
################################################################################################

### Calculate IVCC for metabolic analysis
ivcc_inf <- data.frame()
for(i in seq(from=1,to=length(unique(data_plot$replicate))))
{  replicate_ivcc_inf <- data_plot %>%
  group_by(timepoint) %>%
  filter(replicate==replicate_names[i] & condition !="non" & timepoint !="TP-2" & timepoint !="TP-1" & condition_2=="inf")
ivcc_inf[1,i] <- replicate_ivcc_inf$VCD[1]
for(j in seq(from=1,to=length(unique(data_plot$timepoint))-3))
{ 
  ivcc_inf[j+1,i] <- ivcc_inf[j,i]+mean(replicate_ivcc_inf$VCD[j:c(j+1)])*(replicate_ivcc_inf$time[j+1]-replicate_ivcc_inf$time[j]) 
}
}

ivcc_non <- data.frame()
for(i in seq(from=1,to=length(unique(data_plot$replicate))))
{  replicate_ivcc_non <- data_plot %>%
  group_by(timepoint) %>%
  filter(replicate==replicate_names[i] & condition !="inf" & timepoint !="TP-2" & timepoint !="TP-1" & condition_2=="non")
ivcc_non[1,i] <- replicate_ivcc_non$VCD[1]
for(j in seq(from=1,to=length(unique(data_plot$timepoint))-3))
{ 
  ivcc_non[j+1,i] <- ivcc_non[j,i]+mean(replicate_ivcc_non$VCD[j:c(j+1)])*(replicate_ivcc_non$time[j+1]-replicate_ivcc_non$time[j]) 
}
}

### Change the column and rownames and show data
colnames(ivcc_non) <- c("non_R01", "non_R02", "non_R03", "non_R04")
colnames(ivcc_inf) <- c("inf_R01", "inf_R02", "inf_R03", "inf_R04")
rownames(ivcc_non) <- timepoints[1:5]
rownames(ivcc_inf) <- timepoints[1:5]

print(ivcc_non)
print(ivcc_inf)

### Combine Infected and Non-Infected datasets
ivcc <- data.frame(ivcc_non, ivcc_inf)

### Make a table to read out data
write.table(ivcc, file="process_characterisation/intermediate_data/ivcc_data.txt", append= FALSE, sep="\t")

### Read raw metabolic data into R
metabolites <- read.csv("process_characterisation/Raw_data/metabolites_R.csv",sep=",")

### Define metabolite names
metabolites_names <- c("Glc", "Gln", "Glu", "Lac", "NH3")

### Gather data, bind rows and make a data.frame
ivcc_combined <- gather(bind_cols(ivcc[-1,]), value="ivcc")
metab_ordered <- bind_rows(metabolites %>% arrange(replicate) %>% filter (condition == "non") %>% dplyr::select(Glc, Gln, Glu, Lac, NH3, replicate, condition), metabolites %>% arrange(replicate) %>% filter (condition == "inf") %>% dplyr::select(Glc, Gln, Glu, Lac, NH3, replicate, condition))
metabolites_values <- data.frame(ivcc_combined, metab_ordered)

### Define formula to calculate metabolic rates
Formula <- formula(paste(paste(metabolites_names[1]),"~ivcc", collapse="+"))
gl_t<- lm(Formula, data=metabolites_values %>% group_by(replicate) %>% filter(condition=="inf", replicate==replicate_names[1]))
gl_t$coefficients["ivcc"]

### Calculate metabolic rates
metabrates_non <- data.frame()
metabrates_error_non <- data.frame()
for(i in seq(from=1,to=length(unique(metabolites$replicate))))
{  
  for(j in seq(from=1,to=length(metabolites_names)))
  { 
    metabsub_non <- metabolites_values %>% filter(condition=="non", replicate==replicate_names[i]) %>% dplyr::select(2,c(2+j))
    Formula <- formula(paste(paste(metabolites_names[j]),"~ivcc", collapse="+"))
    metabrate_results_non <- lm(Formula, metabsub_non)
    metabrates_non[i,j] <- coef(summary(metabrate_results_non))["ivcc","Estimate"]
    metabrates_error_non[i,j] <- coef(summary(metabrate_results_non))["ivcc","Std. Error"]
  }
}

metabrates_inf <- data.frame()
metabrates_error_inf <- data.frame()
for(i in seq(from=1,to=length(unique(metabolites$replicate))))
{  
  for(j in seq(from=1,to=length(metabolites_names)))
  { 
    metabsub_inf <- metabolites_values %>% filter(condition=="inf", replicate==replicate_names[i]) %>% dplyr::select(2,c(2+j))
    Formula <- formula(paste(paste(metabolites_names[j]),"~ivcc", collapse="+"))
    metabrate_results_inf <- lm(Formula, metabsub_inf)
    metabrates_inf[i,j] <- coef(summary(metabrate_results_inf))["ivcc","Estimate"]
    metabrates_error_inf[i,j] <- coef(summary(metabrate_results_inf))["ivcc","Std. Error"]
  }
}

### Change row and column names
rownames_metabrates <- c("non","non","non","non","inf","inf","inf","inf")
colnames(metabrates_non) <- metabolites_names
rownames(metabrates_non) <- c("non_R01","non_R02","non_R03","non_R04")
colnames(metabrates_inf) <- metabolites_names
rownames(metabrates_inf) <- c("inf_R01","inf_R02","inf_R03","inf_R04")

### Combined data sets and multiply by 1000
metabrates <- rbind(metabrates_non, metabrates_inf)*1000
metabolic_rates <- cbind(format(round(metabrates, 2), nsmall = 2), condition=rownames_metabrates, replicate=c(replicate_names, replicate_names))

### Make an html table for metabolic rates
table2 <- metabolic_rates %>%
  DT::datatable(
    rownames = T,
    options = list(
      pagelength = T,
      scrollX = T,
      columnDefs = list(list(
        className = "dt-head-center dt-center",
        targets = 1:6
      ))
    )
  )

html2 <- "metabolic_rates_table.html"
saveWidget(table2, paste(intermediatadir, html2, sep = "/"))

### Manipulate data for plotting
plot_data <- data.frame(gather(metabolic_rates, "Glc", "Gln", "Glu", "Lac", "NH3", key=metabolite, value=metabolic_rate))

Glc <- plot_data %>% group_by(condition) %>% filter(metabolite== "Glc")
Gln <- plot_data %>% group_by(condition) %>% filter(metabolite== "Gln")
Glu <- plot_data %>% group_by(condition) %>% filter(metabolite== "Glu")
Lac <- plot_data %>% group_by(condition) %>% filter(metabolite== "Lac")
NH3 <- plot_data %>% group_by(condition) %>% filter(metabolite== "NH3")

### Perform student t-test for significance testing
Glc_test <- t.test(x= as.numeric(Glc$metabolic_rate[1:4]), y= as.numeric(Glc$metabolic_rate[5:8]), alternative="two.sided")
Gln_test <-t.test(x= as.numeric(Gln$metabolic_rate[1:4]), y= as.numeric(Gln$metabolic_rate[5:8]), alternative="two.sided")
Glu_test <-t.test(x= as.numeric(Glu$metabolic_rate[1:4]), y= as.numeric(Glu$metabolic_rate[5:8]), alternative="two.sided")
Lac_test <-t.test(x= as.numeric(Lac$metabolic_rate[1:4]), y= as.numeric(Lac$metabolic_rate[5:8]), alternative="two.sided")
NH3_test <-t.test(x= as.numeric(NH3$metabolic_rate[1:4]), y= as.numeric(NH3$metabolic_rate[5:8]), alternative="two.sided")

### Combine data for plotting
rates_sum <- plot_data %>% group_by(interaction(metabolite, condition)) %>% summarise(n=n(), mean=mean(as.numeric(metabolic_rate)), sd=sd(as.numeric(metabolic_rate)))
rates_sum_plot <- data.frame(rates_sum, condition=c(rep("inf",5),rep("non",5)), metabolites=c(rep(metabolites_names,2)))

position <- rates_sum_plot %>% dplyr::filter(condition=="inf") %>% mutate(position = case_when(
  mean < 0 ~ mean-sd-1,
  mean > 0 ~ mean+sd+1
))

p_values <- data.frame(values=c(Glc_test$p.value, Gln_test$p.value, Glu_test$p.value, Lac_test$p.value, NH3_test$p.value), metabolites=metabolites_names, position=position$position[])

### Set astherics for plotting
p_values <- p_values %>% mutate(signif = case_when(
  values > 0.05  ~ "ns",
  values > 0.01 & values <= 0.05  ~ "*",
  values > 0.001 & values <= 0.01  ~ "**",
  values <= 0.001  ~ "***"
))

### Make barplot for figure 1c
axis_text_y <- bquote(atop(bold('Metabolic rates'),
                           '[nmol.10'^{6}~'cell'^{-1}~'.h'^{-1}~']'))

barplot_rates <- ggplot(rates_sum_plot %>% mutate(condition = fct_relevel(condition, 
                                                                          "non","inf")) , aes(x=condition, y=mean, fill=metabolites))+
  geom_bar(stat="identity", position="dodge", colour="black",width=0.75)+
  geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), position=position_dodge(0.75), width=0.25)+
  theme_script() +
  geom_hline(aes(yintercept = 0))+
  theme(
    axis.text = element_text(size = rel(1.8)),
    axis.title = element_text(face = "bold", size = rel(2.2)),
    legend.text = element_text(size = rel(1.8)),
    legend.title = element_blank(),
    legend.position = "bottom") +
  labs(y =axis_text_y, x = "Condition") +
  scale_fill_brewer(palette = "Greys", name="Metabolites", labels=c("Glc"="Glucose", "Gln"="Glutamine", "Glu"="Glutamate", "Lac"="Lactate","NH3"="Ammonia"))+
  scale_x_discrete(labels=c("inf"="Infected", "non"="Non-infected"))

###Labels with brackets
p_values_new <- data.frame(values=c(Glc_test$p.value, Gln_test$p.value, Glu_test$p.value, Lac_test$p.value, NH3_test$p.value), metabolites=metabolites_names, position=position$position[], group1=c(0.7,0.85,1,1.15,1.3), group2=c(1.7,1.85,2,2.15,2.3), y.position=c(5, 8, 11, 14, 17))
p_values_new <- p_values_new %>% mutate(signif = case_when(
  values > 0.05  ~ "ns",
  values > 0.01 & values <= 0.05  ~ "*",
  values > 0.001 & values <= 0.01  ~ "**",
  values <= 0.001  ~ "***"
))
barplot_final_3 <- barplot_rates + stat_pvalue_manual(p_values_new, label = "signif", size = 5.5) 
barplot_final_3

### Save figure 1c
ggsave(filename= "/barplot_metbrates_plot.tiff", plot= barplot_final_3, path = figsdir, device="tiff", width=21.3, height=12.7, unit = "cm", dpi = 600)

################################################################################################
################### Plot full and empty AAV capsids for intracellular samples ##################
################################################################################################
### Load data
total_data <- read.csv("process_characterisation/Raw_data/total_capsids_R.csv",sep=",")

total_sum <- total_data %>% group_by(interaction(time, option, method)) %>% summarise(n=n(), mean=mean(capsids), sd=sd(capsids), time=unique(time), condition=unique(option), method=unique(method))

labels <- c(elisa = "Total particles", qpcr = "Genomic particles")

# Plot total AAV particles
y_axis_total_particles <- bquote(atop(bold('Intracellular total AAV'),
                                        bold("particles")~"[TP.mL" ^{-1} ~"]"))

elisa_plot <- ggplot(total_sum %>% filter(condition=="intra" & time != "0h" & method =="elisa"),aes(x=time, y=mean, fill=condition))+
  geom_bar(stat="identity", position=position_dodge(), colour="black", size=.3, fill="grey")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, position=position_dodge(.9))+
  theme_script() +
  xlab(bquote(bold('Time after infection') ~ '[hours]'))+
  scale_y_break(c(10000, 10000000))+
  scale_y_log10(name=y_axis_total_particles, limits = c(1, 100000000000), labels = trans_format("log10", label_math()))+
  scale_x_discrete(labels=c("24h" = "24", "48h" = "48", "72h" = "72","96h" = "96"))+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=rel(2)),
    axis.title = element_text(size = rel(2.2)),
    axis.text=element_text(size=rel(2)),
    legend.position = "none",
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank())
  
  

elisa_plot
ggsave(filename= "/total_capsids.tiff", plot= elisa_plot, path = figsdir, device="tiff", width=14, height=12.7, unit = "cm", dpi = 600)


# Plot genomic AAV particles
y_axis_genomic_particles <- bquote(atop(bold('Intracellular genomic AAV'),
                                         bold("particles")~"[VG.mL" ^{-1} ~"]"))

qpcr_plot <- ggplot(total_sum %>% filter(condition=="intra" & time != "0h" & method =="qpcr"),aes(x=time, y=mean, fill=condition))+
  geom_bar(stat="identity", position=position_dodge(), colour="black", size=.3, fill="grey")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2, position=position_dodge(.9))+
  theme_script() +
  xlab(bquote(bold(~"Time after infection") ~ "[hours]"))+
  scale_y_break(c(10000, 10000000))+
  scale_y_log10(name=y_axis_genomic_particles, limits = c(1, 100000000000), labels = trans_format("log10", label_math()))+
  scale_x_discrete(labels=c("24h" = "24", "48h" = "48", "72h" = "72","96h" = "96"))+
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size=rel(2)),
    axis.title = element_text(size = rel(2.2)),
    axis.text=element_text(size=rel(2)),
    legend.position = "none",
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank())+
  geom_text(aes(x = 1, y = 2), label = "n.d.", size = 6)

qpcr_plot

ggsave(filename= "/genomic_capsids.tiff", plot= qpcr_plot, path = figsdir, device="tiff", width=14, height=12.7, unit = "cm", dpi = 600)


# Plot total extracellular baculovirus particles
baculovirus <- read.csv("process_characterisation/Raw_data/baculovirus.csv",sep=",")


y_axis_bacv_particles <- bquote(atop(bold('Extracellular baculovirus particles'),
                                     bold('concentration')~"["~italic("ie1")~"genome copies.mL" ^{-1} ~"]"))

x_axis_bacv_particles <- bquote(bold('Cell viability') ~ '[%]')

bacv_plot <- ggplot(baculovirus,aes(x=viability, y=capsids, fill = time, color = time, size = 15))+
  geom_point(shape = 21, size = 10)+
  theme_script() +
  xlab(bquote(bold(~"Time after infection") ~ "[hours]"))+
  scale_y_log10(name=y_axis_bacv_particles, limits=c(1000000, 100000000000), labels = trans_format("log10", label_math()))+
  scale_x_continuous(name=x_axis_bacv_particles, limits = c(0, 100))+
  scale_fill_manual(values = c("grey10", "grey30", "grey50", "grey70"), name = "Time after\ninfection [h]", labels=c("24h" = "24", "48h" = "48", "72h" = "72","96h" = "96"))+
  scale_color_manual(values = c("grey10", "grey30", "grey50", "grey70"), name = "Time after\ninfection [h]", labels=c("24h" = "24", "48h" = "48", "72h" = "72","96h" = "96"))+
  theme(
    legend.title = element_text(size=rel(2)),
    legend.text = element_text(size=rel(2)),
    axis.title = element_text(size = rel(2.2)),
    axis.text=element_text(size=rel(2)),
    legend.position = "right",
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank())

bacv_plot

ggsave(filename= "/figS3.tiff", plot= bacv_plot, path = figsdir, device="tiff", width=22, height=18, unit = "cm", dpi = 600)

#####End of script

