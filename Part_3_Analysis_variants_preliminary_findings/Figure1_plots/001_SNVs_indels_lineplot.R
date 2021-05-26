#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)
library("ggsci")
require(gridExtra)
library(cowplot)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#data table for barplot
barplot = as.data.table(table(samples_per_mut$num_of_samples_with_mut, samples_per_mut$patient))
barplot = as.data.table(filter(barplot, N >0))
colnames(barplot) = c("num_of_samples_with_mut", "patient", "num_of_muts")
barplot$num_of_samples_with_mut = factor(barplot$num_of_samples_with_mut, levels=unique(barplot$num_of_samples_with_mut))
barplot$patient[barplot$patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
barplot$patient[barplot$patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
barplot$patient[barplot$patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
barplot$patient = factor(barplot$patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
colnames(barplot)[2]="Patient"

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/001_samples_per_mutation_lineplot.pdf",
width=5, height=5)
# Basic barplot
p<-ggline(barplot, x="num_of_samples_with_mut", y="num_of_muts",
palette = c("#00AFBB", "#E7B800", "#FC4E07"), color="Patient")+
xlab("# of samples sharing mutation") + ylab("Mutation count")+
scale_y_continuous(breaks=seq(0,300000,20000))
print(p)
dev.off()
