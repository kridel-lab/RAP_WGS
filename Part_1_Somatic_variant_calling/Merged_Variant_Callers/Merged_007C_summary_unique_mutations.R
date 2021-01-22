#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(tidyverse)
library(tidygraph)
library(ggraph)
library(tidyverse)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

z=which(read_only$Sample == filter(samps,Tissue_Site =="Kidney, NOS 2")$Indiv)
read_only$Tissue_Site[z] = "Kidney, NOS 2"

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

#1. get unique mutations for each patient for each sample 

get_plot = function(patient){

  dat=as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
  unique=as.data.table(table(dat$mut_id, dat$Tissue_Site)) %>% filter(N ==1)
  unique = as.data.table(table(unique$V1)) %>% filter(N==1)
  print(dim(unique))
  unique = as.data.table(filter(dat, mut_id %in% unique$V1))
  print(dim(unique))
  
  #num of unique mutations per sample
  u_per_s=as.data.table(table(unique$Tissue_Site))
  u_per_s =u_per_s[order(-N)]
  colnames(u_per_s) = c("Sample", "N_unique_muts")
  u_per_s$Sample = factor(u_per_s$Sample, levels=unique(u_per_s$Sample))
  
  #make plot
  pdf(paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/", patient, "_summary_unique_mutations_per_sample.pdf", sep=""))
  p<-ggplot(data=u_per_s, aes(x=Sample, y=N_unique_muts)) +
    geom_bar(stat="identity")+theme_minimal()+#+ggtitle("Number of mutations per sample")+
  	theme(axis.text.x = element_text(angle = 60, hjust=1)) + xlab("Sample")+
  	ylab("Number of unique mutations")#+
  	#scale_y_continuous(breaks=seq(0, 400000, by = 25000))
  print(p)
  dev.off()
  
  print("done")
  print(u_per_s)
}

patients = unique(read_only$STUDY_PATIENT_ID)
llply(patients, get_plot)
  
  
  