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
library(ggpubr)

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

mypal = c("#E5DFD9","#EAD286" ,"#D1EB7B", "#96897F" ,"#E5C0A6" ,
  "#72A93B", "#74DAE3" ,"#49B98D" ,"#D97B8F" ,"#70A2A4", "#64709B" ,"#DFBF38" ,"#61EA4F" ,
  "#C7CBE7", "#786DDA",
"#CFA0E0" ,"#67E9D0" ,"#7C9BE1", "#D94753" ,
"#AAE6B0", "#D13BDF" ,"#DEAEC7" ,"#BBE6DF" ,"#B2B47A" ,"#E6ECBA", "#C86ED7",
 "#7BEE95" ,"#6F46E6" ,"#65B9E0", "#C0EC3E",
"#DE8D54" ,"#DF4FA6")

#1. get unique mutations for each patient for each sample

get_plot = function(patient){

  dat=as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))

  #make plot
  pdf(paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/", patient, "_summary_VAFs_of_mutations.pdf", sep=""))

  set.seed(10)
  p=ggdensity(dat, x = "gt_AF",
   color = "Tissue_Site", fill = "Tissue_Site",
   palette = sample(mypal, length(unique(dat$Tissue_Site))))+theme_bw()
  print(p)

  dev.off()

  print("done")
}

patients = unique(read_only$STUDY_PATIENT_ID)
llply(patients, get_plot)
