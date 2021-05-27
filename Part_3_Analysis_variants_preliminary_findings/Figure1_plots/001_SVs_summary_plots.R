#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)
library("ggsci")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

dir.create(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))
setwd(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))

svs = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Manta/2021-04-20_RAP_WGS_all_SVs_heatmap_plot.rds")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#average number of SVs per sample/patient
t=as.data.table(table(svs$STUDY_PATIENT_ID, svs$Indiv))
t=filter(t, N >0)
t %>% group_by(V1) %>% dplyr::summarize(mean = mean(N))
t %>% group_by(V1) %>% dplyr::summarize(sd = sd(N))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

samples_per_mut = as.data.table(table(svs$id, svs$STUDY_PATIENT_ID))
mut_gene = unique(svs[,c("gene", "id", "SVTYPE")])

#summarize number of mutations per sample
z = which((svs$Tissue_Site == "Adrenal gland, NOS") & (svs$STUDY_PATIENT_ID == "LY_RAP_0003"))
svs$Tissue_Site[z] = "Adrenal gland"
z = which((svs$Tissue_Site == "Aorta, ascending, not specified \n\n") & (svs$STUDY_PATIENT_ID == "LY_RAP_0001"))
svs$Tissue_Site[z] = "Aorta, ascending"

muts_per_sample = as.data.table(table(svs$STUDY_PATIENT_ID, svs$Tissue_Site, svs$SVTYPE))
muts_per_sample = as.data.table(filter(muts_per_sample, N >0))
muts_per_sample = muts_per_sample[order(-N)]

colnames(muts_per_sample) = c("Patient", "Sample", "type_SV", "num_of_muts")
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
muts_per_sample$Patient = factor(muts_per_sample$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

#get order of samples for barplot
t = as.data.table(table(svs$STUDY_PATIENT_ID,svs$Tissue_Site))
t = as.data.table(filter(t, N >0))
t = t[order(-N)]

muts_per_sample$Sample = factor(muts_per_sample$Sample, levels=unique(t$V2))
write.table(muts_per_sample, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure1_MAIN_SVs.txt", quote=F, row.names=F, sep="\t")
