#-------------------------------------------------------------------------------
#001_SNVs_indels_summary.R
#This scripts preps mutations summary for Figure 1 barplot
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

#All SNVs--------------------------------------------------------------

#summarize number of mutations per sample
z = which((read_only$Tissue_Site == "Adrenal gland, NOS") & (read_only$STUDY_PATIENT_ID == "LY_RAP_0003"))
read_only$Tissue_Site[z] = "Adrenal gland"
z = which((read_only$Tissue_Site == "Aorta, ascending, not specified \n\n") & (read_only$STUDY_PATIENT_ID == "LY_RAP_0001"))
read_only$Tissue_Site[z] = "Aorta, ascending"

muts_per_sample = as.data.table(table(read_only$STUDY_PATIENT_ID,read_only$Tissue_Site))
muts_per_sample = as.data.table(filter(muts_per_sample, N >0))
muts_per_sample = muts_per_sample[order(-N)]
colnames(muts_per_sample) = c("Patient", "Sample", "num_of_muts")
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
muts_per_sample$Patient = factor(muts_per_sample$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
muts_per_sample$Sample = factor(muts_per_sample$Sample, levels=unique(muts_per_sample$Sample))
muts_per_sample$N_trans = log1p(muts_per_sample$num_of_muts)

sample_order = unique(muts_per_sample$Sample)
write.table(muts_per_sample, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure1_MAIN_SNVs_Indels_ALL.txt", quote=F, row.names=F, sep="\t")
write.table(sample_order, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure1_MAIN_sample_order.txt", quote=F, row.names=F, sep="\t")

#Coding only--------------------------------------------------------------------

#get only SNVs/indels in coding regions
coding_only = filter(read_only, Func.ensGene %in% c("exonic", "splicing", "exonic\\x3bsplicing"),
ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "stopgain", "frameshift_deletion", "frameshift_insertion", "stoploss"))

muts_per_sample = as.data.table(table(coding_only$STUDY_PATIENT_ID,coding_only$Tissue_Site, coding_only$ExonicFunc.ensGene))
muts_per_sample = as.data.table(filter(muts_per_sample, N >0))
muts_per_sample = muts_per_sample[order(-N)]
colnames(muts_per_sample) = c("Patient", "Sample", "mut_type", "num_of_muts")
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
muts_per_sample$Patient = factor(muts_per_sample$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
muts_per_sample$Sample = factor(muts_per_sample$Sample, levels=sample_order)
muts_per_sample$N_trans = log1p(muts_per_sample$num_of_muts)

write.table(muts_per_sample, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure1_MAIN_SNVs_Indels_CODING.txt", quote=F, row.names=F, sep="\t")
