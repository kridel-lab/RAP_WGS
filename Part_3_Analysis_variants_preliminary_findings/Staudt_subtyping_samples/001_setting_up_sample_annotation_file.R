#----------------------------------------------------------------------
#001_setting_up_sample_annotation_file.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/LymphGen")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust",
              "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr",
              "ggExtra", "broom", "ggthemes")

lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare sample annotation file to be able to run the tool by Staudt
#https://www.sciencedirect.com/science/article/pii/S1535610820301550

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#sample annotation
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#structural variant annotation
svs = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_RESULTS/PROCESSED_VCFs/2019-12-16_all_SVs_samples.txt")

#filter out BCL2 and BCL6 translocations
svs = as.data.table(filter(svs, gene %in% c("BCL2", "BCL6"), SVTYPE == "BND"))

unique(svs$SVTYPE)
unique(svs$pat)

#format file needs to be in
#column1 - sample name
#column2 - CNA available (1=yes/0=no)
#column3 - BCL2 translocation (1=yes/0=no)
#column4 - BCL6 translocation (1=yes/0=no)

input_f = as.data.table(samps[,c(Indiv)])
input_f$CNA = 1
colnames(input_f) = c("Sample.ID", "Copy.Number")
input_f$BCL2.transloc = 1
input_f$BCL6.transloc = 0

write.table(input_f, paste(date, "LymphGen_Sample_Annotation_File.txt", sep="_"), quote=F, row.names=F, sep="\t")
