#----------------------------------------------------------------------
#001_setting_up_sample_annotation_file.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)

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

#Our mutations
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

#results files
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/LymphGen")

#result
res = fread("rap_wgs_Result.txt")

#compare
comp = fread("rap_wgs_Compare.txt")
