#----------------------------------------------------------------------
#sciclone_001_make_input_files.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone")

args = commandArgs(trailingOnly = TRUE)
index = args[1]
index = as.numeric(index)
print(index)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom")
lapply(packages, require, character.only = TRUE)
library(RColorBrewer)
library(openxlsx)
library(plotly)
library(readxl)
library(GenomicRanges)
library(sciClone)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare input mutation files for sciclone
#-just protein coding gene mutations
#-all mutations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. SNV list
snv_list = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_snvs_df.rds")

#2. Copy Number List
cnas_list = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_cnas_df.rds")

#3. sample names (use ids for now)
samples = names(snv_list)

#4. annotation file
annotations = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/annotation_file.rds")

#5. loh regions
loh = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_loh_df.rds")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#run sciclone, one sample at a time
sc = sciClone(vafs=snv_list[[index]], useSexChrs=TRUE, regionsToExclude=loh[[index]] ,
         copyNumberCalls=cnas_list[[index]],
         sampleNames=samples[index], minimumDepth=50, annotation=annotations)

#save output
#create output
writeClusterTable(sc, paste(samples[index], "sciclone_output.txt", sep="_"))
sc.plot1d(sc, paste(samples[index], "clusters_results.pdf", sep="_"))
