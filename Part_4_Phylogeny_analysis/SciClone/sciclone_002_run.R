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

#here try running using full mutation matrix produced via bamreadcount
#which should include read counts for mutations not called in cases
#try running sciclone on all samples and all mutations to see if it works
#otherwise might need to reduce number of mutations

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

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#run sciclone, one sample at a time
sc = sciClone(vafs=snv_list, useSexChrs=FALSE,
         copyNumberCalls=cnas_list,
         sampleNames=samples, minimumDepth=50, annotation=annotations)

#------------------------------------
#2d clustering using two samples:
sc = sciClone(vafs=snv_list[4:5],
    copyNumberCalls=cnas_list[4:5],
    sampleNames=samples[c(4,5)], minimumDepth=50)

#create output
writeClusterTable(sc, "clusters2")
sc.plot1d(sc,"clusters2.1d.pdf")
sc.plot2d(sc,"clusters2.2d.pdf")


#save output
#create output
writeClusterTable(sc, paste(samples[index], "sciclone_output.txt", sep="_"))
sc.plot1d(sc, paste(samples[index], "clusters_results.pdf", sep="_"))
