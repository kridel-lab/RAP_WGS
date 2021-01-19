#----------------------------------------------------------------------
#sciclone_001_make_input_files.R
#----------------------------------------------------------------------

date = Sys.Date()

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools", "sciClone")
lapply(packages, require, character.only = TRUE)

setwd("/Users/kisaev/Documents/Hatchet_analysis")
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
snv_list = readRDS("list_snvs_df.rds")

#2. sample names (use ids for now)
samples = names(snv_list)

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#run sciclone, one sample at a time
sc = sciClone(vafs=snv_list, minimumDepth=50,
         sampleNames=samples, maximumClusters=10, verbose=TRUE) #, clusterMethod="binomial.bmm")

#save output
#create output
writeClusterTable(sc, "sciclone_output.txt")
sc.plot2d(sc, "clusters_results.pdf")
