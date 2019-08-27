setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution")

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final")
date = Sys.Date()

#----------------------------------------------------------------------
#load data
#----------------------------------------------------------------------

#bedtools output files
files = list.files(pattern = ".segs.txt")

#all titan output all patients 
all_titan = fread("all_titan_results_pats.txt")







