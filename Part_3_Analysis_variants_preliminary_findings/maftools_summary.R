#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)
library(ccube)
#> Loading required package: foreach
library(dplyr)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final")
date = Sys.Date()

source("/cluster/home/kisaev/scripts/annovar_to_maftools.R")

#----------------------------------------------------------------------
#load data
#----------------------------------------------------------------------

#bedtools output files
files = list.files(pattern = "bedtools_output.bed")

#all titan output all patients 
all_titan = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/all_titan_results_pats.txt")
z = which(all_titan$Sample == "Sample")
all_titan = all_titan[-z,]

#purity 
purity = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution.txt")

#cna vcf sample conversion file
conversion = fread("cna_vcf_sample_conversion.csv")

#annovar txt files
anno_files = list.files(pattern="multianno.txt")

#----------------------------------------------------------------------
#get maf files
#----------------------------------------------------------------------

get_maf = function(file){
	maf = annovarToMaf(file, table="ensGene")
	write.table(maf, file=paste(file, "maf_file.txt", sep="_"), quote=F, row.names=F, sep="\t")
	print("dome")
	return(maf)
}

mafs = llply(anno_files, get_maf) 

#in bash 
#cat *maf* > all_mafs_all_samples.txt


