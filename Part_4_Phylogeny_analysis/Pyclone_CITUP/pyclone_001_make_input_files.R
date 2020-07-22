#----------------------------------------------------------------------
#make_pyclone_input.R
#karin isaev
#last updated: August 30th,2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "ggthemes")
lapply(packages, require, character.only = TRUE)
library(RColorBrewer)
library(openxlsx)
library(plotly)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare input files for pyclone

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. SNV/CNA input data for pyclone

muts = fread(list.files(pattern="PYCLONE_INPUT_MUTS.txt")[length(list.files(pattern="PYCLONE_INPUT_MUTS.txt"))]) #get most recent mutation file
pats = unique(muts$id)
t = as.data.table(table(muts$mut_id))
t=t[order(V1)]

#make sure mutations ordered in the same way in each sample specific file

make_input = function(patient){
  pat_dat = as.data.table(filter(muts, id == patient))
  pat_dat$normal_cn = 2
  pat_dat = unique(pat_dat[,c("mut_id", "Ref_counts", "alt_counts", "normal_cn",
  "MinorCN", "MajorCN", "symbol", "Func.ensGene", "id")])
  print(table(pat_dat$MajorCN))

  z = which(pat_dat$mut_id %in% t$V1)
  pat_dat = pat_dat[z,]

  colnames(pat_dat) = c("mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn", "gene_name", "region", "id")
  print(tail(pat_dat))
  write.table(pat_dat, file=paste(pat_dat$id[1], "pyclone_input.tsv", sep="_"), quote=F, row.names=F, sep="\t")
}

llply(pats, make_input)
