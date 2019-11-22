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
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

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

#summarized snvs and cnvs from 21 sequencing folders 
#here, summarize number of mutations/sample/location
#which genes are mutated across all sites which are unique?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. SNV/CNA input data for pyclone 

muts = fread(list.files(pattern="PYCLONE_INPUT_MUTS.txt")[length(list.files(pattern="PYCLONE_INPUT_MUTS.txt"))]) #get most recent mutation file
pats = unique(muts$Indiv)

make_input = function(patient){
  pat_dat = as.data.table(filter(muts, Indiv == patient))  
  pat_dat$normal_cn = 2
  pat_dat = unique(pat_dat[,c("mut_id", "Ref_counts", "alt_counts", "normal_cn", "Nmin", "Nmaj", "hg19.ensemblToGeneName.value", "Func.ensGene", "id")])
  print(table(pat_dat$Nmaj))
  pat_dat = as.data.table(filter(pat_dat, Nmaj > 0))
  colnames(pat_dat) = c("mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn", "gene_name", "region", "id")
  write.table(pat_dat, file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone/", patient, "pyclone_input.tsv", sep="_"), quote=F, row.names=F, sep="\t")
}

llply(pats, make_input)



