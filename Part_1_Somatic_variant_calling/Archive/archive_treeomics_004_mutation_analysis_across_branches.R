#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Treeomics/Treeomics_WES_wpurities")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust",
              "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr",
              "ggExtra", "broom", "ggthemes", "readxl")

lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#analyze results from treeomics
#which DLBCL mutations or other known cancer gene mutations
#potentially give rise to new region specific tumours

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#sample annotation
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#treeomics mutation location results
tree=fread("STRELKA_VCF_INPUTS_20_mf=0_1_e=0_01_c0=0_5_af=0_05_mps=5_variants.csv")

#change names of patient sample id to regions of the body
tree$samps_w_mut = ""
get_samp_names = function(x){
  pats = unlist(strsplit(x, "__"))
  z = which(samps$Indiv %in% pats)
  pats_new = samps$id[z]
  return(pats_new)
}
tree$samps_w_mut = sapply(tree$BranchName, get_samp_names)

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#which mutations lead to which organ
#start with driver mutations as annotated by Treeomics
filter(tree, Driver == TRUE)[c(1:10, 93)]

#see if other DLBCL genes are here from Morin supplementary
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))

#filter(tree, GeneSymbol %in% morin$Gene)

genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 2))

#which are candidates for convergent evolution
ce = as.data.table(table(tree$GeneSymbol, tree$Phylogeny))
ce = as.data.table(filter(ce, N >0))
ce = as.data.table(table(ce$V1))
ce = as.data.table(filter(ce, N >1))
