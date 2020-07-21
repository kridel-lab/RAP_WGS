#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

#purpose: complete final soft filters
#mainly intersect with copy number file
#generate filters required for Treeomics, PhyloWGS and Pyclone
#so that VCFs can be filtered accordingly

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
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom")
lapply(packages, require, character.only = TRUE)
library(RColorBrewer)
library(openxlsx)
library(plotly)
library(readxl)
library(GenomicRanges)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare input mutation files for treeomics
#-just protein coding gene mutations
#-all mutations and then run with -wes parameter

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. READ-ONLY FILE = FINAL MUTATIONS = KEEP THIS WAY UNLESS MAJOR CHANGE NEEDED
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

#----------------------------------------------------------------------
#SEPERATE FILES INTO FINAL SETS BASED ON TOOL INPUT
#----------------------------------------------------------------------

#+++++ PCG mutations only +++++++++++++++++++++++++++++++++++++++++++++

#2. TREEOMICS INPUT VCF FILES = remove mutations in noncoding genes
#keep all copy number states
treeomics_input = as.data.table(filter(read_only,
!(ExonicFunc.ensGene %in% c("unknown", "."))))

#BED file for intersecting VCF files with this file
#CHROM POS tab seperated
#create one file for each patient sample - sample specific mutations to extract

get_pat_treeomics_muts = function(pat, dat, type){
  print(dim(dat))
  pat_dat = as.data.table(filter(dat, Indiv == pat))
  pat_input = unique(pat_dat[,c("CHROM", "POS")])
  pat_input$CHROM = sapply(pat_input$CHROM, function(x){unlist(strsplit(x, "chr"))[2]})
  file_name = paste(pat, "treeomics_input", type, ".bed", sep="_")
  write.table(pat_input, file=paste(
    "/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Treeomics/input_dat/",
  file_name, sep=""), quote=F, row.names=F, sep="\t")
  print("done")
}

pats = unique(treeomics_input$Indiv)
llply(pats, get_pat_treeomics_muts, treeomics_input, "pcgs_only")

#+++++ ALL mutations +++++++++++++++++++++++++++++++++++++++++++++++++

llply(pats, get_pat_treeomics_muts, read_only, "all_muts")
