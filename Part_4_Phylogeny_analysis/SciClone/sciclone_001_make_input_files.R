#----------------------------------------------------------------------
#sciclone_001_make_input_files.R
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

#define necessary parameters for sciclone
#> table(read_only$Corrected_Call)
#   AMP   GAIN   HETD  HLAMP   HOMD   NEUT
# 71424 113035  23816  30672    451 453867

#SciClone is meant to look at copy neutral regions but we do have many
#amplifications

#vaf data frame:
#1. chromosome
#2. position
#3. reference-supporting read counts
#4. variant-supporting read counts
#5. variant allele fraction (between 0-100)

#copyNumber data frame:
#1. chromosome
#2. segment start position
#3. segment stop position
#4. copy number value for that segment.
#Unrepresented regions are assumed to have a copy number of 2.

#sampleNames: vector of names describing each sample ex: ("Primary
          #Tumor", "Relapse")

#1. READ-ONLY FILE = FINAL MUTATIONS = KEEP THIS WAY UNLESS MAJOR CHANGE NEEDED
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

#2. Copy Number Segment files as used for input to palimpsest and dlbcl classifier
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/TITAN_CNA/results/titan/hmm/optimalClusterSolution_files/titanCNA_ploidy2")






make_input_pyclone = function(input_muts, type){
  muts = fread(input_muts) #get most recent mutation file
  pats = unique(muts$id)
  t = as.data.table(table(muts$mut_id))
  t=t[order(V1)]

  #make sure mutations ordered in the same way in each sample specific file
  #keep only mutations that are copy neutral in all samples
  muts$tot_cn = muts$MajorCN + muts$MinorCN
  muts_keep = as.data.table(filter(muts, tot_cn <= 4))
  muts_sum = filter(as.data.table(table(muts_keep$mut_id)), N ==20)$V1
  muts = as.data.table(filter(muts, mut_id %in% muts_sum))

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
    print(dim(pat_dat))
    write.table(pat_dat, file=paste(pat_dat$id[1], type, "pyclone_input.tsv", sep="_"), quote=F, row.names=F, sep="\t")
  }

  llply(pats, make_input)

}

#make input for all muts
make_input_pyclone("2020-07-30_full_mutations_PYCLONE_INPUT_MUTS.txt", "all_muts")
#make input for some muts
make_input_pyclone("2020-07-30_subset_mutations_PYCLONE_INPUT_MUTS.txt", "subset_muts")
