#----------------------------------------------------------------------
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
              "plyr", "openxlsx",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "ggthemes")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare input files for pyclone-vi
#Input format
#mutation_id - Unique identifier for the mutation. This is free form but should match across all samples.
#sample_id - Unique identifier for the sample.
#ref_counts - Number of reads matching the reference allele.
#alt_counts - Number of reads matching the alternate allele.
#major_cn - Major copy number of segment overlapping mutation.
#minor_cn - Minor copy number of segment overlapping mutation.
#normal_cn - Total copy number of segment in healthy tissue. For autosome this will be two and male sex chromosomes one.
#tumour_content - The tumour content (cellularity) of the sample. Default value is 1.0 if column is not present

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. purity and sample data
purity=fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Palimpsest/input/annotation_data_palimpsest_input.txt")
colnames(purity)[1] = "Indiv"
purity$Gender = NULL

samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")
samps = samps %>% select(id, Indiv)
samps = merge(samps, purity, by = "Indiv")

#2. SNV/CNA input data for pyclone

make_input_pyclone = function(input_muts, type){
  muts = fread(input_muts) #get most recent mutation file
  pats = unique(muts$id)
  t = as.data.table(table(muts$mut_id))
  t=t[order(V1)]

  #make sure mutations ordered in the same way in each sample specific file
  #try analysis assuming mutations are copy neutral
  muts$MajorCN=1
  muts$MinorCN=1
  muts$tot_cn = muts$MajorCN + muts$MinorCN
  muts_keep = as.data.table(filter(muts, tot_cn <= 4))
  muts_sum = filter(as.data.table(table(muts$mut_id)), N ==20)$V1
  muts = as.data.table(filter(muts, mut_id %in% muts_sum))
  muts = merge(muts, samps, by = "id")

  #prepare matrix for pyclone vi input
  muts$normal_cn = 2
  muts = as.data.table(filter(muts, mut_id %in% muts_sum))

  muts = muts %>% select(mut_id, id, Ref_counts, alt_counts, MajorCN, MinorCN, normal_cn, Purity)
  colnames(muts) = c("mutation_id", "sample_id", "ref_counts", "alt_counts", "major_cn", "minor_cn", "normal_cn", "tumour_content")
  print(table(muts$major_cn))

  print(tail(muts))
  print(dim(muts))
  write.table(muts, file=paste("all_samples_pyclonevi", type, "neutral_pyclone_input.tsv", sep="_"), quote=F, row.names=F, sep="\t")

}

#make input for all muts
make_input_pyclone("2020-09-13_full_mutations_PYCLONE_INPUT_MUTS.txt", "all_muts")
#make input for some muts
make_input_pyclone("2020-09-13_subset_mutations_PYCLONE_INPUT_MUTS.txt", "subset_muts")
