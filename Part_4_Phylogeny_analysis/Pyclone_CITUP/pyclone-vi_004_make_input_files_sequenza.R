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

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")
cnas_full = cnas
colnames(cnas_full)[1]="samplename"

cnas = unique(cnas[,c("Sample", "Tissue_Site", "Purity")])
colnames(cnas)[1]="samplename"

#save sample id versus sample name clean
patients= c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")

#all muts with CNA status
all_muts = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/2021-09-25_Mutect2_Strelka_merged_mutations_wCNA_status.rds")

#just keep patient 1 data, mut id, sample, CN major and CN minor
all_muts = unique(all_muts %>% filter(STUDY_PATIENT_ID == "LY_RAP_0001") %>% select(Sample, mut_id, Nmin, Nmaj))
colnames(all_muts)[1] = "samplename"
colnames(all_muts)[2] = "id"

#2. SNV/CNA input data for pyclone

make_input_pyclone = function(input_muts){
  muts = fread(input_muts) #get most recent mutation file
  pats = unique(muts$samplename)
  t = as.data.table(table(muts$mut_id))
  t=t[order(V1)]
  patient = unlist(strsplit(input_muts, "_mutations_PYCLONE_INPUT_MUTS"))[1]
  patient = paste(unlist(strsplit(patient, "_"))[2:4], collapse="_")

  #make sure mutations ordered in the same way in each sample specific file

  if(patient == "LY_RAP_0001"){
  muts_sum = filter(as.data.table(table(muts$mut_id)), N ==3)$V1
  }

  if(patient == "LY_RAP_0002"){
  muts_sum = filter(as.data.table(table(muts$mut_id)), N ==4)$V1}

  if(patient == "LY_RAP_0003"){
  muts_sum = filter(as.data.table(table(muts$mut_id)), N ==20)$V1}

  muts = as.data.table(filter(muts, mut_id %in% muts_sum))
  muts = merge(muts, cnas, by = "samplename")

  #prepare matrix for pyclone vi input
  muts$normal_cn = 2

  muts = muts %>% select(mut_id, samplename, Ref_counts, alt_counts, MajorCN, MinorCN, normal_cn, Purity)
  colnames(muts) = c("mutation_id", "sample_id", "ref_counts", "alt_counts", "major_cn", "minor_cn", "normal_cn", "tumour_content")
  print(table(muts$major_cn))

  mut_rm = unique(muts[which(is.na(muts$major_cn)),]$mutation_id)
  muts = filter(muts, !(mutation_id %in% mut_rm))

  mut_rm = unique(muts[which(muts$major_cn == 0),]$mutation_id)
  muts = filter(muts, !(mutation_id %in% mut_rm))

  print(tail(muts))
  print(dim(muts))
  write.table(muts, file=paste("all_samples_pyclonevi", patient, "pyclone_input.tsv", sep="_"), quote=F, row.names=F, sep="\t")
  print("done!")
}

list.files(pattern="_mutations_PYCLONE_INPUT_MUTS")

#make input for all muts
make_input_pyclone("2021-09-28_LY_RAP_0001_mutations_PYCLONE_INPUT_MUTS.txt")

#make input for all muts
make_input_pyclone("2021-06-03_LY_RAP_0002_mutations_PYCLONE_INPUT_MUTS.txt")

#make input for all muts
make_input_pyclone("2021-06-03_LY_RAP_0003_mutations_PYCLONE_INPUT_MUTS.txt")
