#use mutations evaluted by pyclone-VI

#keep track of ordering of samples by commas

#var_reads: a comma-separated vector of how many whole-number-valued
#genomic reads at this mutation's locus corresponded to the variant allele
#in each tissue sample. Tissue samples may be provided in any order, so long as
#this order is consistent for all mutations. The names of the associated tissue
#samples are given in the .params.json file, detailed below.

#var_read_prob: if a copy-number aberration (CNA) duplicated the reference allele
#in the lineage bearing mutation j prior to j occurring, there will be two
#reference alleles and a single variant allele in all cells bearing j, such that \omega_{js} = 0.3333

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#load everything
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
#library(clonevol)
library("gplots")
library(threadr)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr", "openxlsx",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "ggthemes")
lapply(packages, require, character.only = TRUE)

#set up patients
patients= c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")

#pyclone input files
p001_pyclone_input = fread("all_samples_pyclonevi_LY_RAP_0001_pyclone_input.tsv")
p002_pyclone_input = fread("all_samples_pyclonevi_LY_RAP_0002_pyclone_input.tsv")
p003_pyclone_input = fread("all_samples_pyclonevi_LY_RAP_0003_pyclone_input.tsv")

#pyclone output files
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone/10-06-2021")
p001_pyclone_output = fread("all_samples_pyclonevi_LY_RAP_0001_beta-binomial_rap_wgs_all_muts.tsv")
p002_pyclone_output = fread("all_samples_pyclonevi_LY_RAP_0002_beta-binomial_rap_wgs_all_muts.tsv")
p003_pyclone_output = fread("all_samples_pyclonevi_LY_RAP_0003_beta-binomial_rap_wgs_all_muts.tsv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#prep everything
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

mut_info = unique(read_only[,c("mut_id", "REF", "ALT", "symbol", "STUDY_PATIENT_ID", "Tissue_Site", "Sample",
"Func.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene", "cosmic68", "avsnp142", "biotype")])

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pairtree")

pair_tree_input_ssm = function(py_in, py_out){

  #this is what it needs to look like:
  #id	name	var_reads	total_reads	var_read_prob
  #s0	S_0	54,175,196	1000,1000,1000	0.5,0.5,0.5
  dat = py_in
  dat = as.data.table(filter(dat, mutation_id %in% py_out$mutation_id))

  dat$name = dat$mutation_id
  dat$var_reads = dat$alt_counts
  dat$total_reads = dat$ref_counts + dat$alt_counts
  dat$var_read_prob = dat$major_cn/ (dat$major_cn + dat$minor_cn)
  dat = dat[,c("name", "var_reads", "total_reads", "var_read_prob", "sample_id")]

  #keep only mutations in clusters with at least 20 mutations
  clusts_muts = unique(py_out[,c("cluster_id", "mutation_id")])
  t = as.data.table(table(clusts_muts$cluster_id))
  t = filter(t, N >30)
  muts_keep = filter(py_out, cluster_id %in% t$V1)

  t$new_order = 1:(nrow(t))
  colnames(t)= c("cluster_id", "num_muts_in_cluster", "pairtree_cluster_name")
  t$cluster_id = as.numeric(t$cluster_id)

  patient = as.character(dat$sample_id[1])
  patient = paste(unlist(strsplit(patient,"_"))[1:3], collapse="_")

  #merge with mutation info
  mut_info_pat = filter(mut_info, STUDY_PATIENT_ID == patient)
  colnames(muts_keep)[which(colnames(muts_keep) == "sample_id")] = "Sample"
  colnames(muts_keep)[which(colnames(muts_keep) == "mutation_id")] = "mut_id"

  muts_keep = merge(muts_keep, mut_info_pat, by=c("Sample", "mut_id"))
  muts_keep = merge(muts_keep, t, by="cluster_id")

  #find genes that are convergent (mutated across multiple clones)
  just_coding = filter(muts_keep, biotype == "protein_coding", !(AAChange.ensGene == "."), !(ExonicFunc.ensGene == "synonymous_SNV"))
  tt = as.data.table(table(just_coding$symbol, just_coding$pairtree_cluster_name)) %>% filter(N > 0)
  colnames(tt) = c("gene", "cluster", "N")
  conv_genes = as.data.table(table(tt$gene)) %>% filter(N > 1)
  tt= filter(tt, gene %in% conv_genes$V1)
  just_coding = filter(just_coding, symbol %in% tt$gene)
  muts_keep$convergent = ""
  z = which(muts_keep$mut_id %in% just_coding$mut_id)
  muts_keep$convergent[z] = "convergent"

  #write results
  library(openxlsx)
  write.xlsx(muts_keep, paste(patient, date, "all_muts", 'Pyclone-VI-pairtree-results.xlsx', sep="_"))

  print("done!")

}

pair_tree_input_ssm(p001_pyclone_input, p001_pyclone_output)
pair_tree_input_ssm(p002_pyclone_input, p002_pyclone_output)
pair_tree_input_ssm(p003_pyclone_input, p003_pyclone_output)
