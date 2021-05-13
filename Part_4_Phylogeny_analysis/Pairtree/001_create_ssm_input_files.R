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
library(clonevol)
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
p001_pyclone_input = fread("all_samples_pyclonevi_all_muts_LY_RAP_0001_pyclone_input.tsv")
p002_pyclone_input = fread("all_samples_pyclonevi_all_muts_LY_RAP_0002_pyclone_input.tsv")
p003_pyclone_input = fread("all_samples_pyclonevi_all_muts_LY_RAP_0003_pyclone_input.tsv")

#pyclone output files
p001_pyclone_output = fread("LY_RAP_0001_rap_wgs_all_muts.tsv")
p002_pyclone_output = fread("LY_RAP_0002_rap_wgs_all_muts.tsv")
p003_pyclone_output = fread("LY_RAP_0003_rap_wgs_all_muts.tsv")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#prep everything
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  dat$sample_id = factor(dat$sample_id, levels=unique(dat$sample_id))

  #get all unique mutations
  muts = unique(dat$name)

  concatenate_muts = function(mut){
    print(mut)
    mut_dat = filter(dat, name == mut)
    name=mut
    var_reads = paste(mut_dat$var_reads, collapse=",")
    total_reads = paste(mut_dat$total_reads, collapse=",")
    var_read_prob = paste(mut_dat$var_read_prob, collapse=",")
    samples = paste(mut_dat$sample_id, collapse=",")
    res = c(name, var_reads, total_reads, var_read_prob, samples)
    return(res)
  }

  all_muts = as.data.table(ldply(llply(muts, concatenate_muts)))
  if(length(unique(all_muts$V5)) == 1){
    colnames(all_muts) = c("name", "var_reads", "total_reads", "var_read_prob", "samples")
    all_muts$id = paste("s", 0:(nrow(all_muts)-1), sep="")
    all_muts = all_muts[,c("id", "name", "var_reads", "total_reads", "var_read_prob", "samples")]
    print(head(all_muts))
    return(all_muts)
  }

  print("done!")
}

pair_tree_input_params = function(ssm_file, py_out){

  #samples (required): list of sample names, provided in the same order as
  #entries in the var_reads and total_reads vectors in the .ssm file.
  #Each sample name can be an arbitrary string.

  #clusters (required): list-of-lists denoting variant clusters, each of which
  #corresponds to a subclone. Each sub-list contains strings corresponding to
  #the IDs of all variants that belong to a given cluster. All variants belonging
  #to a cluster are assumed to share the same subclonal frequency in a given
  #tissue sample, with the VAFs of each variant serving as a noisy estimate of
  #this subclonal frequency after using var_read_prob to correct for ploidy.
  #See the section on clustering mutations for instructions on how to use Pairtree
  #to build clusters. Pairtree assumes that all variants within a cluster share
  #the same subclonal frequencies across all cancer samples. Ploidy-corrected
  #estimates of these frequencies can be obtained by dividing a variant's VAF by its variant read probability `omega.

  #this is what it needs to look like:
  #id	name	var_reads	total_reads	var_read_prob
  #s0	S_0	54,175,196	1000,1000,1000	0.5,0.5,0.5
  colnames(py_out)[1] = "name"
  py_out = merge(py_out, ssm_file, by="name")
  samples = unique(py_out$samples)
  samples = unlist(strsplit(samples, ","))

  #make list of lists with which clusters mutations fall into
  muts_clusters = unique(py_out[,c("cluster_id", "id")])
  muts_clusters = muts_clusters[order(cluster_id)]
  muts_list = (split(muts_clusters$id, muts_clusters$cluster_id))

  muts_list_vec=vector(mode="list", length=length(unique(muts_clusters$cluster_id)))
  for(i in 1:length(muts_list)){
    print(i)
    muts_list_vec[[i]]=muts_list[[i]]
  }

  list1 <- vector(mode="list", length=3)
  list1[[1]] <- samples
  list1[[2]] <- muts_list_vec
  list1[[3]] <- ""

  names(list1) = c("samples", "clusters", "garbage")
  return(list1)
}


p1 = pair_tree_input_ssm(p001_pyclone_input, p001_pyclone_output)
p1_samples = unique(p1$samples)
p1_json = pair_tree_input_params(p1, p001_pyclone_output)

p2 = pair_tree_input_ssm(p002_pyclone_input, p002_pyclone_output)
p2_samples = unique(p2$samples)
p2_json = pair_tree_input_params(p2, p002_pyclone_output)

p3 = pair_tree_input_ssm(p003_pyclone_input, p003_pyclone_output)
p3_samples = unique(p3$samples)
p3_json = pair_tree_input_params(p3, p003_pyclone_output)

p1$samples = NULL
p2$samples = NULL
p3$samples = NULL

#save ssm input files in pairtree directory
mainDir="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pairtree"
subDir=paste(date, "input_files", sep="_")

dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

write.table(p1, file="p001_ssm_input.ssm", quote=F, row.names=F, sep="\t")
write_json(p1_json, "p001_input.params.json")
