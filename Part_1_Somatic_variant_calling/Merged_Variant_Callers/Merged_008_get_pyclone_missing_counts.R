#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

date = Sys.Date()
print(date)
options(scipen=999) #no scientific notation

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
source("/cluster/home/kisaev/scripts/bam_readcount_parseline.R")

#load libraries
packages <- c("dplyr", "ggplot2", "tidyr", "data.table", "plyr",
	"stringr")
lapply(packages, require, character.only = TRUE)
library("readxl")
library(GenomicRanges)
library(params)
library(readr)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#used bam readcount to extract counts from all samples
#overlapping mutations that were only found in some samples and not all
#Pyclone input requires these values

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#Our mutations
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone")
#sample info
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#first combine all results from bamreadcount into one matrix
get_bam_readcount = function(file_res){

	#samplename
	samplename = unlist(strsplit(file_res, "_missing_muts"))[1]

	#readcount output
	x = file_res

	#get output
	output_t1reads_t2muts = as.data.table(bam_readcount.parse(x, samplename = samplename, "pyclone_bam_readcount_input.bed"))

	#save output
	return(output_t1reads_t2muts)
}

#get files based on whether it was on all mutations or subset of mutations
#bam readcount values for mutations not found in all patients

get_reads = function(type_analysis){

	#mutations that were only found in one sample
	unique = filter(as.data.table(table(read_only$mut_id)), N == 1)

	if(type_analysis == "full"){
		bamreadcount=list.files(pattern="_missing_muts_all")
		pyclone_input=as.data.table(filter(read_only, !(mut_id %in% unique$V1),
		MajorCN > 0))
		#file with missing variants
		miss_vars = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_muts_pyclone_bam_readcount_input.bed")
		colnames(miss_vars) = c("chr", "start", "end",
																			 "ref_allele", "alt_allele")
	  write.table(miss_vars, "pyclone_bam_readcount_all_muts_input.bed", quote=F, row.names=F, sep="\t")
	}

	if(type_analysis == "subset"){
		bamreadcount=list.files(pattern="_missing_muts_small")
		pyclone_input = as.data.table(filter(read_only, !(mut_id %in% unique$V1),
		Copy_Number >=2, MajorCN > 0))
		#file with missing variants
		miss_vars = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/small_subset_pyclone_bam_readcount_input.bed")
		colnames(miss_vars) = c("chr", "start", "end",
																			 "ref_allele", "alt_allele")
	  write.table(miss_vars, "pyclone_bam_readcount_some_muts_input.bed", quote=F, row.names=F, sep="\t")
	}

	missing_mutations = as.data.table(ldply(llply(bamreadcount, get_bam_readcount)))
	missing_mutations$id = paste(missing_mutations$chr, missing_mutations$start, sep="_")

	#get list of mutations that are present in all
	t = filter(as.data.table(table(pyclone_input$mut_id)), (N==20))
	muts_all = as.data.table(filter(pyclone_input, mut_id %in% t$V1))
	muts_some = as.data.table(filter(pyclone_input, !(mut_id %in% t$V1)))

	#for mutations that are found in only some patients
	#extract mutation info for those mutations across all samples

	get_record = function(mutation){
	  print(mutation)
	  mut_dat = as.data.table(filter(missing_mutations, id == mutation))
		pyclone_dat_mut = as.data.table(filter(pyclone_input, mut_id == mutation))
	  #columns that we need to edit
	  #c("mut_id", "Ref_counts", "alt_counts", "normal_cn", "Nmin", "Nmaj",
		#"hg19.ensemblToGeneName.value", "Func.ensGene", "id")
	  #which patients don't have
		pats = unique(mut_dat$samplename)
		pats_missing = unique(read_only$Indiv)[which(!(unique(read_only$id) %in% pyclone_dat_mut$id))]

	    make_pat = function(pat){
	      pat_dat = filter(mut_dat, samplename == pat)
				pat_dat$mut_id = mutation
				pat_dat$id = read_only$id[which(read_only$Indiv == pat)[1]]

					if(pat %in% pats_missing){
							pat_dat$MajorCN=2
	      			pat_dat$MinorCN=0
							pat_dat$source = "bamreadcount"
							}

					if(!(pat %in% pats_missing)){
						pat_dat$MajorCN=filter(read_only, mut_id==mutation, Indiv == pat)$MajorCN
						pat_dat$MinorCN=filter(read_only, mut_id==mutation, Indiv == pat)$MinorCN
						pat_dat$source = "mutation_called"
					}

				pat_dat$Ref_counts=pat_dat$ref_count
	      pat_dat$alt_counts=pat_dat$alt_count
				pat_dat$symbol = read_only$symbol[which(read_only$mut_id == mutation)[1]]
				pat_dat$Func.ensGene = read_only$Func.ensGene[which(read_only$mut_id == mutation)[1]]
	      return(pat_dat)
	    }

	  all_pats = as.data.table(ldply(llply(pats, make_pat, .progress="text")))
	  return(all_pats)
	}

	all_muts = unique(muts_some$mut_id)
	missing_records = as.data.table(ldply(llply(all_muts, get_record, .progress="text")))

	#combine with mutation data for mutations that were called in all samples
	#need the same columns
	z = which(colnames(muts_all) %in% colnames(missing_records))
	muts_all = muts_all[,..z]

	z = which(colnames(missing_records) %in% colnames(muts_all))
	missing_records = missing_records[,..z]

	#now make them match wtih missing records columns so that can bind them together
	cols  = colnames(missing_records)
	muts_all = muts_all[,..cols]

	all_records = rbind(muts_all, missing_records)
	write.table(all_records, file=paste(date, type_analysis, "mutations", "PYCLONE_INPUT_MUTS.txt", sep="_"), quote=F, row.names=F, sep="\t")

	print(length(unique(all_records$mut_id)))
	print("done")
}

get_reads("full")
get_reads("subset")
