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
	"stringr", "readxl", "GenomicRanges", "params", "readr")
lapply(packages, require, character.only = TRUE)

args = commandArgs(trailingOnly = TRUE) #patient ID
index = args[1]
print(index)
patient = index
print(patient)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#used bam readcount to extract counts from all samples
#overlapping mutations that were only found in some samples and not all
#Pyclone input requires these values

#test
#patient="LY_RAP_0001"

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

source("/cluster/home/kisaev/RAP_WGS/config-file.R")

#save sample id versus sample name clean
patients= c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone/")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

get_major_cn = function(patient, mut){
  print(patient)
  cnas_pat = as.data.table(filter(cnas, Sample == patient))
	snvs = filter(read_only, mut_id == mut)[1,]
	#just save coordiantes of mutations
	snvs = snvs[,c("CHROM", "POS")]

  #make granges objects
  snvs_gr = GRanges(
    seqnames = snvs$CHROM,
    ranges = IRanges(snvs$POS, end = snvs$POS),
    strand = rep("*", length(snvs$POS)),
    score = 1:length(snvs$POS))

  cnas_gr = GRanges(
    seqnames = cnas_pat$CHROM,
    ranges = IRanges(cnas_pat$Start, end = cnas_pat$End),
    strand = rep("*", length(cnas_pat$Start)),
    logR = cnas_pat$depth.ratio)

  #intersect them
  #Then subset the original objects with the negative indices of the overlaps:
  hits <- findOverlaps(snvs_gr, cnas_gr, ignore.strand=TRUE)
  hits_overlap = cbind(snvs[queryHits(hits),], cnas_pat[subjectHits(hits),])
  major_cn = hits_overlap$Nmaj[1]
  return(major_cn)
}

get_minor_cn = function(patient, mut){
  print(patient)
  cnas_pat = as.data.table(filter(cnas, Sample == patient))
	snvs = filter(read_only, mut_id == mut)[1,]
	#just save coordiantes of mutations
	snvs = snvs[,c("CHROM", "POS")]

  #make granges objects
  snvs_gr = GRanges(
    seqnames = snvs$CHROM,
    ranges = IRanges(snvs$POS, end = snvs$POS),
    strand = rep("*", length(snvs$POS)),
    score = 1:length(snvs$POS))

  cnas_gr = GRanges(
    seqnames = cnas_pat$CHROM,
    ranges = IRanges(cnas_pat$Start, end = cnas_pat$End),
    strand = rep("*", length(cnas_pat$Start)),
    logR = cnas_pat$depth.ratio)

  #intersect them
  #Then subset the original objects with the negative indices of the overlaps:
  hits <- findOverlaps(snvs_gr, cnas_gr, ignore.strand=TRUE)
  hits_overlap = cbind(snvs[queryHits(hits),], cnas_pat[subjectHits(hits),])
  minor_cn = hits_overlap$Nmin[1]
  return(minor_cn)
}

#first combine all results from bamreadcount into one matrix
get_bam_readcount = function(file_res, patient){

	print(file_res)

	#samplename
	samplename = unlist(strsplit(file_res, "_missing_muts"))[1]

	#readcount output
	x = file_res

	#mut file
	muts_file=paste(patient, "pyclone_bam_readcount_all_muts_input.bed", sep="_")

	#get output
	output_t1reads_t2muts = as.data.table(bam_readcount.parse(x, samplename = samplename,
		muts_file))

	#save output
	return(output_t1reads_t2muts)
}

#get files based on whether it was on all mutations or subset of mutations
#bam readcount values for mutations not found in all patients

get_reads = function(patient){

	print(patient)

	read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID==patient))

	#for clonal evolution analysis only keep founder mutations that are in DLBCL genes
	if(patient=="LY_RAP_0001"){
		miss_vars_full = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_mutations_LY_RAP_0001_pyclone_bam_readcount_input.bed")
	}

	if(patient=="LY_RAP_0002"){
		miss_vars_full = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_mutations_LY_RAP_0002_pyclone_bam_readcount_input.bed")
	}

	if(patient=="LY_RAP_0003"){
		miss_vars_full = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_mutations_LY_RAP_0003_pyclone_bam_readcount_input.bed")
	}

	bamreadcount=list.files(pattern="_missing_muts_all")
	z = which(str_detect(bamreadcount, patient))
	bamreadcount = bamreadcount[z]

	#file with missing variants
	miss_vars = miss_vars_full
	colnames(miss_vars) = c("chr", "start", "end",
																			 "ref_allele", "alt_allele")
	write.table(miss_vars, paste(patient, "pyclone_bam_readcount_all_muts_input.bed", sep="_"), quote=F, row.names=F, sep="\t")

	missing_mutations = as.data.table(ldply(llply(bamreadcount, get_bam_readcount, patient)))
	missing_mutations$id = paste(missing_mutations$chr, missing_mutations$start, sep="_")
	missing_mutations$patient = patient

	z = which(is.na(missing_mutations$samplename))
	if(!(length(z) == 0)){
		missing_mutations = missing_mutations[-z,]
	}

	#make sure muts appear in all samples now
	t=as.data.table(table(missing_mutations$id))
	print(unique(t$N))

	muts_keep = filter(t, N ==3)$V1
	missing_mutations = filter(missing_mutations, id %in% muts_keep)

	#for mutations that are found in only some patients
	#extract mutation info for those mutations across all samples
	get_record = function(mutation){
	  print(mutation)
	  mut_dat = as.data.table(filter(missing_mutations, id == mutation))
		pyclone_dat_mut = as.data.table(filter(read_only_pat, mut_id == mutation))
		pats = unique(mut_dat$samplename)
		pats_missing = pats[which(!(pats %in% read_only_pat$Indiv[which(read_only_pat$mut_id %in% mutation)]))]

	    make_pat = function(pat){
	      pat_dat = filter(mut_dat, samplename == pat)
				pat_dat$mut_id = mutation
				#pat_dat$id = read_only_pat$id[which(read_only_pat$Indiv == pat)[1]]

					if(pat %in% pats_missing){
							pat_dat$MajorCN=get_major_cn(pat, mutation)
	      			pat_dat$MinorCN=get_minor_cn(pat, mutation)
							pat_dat$source = "bamreadcount"
							}

					if(!(pat %in% pats_missing)){
						pat_dat$MajorCN=filter(read_only_pat, mut_id==mutation, Indiv == pat)$Nmaj
						pat_dat$MinorCN=filter(read_only_pat, mut_id==mutation, Indiv == pat)$Nmin
						pat_dat$source = "mutation_called"
					}

				pat_dat$Ref_counts=pat_dat$ref_count
	      pat_dat$alt_counts=pat_dat$alt_count
				pat_dat$symbol = read_only_pat$symbol[which(read_only_pat$mut_id == mutation)[1]]
				pat_dat$Func.ensGene = read_only_pat$Func.ensGene[which(read_only_pat$mut_id == mutation)[1]]
	      return(pat_dat)
	    }

	  all_pats = as.data.table(ldply(llply(pats, make_pat)))
	  return(all_pats)
	}

	all_muts = unique(missing_mutations$id)
	missing_records = as.data.table(ldply(llply(all_muts, get_record, .progress="text")))
	all_records = missing_records
	write.table(all_records, file=paste(date, patient, "mutations", "PYCLONE_INPUT_MUTS.txt", sep="_"), quote=F, row.names=F, sep="\t")

	print(length(unique(all_records$mut_id)))
	print("done")
}

get_full_missing_data = function(patient){
	print(patient)
	get_reads(patient)
	print("done")
}

#run analysis
get_full_missing_data(patient)
