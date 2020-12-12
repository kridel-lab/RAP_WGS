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
read_only = readRDS(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.rds")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.rds"))])

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone")
#sample info
samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")
colnames(samps)[4] ="Indiv"
z = which(samps$Indiv %in% read_only$Indiv)
samps = samps[z,]

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))

#save sample id versus sample name clean
patients= c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#first combine all results from bamreadcount into one matrix
get_bam_readcount = function(file_res, type_analysis){

	print(file_res)

	#samplename
	samplename = unlist(strsplit(file_res, "_missing_muts"))[1]

	#readcount output
	x = file_res

	#mut file
	if(type_analysis == "full"){
		muts_file="pyclone_bam_readcount_all_muts_input.bed"
	}

	if(type_analysis == "subset"){
		muts_file="pyclone_bam_readcount_some_muts_input.bed"
	}

	#get output
	output_t1reads_t2muts = as.data.table(bam_readcount.parse(x, samplename = samplename,
		muts_file))

	#save output
	return(output_t1reads_t2muts)
}

#get files based on whether it was on all mutations or subset of mutations
#bam readcount values for mutations not found in all patients

get_reads = function(patient, type_analysis){

	print(patient)

	#mutations that were only found in one sample
	unique = filter(as.data.table(table(filter(read_only, STUDY_PATIENT_ID==patient)$mut_id)), N == 1)
	read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID==patient))

	#for clonal evolution analysis only keep founder mutations that are in DLBCL genes
	if(patient=="LY_RAP_0001"){
		founds_keep = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 3, V2=="yes")
		founds_remove = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 3, V2=="")
		miss_vars_full = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_full_LY_RAP_0001_pyclone_bam_readcount_input.bed")
		miss_vars_small = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_small_subset_LY_RAP_0001_pyclone_bam_readcount_input.bed")
	}

	if(patient=="LY_RAP_0002"){
		founds_keep = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 4, V2=="yes")
		founds_remove = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 4, V2=="")
		miss_vars_full = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_full_LY_RAP_0002_pyclone_bam_readcount_input.bed")
		miss_vars_small = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_small_subset_LY_RAP_0002_pyclone_bam_readcount_input.bed")
	}

	if(patient=="LY_RAP_0003"){
		founds_keep = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 20, V2=="yes")
		founds_remove = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 4, V2=="")
		miss_vars_full = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_full_LY_RAP_0003_pyclone_bam_readcount_input.bed")
		miss_vars_small = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_small_subset_LY_RAP_0003_pyclone_bam_readcount_input.bed")
	}

	if(type_analysis == "full"){
		bamreadcount=list.files(pattern="_missing_muts_all")
		z = which(str_detect(bamreadcount, patient))
		bamreadcount = bamreadcount[z]

		pyclone_input = as.data.table(filter(read_only_pat, !(mut_id %in% unique$V1),
		!(mut_id %in% founds_remove$V1), tot_counts >=60, gt_AF >=0.15,
		MajorCN > 0, Copy_Number >=2, Copy_Number <=8))
		print(length(unique(pyclone_input$mut_id)))

		#file with missing variants
		miss_vars = miss_vars_full
		colnames(miss_vars) = c("chr", "start", "end",
																			 "ref_allele", "alt_allele")
	  write.table(miss_vars, "pyclone_bam_readcount_all_muts_input.bed", quote=F, row.names=F, sep="\t")
	}

	if(type_analysis == "subset"){
		bamreadcount=list.files(pattern="_missing_muts_small")
		z = which(str_detect(bamreadcount, patient))
		bamreadcount = bamreadcount[z]

		pyclone_full = as.data.table(filter(read_only_pat, !(mut_id %in% unique$V1),
		!(mut_id %in% founds_remove$V1), tot_counts >=60,gt_AF >=0.15,
		MajorCN > 0, Copy_Number >=2, Copy_Number <=8))

		pyclone_input = as.data.table(filter(pyclone_full,
		!(Func.ensGene %in% c("ncRNA_intronic", "intergenic", "intronic"))))
		length(unique(pyclone_input$mut_id))

		#file with missing variants
		miss_vars = miss_vars_small
		colnames(miss_vars) = c("chr", "start", "end",
																			 "ref_allele", "alt_allele")
	  write.table(miss_vars, "pyclone_bam_readcount_some_muts_input.bed", quote=F, row.names=F, sep="\t")
	}

	missing_mutations = as.data.table(ldply(llply(bamreadcount, get_bam_readcount, type_analysis)))
	missing_mutations$id = paste(missing_mutations$chr, missing_mutations$start, sep="_")
	missing_mutations$patient = patient

	#get list of mutations that are present in all
	if(patient == "LY_RAP_0001"){
		t = filter(as.data.table(table(pyclone_input$mut_id)), (N==3))
		muts_all = as.data.table(filter(pyclone_input, mut_id %in% t$V1))
		muts_some = as.data.table(filter(pyclone_input, !(mut_id %in% t$V1)))
	}
	if(patient == "LY_RAP_0002"){
		t = filter(as.data.table(table(pyclone_input$mut_id)), (N==4))
		muts_all = as.data.table(filter(pyclone_input, mut_id %in% t$V1))
		muts_some = as.data.table(filter(pyclone_input, !(mut_id %in% t$V1)))
	}
	if(patient == "LY_RAP_0003"){
		t = filter(as.data.table(table(pyclone_input$mut_id)), (N==20))
		muts_all = as.data.table(filter(pyclone_input, mut_id %in% t$V1))
		muts_some = as.data.table(filter(pyclone_input, !(mut_id %in% t$V1)))
	}

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
		pats_missing = unique(read_only_pat$Indiv)[which(!(unique(read_only_pat$id) %in% pyclone_dat_mut$id))]

	    make_pat = function(pat){
	      pat_dat = filter(mut_dat, samplename == pat)
				pat_dat$mut_id = mutation
				pat_dat$id = read_only_pat$id[which(read_only_pat$Indiv == pat)[1]]

					if(pat %in% pats_missing){
							pat_dat$MajorCN=2
	      			pat_dat$MinorCN=0
							pat_dat$source = "bamreadcount"
							}

					if(!(pat %in% pats_missing)){
						pat_dat$MajorCN=filter(read_only_pat, mut_id==mutation, Indiv == pat)$MajorCN
						pat_dat$MinorCN=filter(read_only_pat, mut_id==mutation, Indiv == pat)$MinorCN
						pat_dat$source = "mutation_called"
					}

				pat_dat$Ref_counts=pat_dat$ref_count
	      pat_dat$alt_counts=pat_dat$alt_count
				pat_dat$symbol = read_only_pat$symbol[which(read_only_pat$mut_id == mutation)[1]]
				pat_dat$Func.ensGene = read_only_pat$Func.ensGene[which(read_only_pat$mut_id == mutation)[1]]
	      return(pat_dat)
	    }

	  all_pats = as.data.table(ldply(llply(pats, make_pat, .progress="text")))
	  return(all_pats)
	}

	all_muts = unique(muts_some$mut_id)
	missing_records = as.data.table(ldply(llply(all_muts, get_record, .progress="text")))

	#combine with mutation data for mutations that were called in all samples
	#need the same columns, make sure sample column included
	colnames(muts_all)[1] = "samplename"
	z = which(colnames(muts_all) %in% colnames(missing_records))
	muts_all = muts_all[,..z]

	z = which(colnames(missing_records) %in% colnames(muts_all))
	missing_records = missing_records[,..z]

	#now make them match wtih missing records columns so that can bind them together
	cols  = colnames(missing_records)
	muts_all = muts_all[,..cols]

	all_records = rbind(muts_all, missing_records)
	write.table(all_records, file=paste(date, type_analysis, patient, "mutations", "PYCLONE_INPUT_MUTS.txt", sep="_"), quote=F, row.names=F, sep="\t")

	print(length(unique(all_records$mut_id)))
	print("done")
}

get_full_missing_data = function(patient){
	print(patient)
	get_reads(patient, "full")
	get_reads(patient, "subset")
	print("done")
}

#run analysis
llply(patients, get_full_missing_data)
