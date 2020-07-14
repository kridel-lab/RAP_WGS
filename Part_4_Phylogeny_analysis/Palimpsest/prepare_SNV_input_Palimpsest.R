#----------------------------------------------------------------------
#karin isaev
#Nov 1st, 2019
#----------------------------------------------------------------------

date = Sys.Date()

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_RESULTS/strelka_filtered")

#output from TitanCNA, combine all samples into one dataframe
#use optimalClusterSolution.txt file to identify optimal cluster for each sample
#use sample to identifier conversion to get actual sample name

#1]. mut_data: somatic mutation data

#Sample: Sample identifier. Any alphanumeric string.
#Type: Mutation type [SNV/Indel].
#CHROM: Chromosome. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
#POS: Mutation position. A positive integer.
#REF: Reference base(s): Each base must be one of A,C,G,T (upper case). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the string.
#ALT: Alternate base(s): Each base must be one of A,C,G,T (upper case). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the string.
#Tumor_Varcount: Number of variant bases at the position in the tumor sample.
#Tumor_Depth: Tumor sample sequencing depth at the position.
#Normal_Depth: Normal sample sequencing depth at the position.
#Gene_Name: OPTIONAL column for representing mutated gene name.
#Driver: OPTIONAL column indicating the driver events to be annotated in tumor history plots.

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

vcfs = list.files(pattern="no_info")

#args = commandArgs(trailingOnly = TRUE)
#index = as.integer(args[1])

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#final muts
final_muts = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/2019-11-01_final_list_of_mutations_input_palimpsest.rds")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(vcf){

  pat = paste(unlist(strsplit(vcf, "_"))[4:9], collapse="_")
  patient_muts = as.data.table(filter(final_muts, Indiv == pat))

  mutations_T1 =read.vcfR(vcf)
  mutations_T1 = vcfR2tidy(mutations_T1)
  meta = as.data.table(mutations_T1$meta)

  vcf = as.data.table(mutations_T1$fix)
  chr_conv = unique(vcf[,c("ChromKey", "CHROM")])

  gt = as.data.table(mutations_T1$gt)
  gt = merge(gt, chr_conv, by="ChromKey")

  gt$mut_id = paste(gt$CHROM, gt$POS, sep="_")
  vcf$mut_id = paste(vcf$CHROM, vcf$POS, sep="_")

  #1. keep only the ones that passed default mutect2 filter
  vcf = as.data.table(filter(vcf, FILTER=="PASS"))
  print(paste("number of variants that passed filtering=", dim(vcf)[1]))

  #2. combine vcf and gt info
  cols = colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt = merge(gt, vcf, by= cols)

  #4. remove variants from chromosome X and Y
  gt$CHROM = as.numeric(gt$CHROM)
  gt = as.data.table(filter(gt, (CHROM %in% c(1:22))))
  gt = as.data.table(filter(gt, mut_id %in% patient_muts$mut_id))

  #Tumor_Varcount: Number of variant bases at the position in the tumor sample.
  #Tumor_Depth: Tumor sample sequencing depth at the position.
  #Normal_Depth: Normal sample sequencing depth at the position.
  gt$CHROM = paste("chr", gt$CHROM, sep="")
  gt$Type = "SNV"
  gt$patient = pat

  get_mut_dat = function(mut){
	mut_dat = gt[,c("Indiv", "Type", "CHROM", "POS", "REF", "ALT", "mut_id", "gt_AU", "gt_CU", "gt_DP", "gt_GU", "gt_TU", "patient")]
	mut_dat = as.data.frame(t(as.data.table(filter(mut_dat, mut_id == mut))))
	colnames(mut_dat) = mut_dat[1,]
	mut_dat = mut_dat[-1,]
	Tumor_Depth = mut_dat$TUMOR[which(rownames(mut_dat) == "gt_DP")]
	#get tumour varcount and normal depth
	Normal_Depth = mut_dat$NORMAL[which(rownames(mut_dat) == "gt_DP")]
	alt = mut_dat$TUMOR[which(rownames(mut_dat) == "ALT")]
	alt = paste("gt_", alt, "U", sep="")
	Tumor_Varcount = unlist(strsplit(mut_dat[which(rownames(mut_dat)==alt),which(colnames(mut_dat)=="TUMOR")], ","))[1]
	mut_dat = gt[,c("patient", "Indiv", "Type", "CHROM", "POS", "REF", "ALT", "mut_id", "gt_AU", "gt_CU", "gt_DP", "gt_GU", "gt_TU")]
	mut_dat = as.data.frame(as.data.table(filter(mut_dat, mut_id == mut)))
	mut_dat$Tumor_Varcount = Tumor_Varcount
	mut_dat$Tumor_Depth = Tumor_Depth
	mut_dat$Normal_Depth = Normal_Depth
	mut_dat = as.data.table(filter(mut_dat, Indiv == "TUMOR"))
	mut_dat = mut_dat[,c("patient", "Type", "CHROM", "POS", "REF", "ALT", "Tumor_Varcount", "Tumor_Depth", "Normal_Depth")]
	colnames(mut_dat)[1] = "Sample"
	return(mut_dat)
  }

  mut_info = as.data.table(ldply(llply(unique(gt$mut_id), get_mut_dat, .progress="text")))
  return(mut_info)
}

all_muts = as.data.table(ldply(llply(vcfs, clean_up_001, .progress="text")))
saveRDS(all_muts, file="SNV_input_Palimpsest.rds")
