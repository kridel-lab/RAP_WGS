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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/TITAN_CNA/results/titan/hmm/optimalClusterSolution_files/titanCNA_ploidy2")

#output from TitanCNA, combine all samples into one dataframe 
#use optimalClusterSolution.txt file to identify optimal cluster for each sample
#use sample to identifier conversion to get actual sample name 

#2]. cna_data: copy number alteration data

#Sample: Sample identifier. Any alphanumeric string.
#CHROM: Chromosome. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
#POS_START: Start position of segmented chromosome.
#POS_END: End position of segmented chromosome.
#LogR: LogR information.
#Nmin: Minor allele copy number.
#Nmaj: Major allele copy number.
#ntot: Total copy number of segmented chromosome.
#Ploidy: Tumor ploidy.

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#sample conversion
samples = fread("cna_vcf_sample_conversion.csv")
colnames(samples) = c("barcode", "Sample")

#optimal clusters 
clusters = fread("optimalClusterSolution.txt")
clusters= merge(samples, clusters, by = "barcode")

#titanCNA results 
files = list.files(pattern="segs.txt")
files = files[sapply(clusters$id, function(x){which(str_detect(files, x))})]

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#make CNA input for phyloWGS 
make_cna = function(file){
	
	#read in data files 
	all_cnas = as.data.table(fread(file))
	all_cnas$End_Position.bp. = all_cnas$End_Position.bp. + 1
	print(head(all_cnas))
	colnames(all_cnas)[3] = "Start_Position(bp)"
	colnames(all_cnas)[4] = "End_Position(bp)"
	colnames(all_cnas)[14] = "Clonal_Frequency"
	sample = all_cnas$Sample[1]
	print(sample)
	all_cnas$Sample = clusters$Sample[which(clusters$barcode == sample)]
	write.table(all_cnas, file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/CNAs_titan/", all_cnas$Sample[1], "cna_datat_input_phyloWGS_CNV_parser.txt", sep="_"), quote=F, row.names=F, sep="\t")
	print("done")
}

llply(files, make_cna)

