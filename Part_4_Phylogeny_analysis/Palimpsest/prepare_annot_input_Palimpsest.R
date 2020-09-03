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
clusters$Gender = "M"
clusters = clusters[,c("Sample", "Gender", "purity")]
colnames(clusters)[3] = "Purity"

write.table(clusters, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Palimpsest/input/annotation_data_palimpsest_input.txt", quote=F, row.names=F, sep="\t")
