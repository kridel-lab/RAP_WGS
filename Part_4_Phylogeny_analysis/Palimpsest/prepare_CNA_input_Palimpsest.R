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

#save sample cluster and purity info and upload to files
cna_save = clusters[,c("Sample", "numClust", "cellPrev",
"purity", "norm", "ploidy")]

write.table(cna_save, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/TitanCNA_summary.txt",
quote=F, row.names=F, sep=";")

#titanCNA results
files = list.files(pattern="seg.txt")
files = files[sapply(clusters$id, function(x){which(str_detect(files, x))})]

#read in data files
all_cnas = as.data.table(ldply(llply(files, function(x){fread(x)})))
colnames(all_cnas)[1] = "barcode"
all_cnas = merge(all_cnas, clusters, by = "barcode")
all_cnas$CHROM = paste("chr", all_cnas$Chromosome, sep="")

#save for filtering SNVs
all_cnas = all_cnas[,c("Sample", "CHROM", "Start", "End",
"logR_Copy_Number", "MinorCN", "MajorCN", "TITAN_state",
"Copy_Number", "Corrected_Call", "ploidy")]
#save full dataset
saveRDS(all_cnas, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_TITAN.rds")

all_cnas_palimpsest = all_cnas[,c("Sample", "CHROM", "Start", "End",
"logR_Copy_Number", "MinorCN", "MajorCN",
"Copy_Number", "ploidy")]
colnames(all_cnas_palimpsest) = c("Sample", "CHROM", "POS_START",
"POS_END", "LogR", "Nmin", "Nmaj", "ntot", "Ploidy")
write.table(all_cnas_palimpsest, file="copy_number_alteration_data_palimpsest_input.txt", quote=F, row.names=F, sep="\t")
