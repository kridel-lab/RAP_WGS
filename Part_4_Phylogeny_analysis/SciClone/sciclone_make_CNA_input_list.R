#----------------------------------------------------------------------
#karin isaev
#make CNA list of dataframes as input for sciclone
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

#copyNumber data frame:
#1. chromosome
#2. segment start position
#3. segment stop position
#4. copy number value for that segment.
#Unrepresented regions are assumed to have a copy number of 2.

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

#isolate LOH regions
#hemizygous deletion LOH (DLOH),
#copy neutral LOH (NLOH),
#amplified LOH (ALOH)
z = which(all_cnas$TITAN_call %in% c("NLOH"))
loh = all_cnas[z,]
all_cnas = all_cnas[-z,]

#save
all_cnas = all_cnas[,c("Sample", "CHROM", "Start", "End",
"Copy_Number")]

#split dataframe into list by sample
#save order of samples - will need to keep order same for other sciclone input
list_df <- split(all_cnas, f = all_cnas$Sample)

#for each element in list remove original sample column
for(i in 1:length(list_df)){
  list_df[[i]] = as.data.frame(list_df[[i]])
  list_df[[i]]$Sample = NULL
}

#save full dataset
saveRDS(list_df, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_cnas_df.rds")

#prepare list of regions to exclude = LOH
# 1. chromosome, 2. window start position 3. window stop position;
loh = unique(loh[,c("Sample", "CHROM", "Start", "End")])

list_loh <- split(loh, f = loh$Sample)

#for each element in list remove original sample column
for(i in 1:length(list_loh)){
  list_loh[[i]] = as.data.frame(list_loh[[i]])
  list_loh[[i]]$Sample = NULL
}

#save full dataset as list
saveRDS(list_loh, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_loh_df.rds")

#save full as dataset as one file
loh$Sample = NULL
loh = unique(loh)
loh = as.data.frame(loh)
saveRDS(loh, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/loh_full_df.rds")
