#----------------------------------------------------------------------
#karin isaev
#make SNV list of dataframes as input for sciclone
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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#vaf data frame:
#1. chromosome
#2. position
#3. reference-supporting read counts
#4. variant-supporting read counts
#5. variant allele fraction (between 0-100)

#sampleNames: vector of names describing each sample ex: ("Primary
          #Tumor", "Relapse")

#1. Mutation file produced using bam readcount - mutation accounted for every patient even if not called
muts_pyclone = fread("2020-09-13_full_mutations_PYCLONE_INPUT_MUTS.txt")

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone")

#split dataframe into list by sample
#save order of samples - will need to keep order same for other sciclone input
list_df <- split(muts_pyclone, f = muts_pyclone$id)

#save order of samples
samps_order = names(list_df)
saveRDS(samps_order, "snvs_sciclone_list_sample_order.rds")

#for each element in list remove original sample column
for(i in 1:length(list_df)){
  list_df[[i]] = as.data.frame(list_df[[i]])
  list_df[[i]]$id = NULL
	list_df[[i]]$vaf = list_df[[i]]$alt_counts/(list_df[[i]]$Ref_counts + list_df[[i]]$alt_counts)*100
	list_df[[i]] = list_df[[i]][,c("chr", "start", "Ref_counts", "alt_counts", "vaf")]
}

#save full dataset
saveRDS(list_df, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_snvs_df.rds")

#create annotation file
#annotation:
#1) chromosome
#2) position
#3) gene name.

anno = unique(muts_pyclone[,c("chr", "start", "symbol")])
saveRDS(anno, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/annotation_file.rds")
