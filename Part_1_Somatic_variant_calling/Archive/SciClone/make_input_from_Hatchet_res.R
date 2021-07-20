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
	"ggrepel", "stringr", "maftools", "sciClone")
lapply(packages, require, character.only = TRUE)

setwd("/Users/kisaev/Documents/Hatchet_analysis")

muts = fread("output.txt")
list_df <- split(muts, f = muts$'PATIENT-SAMPLE')

#save order of samples
samps_order = names(list_df)
saveRDS(samps_order, "snvs_sciclone_list_sample_order.rds")

#for each element in list remove original sample column
for(i in 1:length(list_df)){
  list_df[[i]] = as.data.frame(list_df[[i]])
  list_df[[i]]$id = NULL
	list_df[[i]]$vaf = (list_df[[i]]$CCF/2) * 100
  list_df[[i]]$Ref_counts = as.numeric(sapply(list_df[[i]]$COUNTS, function(x){unlist(strsplit(x, ","))[1]}))
  list_df[[i]]$alt_counts = as.numeric(sapply(list_df[[i]]$COUNTS, function(x){unlist(strsplit(x, ","))[2]}))
	list_df[[i]] = list_df[[i]][,c('#CHR', "POS", "Ref_counts", "alt_counts", "vaf")]
}

#save full dataset
saveRDS(list_df, file="list_snvs_df.rds")
