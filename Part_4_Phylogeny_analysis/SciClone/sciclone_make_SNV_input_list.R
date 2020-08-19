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

#1. READ-ONLY FILE = FINAL MUTATIONS = KEEP THIS WAY UNLESS MAJOR CHANGE NEEDED
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])
read_only = read_only[,c("Sample", "CHROM", "POS", "Ref_counts", "alt_counts", "gt_AF")]
read_only = read_only[,c("Sample", "CHROM", "POS", "Ref_counts", "alt_counts", "gt_AF")]

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone")

#split dataframe into list by sample
#save order of samples - will need to keep order same for other sciclone input
list_df <- split(read_only, f = read_only$Sample)

#for each element in list remove original sample column
for(i in 1:length(list_df)){
  list_df[[i]] = as.data.frame(list_df[[i]])
  list_df[[i]]$Sample = NULL
}

#save full dataset
saveRDS(list_df, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_snvs_df.rds")
