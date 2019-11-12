#----------------------------------------------------------------------
#karin isaev
#oct 30, 2019
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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#here:
#1. remove blacklisted region mutations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

files = list.files(pattern="vcf_file_filtered.bed")

#read all files and rbind 
all_muts = as.data.table(ldply(llply(files, function(x){fread(x)})))
all_muts_cords = unique(all_muts[,c("CHROM", "POS", "mut_id")])
all_muts_cords$POSend = all_muts_cords$POS
all_muts_cords$CHROM = paste("chr", all_muts_cords$CHROM, sep="")
all_muts_cords = all_muts_cords[,c("CHROM", "POS", "POSend", "mut_id")]
all_muts_cords=all_muts_cords[order(CHROM, POS)]
write.table(all_muts_cords, file="all_muts_pre_blacklist_filter.bed", quote=F, row.names=F, sep="\t", col.names=F)

#DO NOT RUN 
#module load bedtools
#bedtools intersect -wa -wb \
#    -a all_muts_pre_blacklist_filter.bed \
#    -b hg19-blacklist.v2.bed -f 1.0 > all_muts_post_blacklist_filter.bed \

blacklist_muts = fread("all_muts_post_blacklist_filter.bed")
colnames(blacklist_muts)[1:4] = c("CHROM", "POS", "POS_end", "mut_id")
z = which(blacklist_muts$mut_id %in% all_muts$mut_id) #643 mutations 
z = which(all_muts$mut_id %in% blacklist_muts$mut_id) #643 mutations 
all_muts = all_muts[-z,]

#other filters to apply to these variants? 

#split into num ref and alt allele counts 
all_muts = all_muts %>% separate(gt_AD, c("Ref_counts", "alt_counts"))
all_muts$Ref_counts = as.numeric(all_muts$Ref_counts)
all_muts$alt_counts = as.numeric(all_muts$alt_counts)
#remove empty columns
all_muts$gt_DP = NULL
all_muts$gt_GQ = NULL
all_muts$gt_PGT = NULL
all_muts$gt_PID = NULL
all_muts$gt_PL = NULL
all_muts$ID = NULL
all_muts$QUAL = NULL
all_muts$AC = NULL
all_muts$AF = NULL
all_muts$AN = NULL
all_muts$IN_PON = NULL
all_muts$NLOD = NULL
all_muts$N_ART_LOD = NULL
all_muts$RPA = NULL
all_muts$RU = NULL
all_muts$STR = NULL
all_muts$ALLELE_END = NULL

#TLOD - min 10 - confidence score for variant being really somatic 
all_muts = as.data.table(filter(all_muts, TLOD > 10)) #64970 unique mutations 

#min counts for alternative reads 
all_muts = as.data.table(filter(all_muts, alt_counts > 10))
length(unique(all_muts$mut_id))

#save file so far --> download locally and plot summary stats for various variable 
length(unique(all_muts$mut_id))
saveRDS(all_muts, file=paste(date, "all_muts_merged_mut_calls.rds", sep="_")) #this file is used for downstream exploratory analysis 






