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
	"ggrepel", "stringr", "maftools", "readxl")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare mutations in SNV format required for palimpsest to run

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#Our mutations
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

#sample info
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#prepare SNV input file:

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

input = read_only %>%
	mutate(Type = "SNV", Tumor_Varcount = alt_counts,
Tumor_Depth = DP, Normal_Depth = Ref_counts, Gene_Name = symbol, Driver = "") %>%
	select(Sample, Type, CHROM, POS, REF, ALT, Tumor_Varcount, Tumor_Depth,  Normal_Depth, Gene_Name, Driver)

input$Driver[which(input$Gene_Name %in% reddy$Gene)] = input$Gene_Name[which(input$Gene_Name %in% reddy$Gene)]

saveRDS(input, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Palimpsest/SNV_input_Palimpsest.rds")
