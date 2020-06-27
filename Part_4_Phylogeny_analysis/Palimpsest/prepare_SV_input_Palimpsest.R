#----------------------------------------------------------------------
#karin isaev
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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Palimpsest")

#using output from MANRA

#4]. sv_data: structural variant data

#    Sample: Sample identifier. Any alphanumeric string.
#    Type: Type of structural variant: INV/DEL/DUP/BND.
#    CHROM_1: Chromosome of the first breakpoint. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
#    POS_1: Position of the first breakpoint. A positive integer.
#    CHROM_2: Chromosome of the second breakpoint. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
#    POS_2: Position of the second breakpoint. A positive integer.
#    Tumor_Varcount: Column for variant allele count information.
#    Tumor_Depth: Column for tumour sequencing depth information.
#    Normal_Depth: Column for normal sequencing depth information.
#    Driver: OPTIONAL column indicating the driver events to be annotated in tumour history plots.

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

svs = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/2019-12-17_RAP_WGS_MANTA_masterlist.csv")

svs_palimpsest = svs[,c("patient", "sv1_SVTYPE", "SV_CHR", "SV_start",
"CHROM_2", "POS_2", "Tumor_Varcount",
"Tumor_Depth", "Normal_Depth", "gene")]

colnames(svs_palimpsest) = c("Sample", "Type", "CHROM_1", "POS_1",
"CHROM_2", "POS_2", "Tumor_Varcount",
"Tumor_Depth", "Normal_Depth", "Driver")]

write.table(svs_palimpsest, file="copy_number_alteration_data_palimpsest_input.txt", quote=F, row.names=F, sep="\t")
