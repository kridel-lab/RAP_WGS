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
library("readxl")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Sequenza")

#[2]. cna_data: copy number alteration data

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

all_segs = fread("all_segments_files_sequenza.txt", header=F)$V1 #27 samples results

read_files = function(seq_file){
	segs = fread(seq_file)
	sample = unlist(strsplit(seq_file, "/"))[1]
	segs$Sample = sample
	return(segs)
}

all_segs_dt = as.data.table(ldply(llply(all_segs, read_files)))

samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")
colnames(samps)[4] ="Sample"

all_segs_dt = merge(all_segs_dt, samps, by="Sample")
print(head(all_segs_dt))

#purity and ploidy info
purities = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Purity_Ploidy_Results_Hatchet_Sequenza.xlsx")) %>%
	filter(Tool == "Sequenza") %>% select(Patient, Sample, Purity, Ploidy)
purities$Sample = paste(purities$Patient, purities$Sample, sep="_")
all_segs_dt = merge(all_segs_dt, purities, by="Sample")

all_segs_dt$CHROM = paste("chr", all_segs_dt$chromosome, sep="")
print(head(all_segs_dt))

#save for filtering SNVs
all_cnas = all_segs_dt[,c("Sample", "Tissue_Site", "CHROM", "start.pos", "end.pos",
"B", "A", "Ploidy", "Purity", "CNt", "depth.ratio")]
#save full dataset
colnames(all_cnas)[4:5] = c("Start", "End")
colnames(all_cnas)[6:7] = c("Nmin", "Nmaj")
colnames(all_cnas)[10] = c("ntot")
print(head(all_cnas))

saveRDS(all_cnas, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")

all_cnas_palimpsest = all_cnas[,c("Sample", "CHROM", "Start", "End",
"depth.ratio", "Nmin", "Nmaj",
"ntot", "Ploidy")]
colnames(all_cnas_palimpsest) = c("Sample", "CHROM", "POS_START",
"POS_END", "LogR", "Nmin", "Nmaj", "ntot", "Ploidy")
write.table(all_cnas_palimpsest, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Palimpsest/input/copy_number_alteration_data_via_Sequenza_palimpsest_input.txt", quote=F, row.names=F, sep="\t")
