#----------------------------------------------------------------------
#summarize pipeline stats as provided by sickkids
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools", "biomaRt", "openxlsx")
lapply(packages, require, character.only = TRUE)

date = Sys.Date()

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#these are processed variant files coming from TCAG

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

folders = list.files()[which(str_detect(list.files(), "files"))]

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. get dnaseq pipeline stats
read_dat = function(folder){
	files = list.files(folder)
	stats = files[which(str_detect(files, "pipeline_stats.tsv"))]
	stats = fread(paste(folder, stats, sep="/"))
	type = paste(unlist(strsplit(folder, "_"))[4:5], collapse="_") #get Aut_FzT status for example..
	stats$type = type
	return(stats)
}

all_dnaseq_stats = llply(folders, read_dat)
all_dnaseq_stats = ldply(all_dnaseq_stats)

file = paste(date, "summary_dnaseq_pipeline_stats_RAP_WGS_samples.txt", sep="_")
write.table(all_dnaseq_stats, file=file, quote=F, row.names=F, sep="\t")
file = paste(date, "summary_dnaseq_pipeline_stats_RAP_WGS_samples.xlsx", sep="_")
write.xlsx(all_dnaseq_stats, file)

#2. summarize mutations for each sample

get_muts = function(folder){
	files = list.files(paste(folder, "annot", sep="/"))
	anno = files[which(str_detect(files, "SUBSET"))]
	anno = fread(paste(folder, "annot", anno, sep="/"))
	type = paste(unlist(strsplit(folder, "_"))[4:5], collapse="_") #get Aut_FzT status for example..
	anno$type = type
	sample = paste(unlist(strsplit(folder, "_"))[1:6], collapse="_")
	anno$sample = sample

	#remove dbSNPs
	z = which(str_detect(anno$dbsnp_region, "rs"))
	anno = anno[-z,]

	#keep only those that passed all filters GATK
	anno = as.data.table(filter(anno, FILTER == "PASS"))

	#filter by overall GNOMAD allele frequency
	#anno = as.data.table(filter(anno, gnomAD_genome_ALL < 0.01))
	z = which(str_detect(colnames(anno), sample))
	colnames(anno)[z] = sapply(colnames(anno)[z], function(x){unlist(strsplit(x, ":"))[2]})
	colnames(anno)[z]
	colnames(anno)[1] = "chrom"

	return(anno)
}

all_muts_sum = llply(folders, get_muts, .progress="text")
all_muts_sum = as.data.table(ldply(all_muts_sum))
file = paste(date, "summary_dnaseq_mutations_RAP.txt", sep="_")
write.table(all_muts_sum, file=file, quote=F, row.names=F, sep="\t")

#3. summarize CNVs for each sample

get_cnvs = function(folder){
	files = list.files(paste(folder, "cnv", sep="/"))
	anno = files[which(str_detect(files, "tagged"))]
	anno = fread(paste(folder, "cnv", anno, sep="/"))
	type = paste(unlist(strsplit(folder, "_"))[4:5], collapse="_") #get Aut_FzT status for example..
	anno$type = type
	sample = paste(unlist(strsplit(folder, "_"))[1:6], collapse="_")
	anno$sample = sample

	#keep only those that passed all filters GATK
	anno = as.data.table(filter(anno, FILTER == "PASS"))

	#there is a column called imprecise but not sure how it is defined so not going to filter by it for now
	return(anno)
}

all_cnvs_sum = llply(folders, get_cnvs, .progress="text")
all_cnvs_sum = as.data.table(ldply(all_cnvs_sum))
file = paste(date, "summary_dnaseq_CNVs_RAP.txt", sep="_")
write.table(all_cnvs_sum, file=file, quote=F, row.names=F, sep="\t")

#summarize genes and CNVs across frozen, FFPE and control
#add site
#RSTUDIO locally
