#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("~/Documents/RAP_analysis")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats")
lapply(packages, require, character.only = TRUE)

# get all table files 
#find -L . -name "*cram.pileups.table*" > cram_tables #480

#Data 
all_files = fread("cram_tables", header=F)

#function read files and combine for each patient 
patients = unique(sapply(all_files$V1, function(x){paste(unlist(strsplit(x, "_"))[2:6], collapse="_")}))
patients = unique(sapply(patients, function(x){paste("LY_", x, sep="")}))
write.table(patients, file="all_patient_ids.txt", quote=F, row.names=F, sep="\t")

get_files = function(patient){
	z = which(str_detect(all_files$V1, patient))
	pats_files = all_files[z]
	pats_files$chr = ""
	chrs = unique(sapply(pats_files$V1, function(x){unlist(strsplit(x, "_"))[7]}))
	chrs = unique(sapply(chrs, function(x){unlist(strsplit(x, "\\."))[1]}))
	pats_files$chr = chrs
	pats_files$chr = as.numeric(pats_files$chr)
	pats_files = pats_files[order(chr)]
	pats_files$chr[23:24] = c("X", "Y")
	pats_files$patient = patient
	
	dat = as.data.table(matrix(ncol=6))
	colnames(dat) = c( "contig","position","ref_count","alt_count",   
	"other_alt_count","allele_frequency")
	
	for(i in 1:length(pats_files$V1)){
		f=fread(pats_files$V1[i])
		dat = rbind(dat, f)
	}

	dat = dat[-1,]
	file = paste(patient, "pileups_final.table", sep="_")
	write.table(dat, file, row.names=F, quote=F, sep="\t")

}

#llply(patients, get_files, .progress="text")