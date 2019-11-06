#----------------------------------------------------------------------
#make_pyclone_input.R
#karin isaev
#last updated: August 30th,2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("~/Documents/RAP_analysis")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats", "ggthemes")
lapply(packages, require, character.only = TRUE)
library(cowplot)

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9, "Set1")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarized snvs and cnvs from 21 sequencing folders 
#here, summarize number of mutations/sample/location
#which genes are mutated across all sites which are unique?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data 
#muts = readRDS("final_somatic_mutations_RAP_WGS_20_samples.rds")
#muts = readRDS("2019-08-28_all_soft_filtered_SNVs_overlapping_titan_cna_calls.rds")
#clean up file
#muts = muts[,c(1:11, 13, 38, 47,48, 49:50, 57:61, 65, 69, 71, 73)]
#write.table(muts, paste(date, "mutect2_filtered_variants_wTITAN_CNA_calls.txt", sep="_"), quote=F, row.names=F, sep="\t")
muts = fread("2019-08-29_mutect2_filtered_variants_wTITAN_CNA_calls.txt")

#2. Sample summary 
dna = fread("RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
dna = merge(dna, biops, by="barcode")
colnames(dna)[2] = "Indiv"
colnames(dna)[7] = "Tissue_Site"
colnames(dna)[8] = "Specimen_Type" 
dna$Specimen_Type = "FT"

dna = as.data.table(filter(dna, Indiv %in% muts$vcf_sample))
ffpe = as.data.table(matrix(ncol=ncol(dna), nrow=3))
colnames(ffpe) = colnames(dna)
ffpe = as.data.frame(ffpe)

#ffpe$Indiv = unique(muts$vcf_sample[which(!(muts$vcf_sample %in% dna$Indiv))])
ffpe$Indiv = c("LY_RAP_0003_Dia_FoT_05", "LY_RAP_0003_Dia_FoT_01" ,"LY_RAP_0003_Dia_FoT_03")
ffpe$barcode =c("15:S12966E", "15:S12966A", "15:S12966C")
ffpe$Tissue_Site = c("left_breast", "right_neck_LN", "left_axilla_LN")
ffpe$Specimen_Type = "FFPE"
ffpe$DNA = "DNA"
ffpe$STUDY_PATIENT_ID = "LY_RAP_0003"
dna = rbind(dna, ffpe)

colnames(dna)[2] = "vcf_sample"
muts = merge(muts, dna, by="vcf_sample", all=TRUE)

#Mutation data from WGS Morin et al 2013
tables2 = fread("TableS2.csv")
muts$id = paste(muts$Specimen_Type, muts$Tissue_Site, muts$barcode, sep="_")
muts$Cellular_Prevalence = as.numeric(muts$Cellular_Prevalence)

#for each patient make input file
patients = unique(muts[,c("vcf_sample", "purity")])
pats = patients$vcf_sample
purities = patients$purity

make_input = function(patient){
  pat_dat = as.data.table(filter(muts, vcf_sample == patient))  
  pat_dat$normal_cn = 2
  pat_dat = unique(pat_dat[,c("mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn", "Symbol", "region", "id")])
  pat_dat = as.data.table(filter(pat_dat, major_cn >0))
  write.table(pat_dat, file=paste(patient, "pyclone_input.tsv", sep="_"), quote=F, row.names=F, sep="\t")
  print(table(pat_dat$major_cn))
}

llply(pats, make_input)



