#----------------------------------------------------------------------
#make_pyclone_input.R
#karin isaev
#last updated: August 30th,2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "ggthemes")
lapply(packages, require, character.only = TRUE)
library(RColorBrewer)
library(openxlsx)
library(plotly)

#2. Sample summary
dna = fread("RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
dna = merge(dna, biops, by="barcode")
colnames(dna)[2] = "Indiv"
colnames(dna)[7] = "Tissue_Site"
colnames(dna)[8] = "Specimen_Type"
dna$Specimen_Type = "FT"

dna = as.data.table(filter(dna, Indiv %in% muts$Indiv))
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
dna$id = paste(dna$Specimen_Type, dna$Tissue_Site, dna$barcode, sep="_")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarized snvs and cnvs from 21 sequencing folders
#here, summarize number of mutations/sample/location
#which genes are mutated across all sites which are unique?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. SNV/CNA input data for pyclone

muts = fread(list.files(pattern="PYCLONE_INPUT_MUTS.txt")[length(list.files(pattern="PYCLONE_INPUT_MUTS.txt"))]) #get most recent mutation file
pats = unique(muts$id)
t = as.data.table(table(muts$mut_id))
t=t[order(V1)]

#make sure mutations ordered in the same way in each sample specific file

make_input = function(patient){
  pat_dat = as.data.table(filter(muts, id == patient))
  pat_dat$normal_cn = 2
  pat_dat$Indiv = NULL
  pat_dat = merge(pat_dat, dna, by = c("id"))

  pat_dat = unique(pat_dat[,c("mut_id", "Ref_counts", "alt_counts", "normal_cn", "Nmin", "Nmaj", "hg19.ensemblToGeneName.value", "Func.ensGene", "Indiv")])
  print(table(pat_dat$Nmaj))

  z = which(pat_dat$mut_id %in% t$V1)
  pat_dat = pat_dat[z,]

  colnames(pat_dat) = c("mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn", "gene_name", "region", "id")
  print(tail(pat_dat))
  write.table(pat_dat, file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone/", pat_dat$id[1], "pyclone_input.tsv", sep="_"), quote=F, row.names=F, sep="\t")
}

llply(pats, make_input)
