#----------------------------------------------------------------------
#processing_annovar_results.R
#karin isaev
#July 11th, 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/burst2/MUTECT2_selected_VCFs")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 was run on paired mode compaing cns to diagnostic tumour
#now it's time to:
#summarize cns specific mutations
#but first should still filter out false positives (note, these are unfilitered variants)
#see how many appear in multiple comparisons (n=5 total)

#note these vcf files have been normalized and fed through annovar
#for annotations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

vcfs = list.files(pattern=".indels.selected.normalized.vcf.gz") #normalized vcf files

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(vcf){

  pat = unlist(strsplit(vcf, "\\."))[1]

  mutations_T1 =read.vcfR(vcf)
  mutations_T1 = vcfR2tidy(mutations_T1)
  meta = as.data.table(mutations_T1$meta)

  vcf = as.data.table(mutations_T1$fix)
  chr_conv = unique(vcf[,c("ChromKey", "CHROM")])

  gt = as.data.table(mutations_T1$gt)
  gt = merge(gt, chr_conv, by="ChromKey")

  gt$mut_id = paste(gt$CHROM, gt$POS, sep="_")
  vcf$mut_id = paste(vcf$CHROM, vcf$POS, sep="_")

  #1. keep only the ones that passed default mutect2 filter
  vcf = as.data.table(filter(vcf, FILTER=="PASS"))
  print(paste("number of variants that passed filtering=", dim(vcf)[1]))

  #2. combine vcf and gt info
  cols = colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt = merge(gt, vcf, by= cols)

  #3. only tumour remove normal records
  z = which(str_detect(gt$Indiv, "Ctl"))
  gt = gt[-z,]

  #4. remove variants from chromosome X and Y
  gt$CHROM = as.numeric(gt$CHROM)
  gt = as.data.table(filter(gt, (CHROM %in% c(1:22))))

  #5. generate bed file - summary of mutation and coordinates to intersect with cnvkit output
  file = paste(pat, "filtered_mutect2_calls_indels.bed", sep="_")
	gt$patient = pat
	table(gt$CHROM)
  write.table(gt, file, quote=F, row.names=F, sep="\t", col.names=T)
  print("done")
}

llply(vcfs, clean_up_001, .progress="text")
