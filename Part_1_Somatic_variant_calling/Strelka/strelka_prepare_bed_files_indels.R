#----------------------------------------------------------------------

#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/strelka_filtered")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

vcf_files = list.files(pattern="indels.vcf.gz.normalized.vcf.gz")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(vcf){

  pat = unlist(strsplit(vcf, "_strelka_s"))[1]

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
	gt = as.data.table(filter(gt, Indiv == "TUMOR"))

  #4. remove variants from chromosome X and Y
  gt$CHROM = as.numeric(gt$CHROM)
  gt = as.data.table(filter(gt, (CHROM %in% c(1:22))))

  #5. generate bed file - summary of mutation and coordinates to intersect with cnvkit output
  file = paste(pat, "filtered_strelka_calls_indels.bed", sep="_")
	gt$patient = pat
	table(gt$CHROM)
  write.table(gt, file, quote=F, row.names=F, sep="\t", col.names=T)
  print("done")
}

llply(vcf_files, clean_up_001, .progress="text")
