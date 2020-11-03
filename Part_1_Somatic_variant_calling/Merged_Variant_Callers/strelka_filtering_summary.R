#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

date = Sys.Date()
print(date)

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)
library(VariantAnnotation)
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#library(BSgenome.Hsapiens.UCSC.hg19)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/strelka_filtered")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

paired = list.files(pattern=".vcf.gz.normalized.vcf.gz")

#sample data
samp_dat = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

clean_up_001 = function(paired_vcf){

  mutations_T1 =read.vcfR(paired_vcf)
  mutations_T1 = vcfR2tidy(mutations_T1)
  meta = as.data.table(mutations_T1$meta)

  vcf = as.data.table(mutations_T1$fix)
  chr_conv = unique(vcf[,c("ChromKey", "CHROM")])

  gt = as.data.table(mutations_T1$gt)
  gt = merge(gt, chr_conv, by="ChromKey")

  gt$mut_id = paste(gt$CHROM, gt$POS, sep="_")
  vcf$mut_id = paste(vcf$CHROM, vcf$POS, sep="_")

  #1. keep only the ones that passed default hard filter
  vcf = as.data.table(filter(vcf, FILTER=="PASS"))
  print(paste("number of variants that passed filtering=", dim(vcf)[1]))

  #2. filter by coverage
  #try 60X
  vcf = as.data.table(filter(vcf, DP >=60))
  print(paste("number of variants that passed coverage=", dim(vcf)[1]))

  #3. combine vcf and gt info
  cols = colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt = merge(gt, vcf, by= cols)

	#4. only tumour remove normal records
	gt = as.data.table(filter(gt, Indiv == "TUMOR"))

	gt$CHROM = as.numeric(gt$CHROM)
  gt = as.data.table(filter(gt, (CHROM %in% c(1:22))))

  #5. remove variants from chromosome X and Y
  gt = as.data.table(filter(gt, !(CHROM %in% c("X", "Y"))))
  print(paste("number of variants that passed X Y=", dim(gt)[1]))

	#6. save file
  gt$pat = unlist(strsplit(paired_vcf, "_strelka"))[1]
  print("done")

	return(gt)
}

all_sampls = as.data.table(ldply(llply(paired, clean_up_001, .progress="text")))
all_sampls = as.data.table(filter(all_sampls, SomaticEVS >=17))
all_sampls = all_sampls[order(pat)]
all_sampls$pat = factor(all_sampls$pat, levels=unique(all_sampls$pat))

#brief barplot
pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/strelka_mutations_summary.pdf")

# Basic barplot
ggplot(all_sampls, aes(x=factor(pat)))+
  geom_bar(stat="count", width=0.7, fill="steelblue")+
  theme_minimal()+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
