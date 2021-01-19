#-------------------------------------------------------------------------------
#prep_SNVs_for_CCF_calculation.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
library(data.table)
library(vcfR)

#-------------------------------------------------------------------------------
#Purpose
#-------------------------------------------------------------------------------

#deBAF in hatchet can run faster is SNP positions are provided
#-L, --snps 	Path to file of SNPs in the format #CHR POS

#-------------------------------------------------------------------------------
#Analysis
#-------------------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS")

vcf_file = "af-only-gnomad.raw.sites.b37.vcf.gz"
#vcf_file = "test_snps.vcf.gz"
#vcf <- read.vcfR(vcf_file, verbose = TRUE)
vcf <- fread(vcf_file)
vcf = vcf[,1:2]
colnames(vcf)=c('#CHR', 'POS')
write.table(vcf, file="deBAF_input_gnomad_snps.txt", quote=F, row.names=F, sep="\t")
