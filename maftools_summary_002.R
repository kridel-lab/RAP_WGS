#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)
library(ccube)
#> Loading required package: foreach
library(dplyr)
library(maftools)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final")
date = Sys.Date()

source("/cluster/home/kisaev/scripts/annovar_to_maftools.R")

#----------------------------------------------------------------------
#load data
#----------------------------------------------------------------------

#bedtools output files
files = list.files(pattern = "bedtools_output.bed")

#all titan output all patients 
all_titan = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/all_titan_results_pats.txt")
z = which(all_titan$Sample == "Sample")
all_titan = all_titan[-z,]

#purity 
purity = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution.txt")

#cna vcf sample conversion file
conversion = fread("cna_vcf_sample_conversion.csv")

#annovar txt files
anno_files = list.files(pattern="multianno.txt")

#----------------------------------------------------------------------
#get maf files
#----------------------------------------------------------------------

all_mafs = "all_mafs_all_samples.txt"
#f = fread(all_mafs)
#f = as.data.table(filter(f, !(Tumor_Sample_Barcode == "Tumor_Sample_Barcode")))
#write.table(f, file=all_mafs, quote=F, row.names=F, sep="\t") #done

mafs = read.maf(maf = all_mafs)
getSampleSummary(mafs)
getGeneSummary(mafs)
getFields(mafs)
write.mafSummary(maf = mafs, basename = 'RAP_WGS')

pdf("maftools_summary.pdf")

plotmafSummary(maf = mafs, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = mafs, top = 20)
mafs_p = titv(maf = mafs, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = mafs_p)
#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = mafs, top = 25, pvalue = c(0.05, 0.1))
dgi = drugInteractions(maf = mafs, fontSize = 0.75)
OncogenicPathways(maf = mafs)
PlotOncogenicPathways(maf = mafs, pathways = "MYC")
library('NMF')
laml.tnm = trinucleotideMatrix(maf = mafs, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
laml.sign = extractSignatures(mat = laml.tnm, nTry = 6, plotBestFitRes = FALSE)
plotSignatures(laml.sign, title_size = 0.8, )



dev.off()





