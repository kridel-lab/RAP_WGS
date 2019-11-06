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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final")
date = Sys.Date()

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

#----------------------------------------------------------------------
#make input for run_ccube
#----------------------------------------------------------------------

get_input = function(file){
	#mutation_id, minor_cn, major_cn, total_cn, total_counts, var_counts, ref_counts, purity 
	f = fread(file)
	colnames(f) = c("chr_cna", "start_cna", "end_cna", "chr_snv", "position_snv", 
		"nothing", "ref", "alt", "nothingtwo", "filter", "details", "chars", "ad_vals")
	f = f %>% separate(ad_vals, c("GT", "AD", "AF"), sep=":") %>% separate(AD, c("ref_counts", "var_counts"))
	#keep only AF > 0.1
	f = as.data.table(filter(f, AF >= 0.1, AF <0.9))
	f$mutation_id = paste(f$chr_snv, f$position_snv, sep="_")
	#get major minor status 
	sample = paste(unlist(strsplit(file, "_"))[1:6], collapse="_")
	sample = unlist(strsplit(sample, "\\{"))[2]
	vcf_sample = sample
	#get cna sample
	sample = conversion$CNA_sample[conversion$VCF_sample==sample]
	sam_cna_dat = as.data.table(filter(all_titan, Sample == sample))
	colnames(sam_cna_dat)[2:4] = c("chr_cna", "start_cna", "end_cna")
	sam_cna_dat$chr_cna = as.numeric(sam_cna_dat$chr_cna)
	sam_cna_dat$start_cna = as.numeric(sam_cna_dat$start_cna)
	sam_cna_dat$end_cna = as.numeric(sam_cna_dat$end_cna)

	f = merge(f, sam_cna_dat, by=c("chr_cna", "start_cna", "end_cna"))
	f$total_counts = as.numeric(f$ref_counts) + as.numeric(f$var_counts)
	f = as.data.table(filter(f, total_counts >=100))
	colnames(f)[colnames(f)=="MinorCN"] = "minor_cn"
	colnames(f)[colnames(f)=="MajorCN"] = "major_cn"
	colnames(f)[colnames(f)=="Copy_Number"] = "total_cn"
	f$purity = filter(purity, barcode==sample)$purity

	f = f[,c("mutation_id", "minor_cn", "major_cn", "total_cn", "total_counts", "var_counts", "ref_counts", "purity")]
	f$minor_cn = as.numeric(f$minor_cn)
	f$major_cn = as.numeric(f$major_cn)
	f$total_cn = as.numeric(f$total_cn)
	f$total_counts = as.numeric(f$total_counts)
	f$var_counts = as.numeric(f$var_counts)
	f$ref_counts = as.numeric(f$ref_counts)
	f$purity = as.numeric(f$purity)

	#run Ccube pipeline
	numOfClusterPool = 1:6
	numOfRepeat = 1
	results <- RunCcubePipeline(ssm = f, 
                            numOfClusterPool = numOfClusterPool, 
                            numOfRepeat = numOfRepeat,
                            runAnalysis = T, 
                            runQC = T)
	print(MakeCcubeStdPlot(ssm = results$ssm, res = results$res, printPlot = F))
	ccfs = as.data.table(results$ssm)
	ccfs$cna_sample = sample
	ccfs$vcf_sample = vcf_sample

	return(ccfs)
}

pdf("all_Ccube_plots_clones.pdf")
all_ccfs = llply(files, get_input, .progress="text")
dev.off()

#complete input for calder 






