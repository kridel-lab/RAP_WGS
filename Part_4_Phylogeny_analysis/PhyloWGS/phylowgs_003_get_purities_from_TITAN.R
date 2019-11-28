#----------------------------------------------------------------------
#karin isaev
#Nov 1st, 2019
#----------------------------------------------------------------------

date = Sys.Date()

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/TITAN_CNA/results/titan/hmm/optimalClusterSolution_files/titanCNA_ploidy2")

#output from TitanCNA, combine all samples into one dataframe 
#use optimalClusterSolution.txt file to identify optimal cluster for each sample
#use sample to identifier conversion to get actual sample name 

#2]. cna_data: copy number alteration data

#Sample: Sample identifier. Any alphanumeric string.
#CHROM: Chromosome. Between chr1 and chr22 or the chrX or chrY ('chr' prefix required).
#POS_START: Start position of segmented chromosome.
#POS_END: End position of segmented chromosome.
#LogR: LogR information.
#Nmin: Minor allele copy number.
#Nmaj: Major allele copy number.
#ntot: Total copy number of segmented chromosome.
#Ploidy: Tumor ploidy.

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#sample conversion
samples = fread("cna_vcf_sample_conversion.csv")
colnames(samples) = c("barcode", "Sample")

#optimal clusters 
clusters = fread("optimalClusterSolution.txt")
clusters= merge(samples, clusters, by = "barcode")

#titanCNA results 
files = list.files(pattern="seg.txt")
files = files[sapply(clusters$id, function(x){which(str_detect(files, x))})]

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#make CNA input for phyloWGS 
make_cna = function(file){
	
	#read in data files 
	all_cnas = as.data.table(fread(file))
	colnames(all_cnas)[1] = "barcode"
	all_cnas = merge(all_cnas, clusters, by = "barcode")
	all_cnas$CHROM = paste("chr", all_cnas$Chromosome, sep="")
	sample_id = all_cnas$Sample[1]
	print(sample_id)
	purity = unique(all_cnas$purity)
	#cols needed
	#chromosome	start	end	copy_number	minor_cn	major_cn	cellular_prevalence
	all_cnas = all_cnas[,c("CHROM", "Start", "End", "Copy_Number", "MinorCN", "MajorCN", "Cellular_Prevalence")]
	colnames(all_cnas) = c("chromosome", "start", "end", "copy_number",	"minor_cn",	"major_cn", "cellular_prevalence")
	all_cnas$end = all_cnas$end + 1
	print(head(all_cnas))
	write.table(all_cnas, file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/CNAs_titan/", sample_id, "cna_datat_input_phyloWGS_CNV_parser.txt", sep="_"), quote=F, row.names=F, sep="\t")
	print("done")
	return(purity)
}

purities = as.data.table(ldply(llply(files, make_cna)))
write.table(purities, file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/CNAs_titan/", "phylogwgs_CNV_parser_purities.txt", sep="_"), quote=F, row.names=F, sep="\t", col.names=F)

