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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/vcfs_annovar_annotated")
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
conversion = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/cna_vcf_sample_conversion.csv")

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2:3] = c("Gene.ensGene", "Symbol")

#----------------------------------------------------------------------
#make input for run_ccube
#----------------------------------------------------------------------

get_input = function(file){
	#mutation_id, minor_cn, major_cn, total_cn, total_counts, var_counts, ref_counts, purity 
	f = fread(file)
	colnames(f)[1:13] = c("chr_cna", "start_cna", "end_cna", "chr_snv", "start_snv", "end_snv", 
		"ref", "alt", "region", "gene_name", "nothing", "type_mut", "transcript")

	colnames(f)[45] = "ad_vals"

	f = f %>% separate(ad_vals, c("GT", "AD", "AF"), sep=":") %>% separate(AD, c("ref_counts", "var_counts"))
	#keep only AF > 0.1
	f = as.data.table(filter(f, AF >= 0.1, AF <0.9))
	colnames(f)[32] = "population_var"
	z = which(str_detect(f$population_var, "rs"))
	f = f[-z,]
	colnames(f)[37] = "position_snv"
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

	f = as.data.table(filter(f, total_counts >=60))

	colnames(f)[colnames(f)=="MinorCN"] = "minor_cn"
	colnames(f)[colnames(f)=="MajorCN"] = "major_cn"
	colnames(f)[colnames(f)=="Copy_Number"] = "total_cn"
	f$purity = filter(purity, barcode==sample)$purity

	#for phylowgs save varaints (smaller list)
	file_name = paste(vcf_sample, "variant_for_phylowgs_input.txt", sep="_")

	#CHROM, POS, ALT, and REF. 
	input = f[,c("chr_snv", "position_snv",  "ref", "alt")]
	write.table(input, file=file_name, sep="\t", row.names=F, col.names=F, quote=F)

	f = f[,c("mutation_id", "minor_cn", "major_cn", "total_cn", "total_counts", "var_counts", "ref_counts", "purity")]
	#for calder keep only varaints in copy neutral regions 
	f$minor_cn = as.numeric(f$minor_cn)
	f$major_cn = as.numeric(f$major_cn)
	f$total_cn = as.numeric(f$total_cn)
	f$total_counts = as.numeric(f$total_counts)
	f$var_counts = as.numeric(f$var_counts)
	f$ref_counts = as.numeric(f$ref_counts)
	f$purity = as.numeric(f$purity)
	f$cna_sample = sample
	f$vcf_sample = vcf_sample

	#run Ccube pipeline
	#numOfClusterPool = 1:4
	#numOfRepeat = 2
	#results <- RunCcubePipeline(ssm = f, 
    #                        numOfClusterPool = numOfClusterPool, 
    #                        numOfRepeat = numOfRepeat,
    #                        runAnalysis = T, 
    #                        runQC = T)

	#print(MakeCcubeStdPlot(ssm = results$ssm, res = results$res, printPlot = F))
	#ccfs = as.data.table(results$ssm)
	#ccfs$cna_sample = sample
	#ccfs$vcf_sample = vcf_sample
	print("done")
	return(f)
}

all_ccfs = llply(files, get_input, .progress="text")
all_ccfs = as.data.table(ldply(all_ccfs))
all_ccfs$gene_name = unlist(all_ccfs$gene_name)
get_gene = function(gene){
	if(length(unlist(strsplit(gene, ";"))) > 1){
		gene_keep = unlist(strsplit(gene, ";"))[1]
		gene = gene_keep
	}
	return(gene)
}

all_ccfs$Gene.ensGene = unlist(sapply(all_ccfs$gene_name, get_gene))

all_ccfs = merge(all_ccfs, genes, by="Gene.ensGene")
saveRDS(all_ccfs, file=paste(date, "all_soft_filtered_SNVs_overlapping_titan_cna_calls.rds", sep="_"))


















#variants='LY_RAP_0003_Aut_FzT_07_variant_for_phylowgs_input.txt'
#vcf_in='/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/LY_RAP_0003_Aut_FzT_07_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz'
#vcf_out='variants_of_interest.vcf'

#bcftools view -O v -R "$variants" "$vcf_in" \
# | grep -Ef <(awk 'BEGIN{FS=OFS="\t";print "#"};{print "^"$1,$2,"[^\t]+",$3,$4"\t"}' "$variants") \
# > "$vcf_out"


















