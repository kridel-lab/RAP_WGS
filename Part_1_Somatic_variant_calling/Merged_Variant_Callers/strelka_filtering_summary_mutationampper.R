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
library(NMF)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(biomaRt)
library(GenomicRanges)

#regulatory <- useEnsembl(biomart="regulation",
#                          dataset="hsapiens_regulatory_feature",
#                          GRCh = 37)

## Download the regulatory CTCF binding sites and convert them to
## a GRanges object.
#CTCF <- getBM(attributes = c('chromosome_name',
#                             'chromosome_start',
#                             'chromosome_end',
#                             'feature_type_name'),
#              filters = "regulatory_feature_type_name",
#              values = "CTCF Binding Site",
#              mart = regulatory)

#CTCF_g <- reduce(GRanges(CTCF$chromosome_name,
#                 IRanges(CTCF$chromosome_start,
#                 CTCF$chromosome_end)))

#saveRDS(CTCF_g, file="CTCF_g_biomart.rds") #"/cluster/home/kisaev"
CTCF_g = readRDS("/cluster/home/kisaev/CTCF_g_biomart.rds")

## Download the promoter regions and convert them to a GRanges object.
#promoter = getBM(attributes = c('chromosome_name', 'chromosome_start',
#                                 'chromosome_end', 'feature_type_name'),
#                  filters = "regulatory_feature_type_name",
#                  values = "Promoter",
#                  mart = regulatory)

#promoter_g = reduce(GRanges(promoter$chromosome_name,
#                     IRanges(promoter$chromosome_start,
#                             promoter$chromosome_end)))
#saveRDS(promoter_g, file="promoter_g_biomart.rds") #"/cluster/home/kisaev"
promoter_g = readRDS("/cluster/home/kisaev/promoter_g_biomart.rds")

## Download the promoter flanking regions and convert them to a GRanges object.
#flanking = getBM(attributes = c('chromosome_name',
#                                 'chromosome_start',
#                                 'chromosome_end',
#                                 'feature_type_name'),
#                  filters = "regulatory_feature_type_name",
#                  values = "Promoter Flanking Region",
#                  mart = regulatory)

#flanking_g = reduce(GRanges(
#                        flanking$chromosome_name,
#                        IRanges(flanking$chromosome_start,
#                        flanking$chromosome_end)))

#saveRDS(flanking_g, file="flanking_g_biomart.rds") #"/cluster/home/kisaev"
flanking_g = readRDS("/cluster/home/kisaev/flanking_g_biomart.rds")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/strelka_filtered")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

vcf_files = list.files(pattern=".vcf.gz.normalized.vcf.gz")

#sample data
samp_dat = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")
patient <- unique(as.character(sapply(vcf_files, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse= "_")})))

get_patient_muts = function(pat){
	z = which(str_detect(vcf_files, pat))
	pats_files = vcf_files[z]
	sample_names = as.character(sapply(pats_files, function(x){unlist(strsplit(x, "_strelka"))[1]}))
	sample_type = as.character(sapply(sample_names, function(x){unlist(strsplit(x, "_"))[5]}))

	grl = read_vcfs_as_granges(pats_files, sample_names, ref_genome)
	type_occurrences <- mut_type_occurrences(grl, ref_genome)
	type_occurrences

	pdf(paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/", pat, "_strelka_mutation_mapper_summary.pdf", sep=""), width=15)

		p1 <- plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE)
		print(p1)

		#96 mutational profile
		#First you should make a 96 trinucleodide mutation count matrix.

		mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
		head(mut_mat)
		p2 = plot_96_profile(mut_mat,condensed=T)
		print(p2)

		#Cosmic signatures
		#To do this you first need to read in some already existing signatures.
		#Here we will use signatures from COSMIC (v3.1) (Alexandrov et al. 2020).
		signatures = get_known_signatures()
		fit_res <- fit_to_signatures(mut_mat, signatures)
		strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
		fit_res_strict <- strict_refit$fit_res
		p3 = plot_contribution(fit_res_strict$contribution,
		  coord_flip = TRUE,
		  mode = "relative")

		merged_signatures <- merge_signatures(signatures, cos_sim_cutoff = 0.8)

		strict_refit <- fit_to_signatures_strict(mut_mat, merged_signatures, max_delta = 0.004)
		fit_res_strict <- strict_refit$fit_res
		p4 = plot_contribution(fit_res_strict$contribution,
			coord_flip = TRUE,
			mode = "relative")

		print(p3)
		print(p4)

		contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
	  merged_signatures,
	  n_boots = 100,
	  method = "strict")

		p5 = plot_bootstrapped_contribution(contri_boots,
		                               mode = "relative",
		                               plot_type = "dotplot")

		cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, merged_signatures)
		p6 = plot_cosine_heatmap(cos_sim_samples_signatures,
		                    cluster_rows = TRUE, cluster_cols = TRUE, method="ward.D2")

		print(p5)
		print(p6)

		#strand bias analysis
		#For the mutations within genes it can be determined whether the
		#mutation is on the transcribed or non-transcribed strand, which
		#can be used to evaluate the involvement of transcription-coupled repair.
		genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
		strand <- mut_strand(grl[[1]], genes_hg19)
		head(strand, 10)

		mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19)
		p7 = plot_192_profile(mut_mat_s)

		strand_counts <- strand_occurrences(mut_mat_s, by = sample_names)
		head(strand_counts)

		#Next, you can use these counts to perform a Poisson test for strand asymmetry.
		#Multiple testing correction is also performed.
		strand_bias <- strand_bias_test(strand_counts)
		head(strand_bias)

		p8 <- plot_strand(strand_counts, mode = "relative")

		#Plot the effect size (log2(untranscribed/transcribed) of the strand bias.
		#Asteriks indicate significant strand bias. Here we use p-values to plot asterisks.
		#By default fdr is used.
		p9 <- plot_strand_bias(strand_bias, sig_type = "p")

		#genomic distribution
		chromosomes <- seqnames(get(ref_genome))[1:22]

		print(p7)
		print(p8)
		print(p9)

		regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
		names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
		seqlevelsStyle(regions) <- "UCSC"

	dev.off()

}#end function

llply(patient, get_patient_muts, .progress="text")
