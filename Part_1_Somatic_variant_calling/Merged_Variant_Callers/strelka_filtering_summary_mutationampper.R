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
		p2=plot_96_profile(mut_mat)

		mut_mat_ext_context <- mut_matrix(grl, ref_genome, extension = 2)
		p3 = plot_profile_heatmap(mut_mat_ext_context, by = sample_names)
		p4 = plot_river(mut_mat_ext_context)

		print(p2)
		print(p3)
		print(p4)

		#Cosmic signatures
		#To do this you first need to read in some already existing signatures.
		#Here we will use signatures from COSMIC (v3.1) (Alexandrov et al. 2020).
		signatures = get_known_signatures()
		fit_res <- fit_to_signatures(mut_mat, signatures)

		p5 = plot_contribution(fit_res$contribution,
  	coord_flip = FALSE,
  	mode = "relative")

		p6 = plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed,
                               y_intercept = 0.95)

		strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
		fit_res_strict <- strict_refit$fit_res
		p7 = plot_contribution(fit_res_strict$contribution,
		  coord_flip = FALSE,
		  mode = "relative")

		merged_signatures <- merge_signatures(signatures, cos_sim_cutoff = 0.8)

		strict_refit <- fit_to_signatures_strict(mut_mat, merged_signatures, max_delta = 0.004)
		fit_res_strict <- strict_refit$fit_res
		p8 = plot_contribution(fit_res_strict$contribution,
			coord_flip = FALSE,
			mode = "relative")

		print(p5)
		print(p6)
		print(p7)
		print(p8)

		contri_boots <- fit_to_signatures_bootstrapped(mut_mat,
	  merged_signatures,
	  n_boots = 100,
	  method = "strict")

		p9 = plot_bootstrapped_contribution(contri_boots)

		p10 = plot_bootstrapped_contribution(contri_boots,
		                               mode = "relative",
		                               plot_type = "dotplot")

		cos_sim_samples_signatures <- cos_sim_matrix(mut_mat, signatures)
		p11 = plot_cosine_heatmap(cos_sim_samples_signatures,
		                    cluster_rows = TRUE, cluster_cols = TRUE)

		cos_sim_samples <- cos_sim_matrix(mut_mat, mut_mat)
		p12 = plot_cosine_heatmap(cos_sim_samples, cluster_rows = TRUE, cluster_cols = TRUE)

		print(p9)
		print(p10)
		print(p11)
		print(p12)

		#strand bias analysis
		#For the mutations within genes it can be determined whether the
		#mutation is on the transcribed or non-transcribed strand, which
		#can be used to evaluate the involvement of transcription-coupled repair.
		genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
		strand <- mut_strand(grl[[1]], genes_hg19)
		head(strand, 10)

		mut_mat_s <- mut_matrix_stranded(grl, ref_genome, genes_hg19)
		p13 = plot_192_profile(mut_mat_s)

		strand_counts <- strand_occurrences(mut_mat_s, by = sample_names)
		head(strand_counts)

		#Next, you can use these counts to perform a Poisson test for strand asymmetry.
		#Multiple testing correction is also performed.
		strand_bias <- strand_bias_test(strand_counts)
		head(strand_bias)

		p14 <- plot_strand(strand_counts, mode = "relative")

		#Plot the effect size (log2(untranscribed/transcribed) of the strand bias.
		#Asteriks indicate significant strand bias. Here we use p-values to plot asterisks.
		#By default fdr is used.
		p15 <- plot_strand_bias(strand_bias, sig_type = "p")

		#genomic distribution
		chromosomes <- seqnames(get(ref_genome))[1:22]
		for(i in 1:length(sample_names)){
			print(i)
			# Make a rainfall plot
			p = plot_rainfall(grl[[i]],
  		title = names(grl[i]),
  		chromosomes = chromosomes, cex = 1, ylim = 1e+09)
			print(p)
		}

	dev.off()

}#end function

llply(patient, get_patient_muts, .progress="text")
