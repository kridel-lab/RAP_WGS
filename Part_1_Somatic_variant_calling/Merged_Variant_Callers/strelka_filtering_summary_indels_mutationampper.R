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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

vcf_files = fread("strelka_indels.txt", header=F)
colnames(vcf_files)[1] = "file"
vcf_files = as.character(vcf_files$file)

#sample data
samp_dat = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")
patient <- unique(as.character(sapply(vcf_files, function(x){unlist(strsplit(x, "/results/variants/somatic.indels.vcf.gz"))})))
patient <- unique(as.character(sapply(patient, function(x){unlist(strsplit(x, "STRELKA_WORKDIR_"))[2]})))
patient <- unique(as.character(sapply(patient, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})))

get_patient_muts = function(pat){
	z = which(str_detect(vcf_files, pat))
	pats_files = vcf_files[z]
	sample_names <- unique(as.character(sapply(pats_files, function(x){unlist(strsplit(x, "/results/variants/somatic.indels.vcf.gz"))})))
	sample_names <- unique(as.character(sapply(sample_names, function(x){unlist(strsplit(x, "STRELKA_WORKDIR_"))[2]})))

	#sample_type = as.character(sapply(sample_names, function(x){unlist(strsplit(x, "_"))[5]}))
	indel_grl = read_vcfs_as_granges(pats_files, sample_names, ref_genome, type = "indel")
	indel_grl <- get_indel_context(indel_grl, ref_genome)
	indel_counts <- count_indel_contexts(indel_grl)
	head(indel_counts)

	pdf(paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/", pat, "_strelka_mutation_mapper_INDLS_summary.pdf", sep=""), width=15)
	p1 = plot_indel_contexts(indel_counts, condensed = TRUE)
	p2 = plot_main_indel_contexts(indel_counts)
	print(p1)
	print(p2)
	signatures_indel = get_known_signatures(muttype = "indel")
	signatures_artifacts = get_known_signatures(incl_poss_artifacts = TRUE)
	dim(signatures_artifacts)
	signatures_signal = get_known_signatures(source = "SIGNAL")
	dev.off()

}#end function

llply(patient, get_patient_muts, .progress="text")
