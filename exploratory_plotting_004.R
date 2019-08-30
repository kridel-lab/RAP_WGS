#----------------------------------------------------------------------
#exploratory_plotting_004.R
#karin isaev
#last updated: August 30th,2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("/Users/kisaev/treeomics/src/input/RAP_vcfs")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats", "ggthemes")
lapply(packages, require, character.only = TRUE)
library(cowplot)
library(MutationalPatterns)
library(BSgenome)
library(RColorBrewer)
library(gridExtra)
display.brewer.all()
display.brewer.pal(9, "Set1")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#here we will look at mutation signatures that my be present across the samples

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Sample summary 
dna = fread("~/Documents/RAP_analysis/RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("~/Documents/RAP_analysis/RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
dna = merge(dna, biops, by="barcode")
colnames(dna)[2] = "Indiv"
colnames(dna)[7] = "Tissue_Site"
colnames(dna)[8] = "Specimen_Type" 
dna$Specimen_Type = "FT"

dna = as.data.table(filter(dna, Indiv %in% muts$vcf_sample))
ffpe = as.data.table(matrix(ncol=ncol(dna), nrow=3))
colnames(ffpe) = colnames(dna)
ffpe = as.data.frame(ffpe)

#ffpe$Indiv = unique(muts$vcf_sample[which(!(muts$vcf_sample %in% dna$Indiv))])
ffpe$Indiv = c("LY_RAP_0003_Dia_FoT_05", "LY_RAP_0003_Dia_FoT_01" ,"LY_RAP_0003_Dia_FoT_03")
ffpe$barcode =c("15:S12966E", "15:S12966A", "15:S12966C")
ffpe$Tissue_Site = c("left_breast", "right_neck_LN", "left_axilla_LN")
ffpe$Specimen_Type = "FFPE"
ffpe$DNA = "DNA"
ffpe$STUDY_PATIENT_ID = "LY_RAP_0003"
dna = rbind(dna, ffpe)
dna$id = paste(dna$Specimen_Type, dna$Tissue_Site, dna$barcode, sep="_")

#2. set up Mutational Pattern
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

#Locate the VCF files of the example data:
vcf_files <- list.files(pattern = ".vcf", full.names = TRUE)

#Define corresponding sample names for the VCF files:
sample_names <- as.character(sapply(vcf_files, function(x){paste(unlist(strsplit(x, "_"))[1:6], collapse="_")}))
sample_names <- as.character(sapply(sample_names, function(x){unlist(strsplit(x, "/"))[2]}))
saples = as.data.table(sample_names) ; colnames(saples)="Indiv"
saple = merge(saples, dna, by="Indiv")

#Load the VCF files into a GRangesList:
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
summary(vcfs)
tissue <- saple$id

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#-------------------------
#Mutation characteristics#
#-------------------------

#Base substitution types
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
type_occurrences

#-------------------------
#Mutation spectrum       #
#-------------------------

p1 <- plot_spectrum(type_occurrences) # Error bars indicate standard deviation
#over all samples. The total number of mutations is indicated.

#Plot the mutation spectrum with distinction between C>T at CpG sites and other sites:
p2 <- plot_spectrum(type_occurrences, CT = TRUE)
#Plot spectrum without legend:
p3 <- plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)

grid.arrange(p1, p2, p3, ncol=3, widths=c(3,3,1.75))

#each person sepperatley
p4 <- plot_spectrum(type_occurrences, by = tissue, CT = TRUE, legend = TRUE)

#-------------------------
#96 mutational profile   #
#-------------------------

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
head(mut_mat)
plot_96_profile(mut_mat[,c(18:20)],condensed = TRUE)

#--------------------------------------------------
#De novo mutational signature extraction using NMFs#
#--------------------------------------------------

mut_mat <- mut_mat + 0.0001
#Use the NMF package to generate an estimate rank plot:
library("NMF")
estimate <- nmf(mut_mat, rank=2:5, method="brunet", nrun=10, seed=123456)
plot(estimate)

nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 50)

#Assign signature names:
colnames(nmf_res$signatures) <- c("Signature A", "Signature B")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B")

#Plot the 96-profile of the signatures:
plot_96_profile(nmf_res$signatures, condensed = TRUE)

#Visualize the contribution of the signatures in a barplot:
pc1 <- plot_contribution(nmf_res$contribution, nmf_res$signature,
                             mode = "relative")
#Visualize the contribution of the signatures in absolute number of mutations:
pc2 <- plot_contribution(nmf_res$contribution, nmf_res$signature,
                             mode = "absolute")
#Combine the two plots:
grid.arrange(pc1, pc2)

plot_contribution(nmf_res$contribution, nmf_res$signature,
                    mode = "absolute", coord_flip = TRUE)

pch1 <- plot_contribution_heatmap(nmf_res$contribution,
                                    sig_order = c("Signature B", "Signature A"))
#Plot signature contribution as a heatmap without sample clustering:
pch2 <- plot_contribution_heatmap(nmf_res$contribution, cluster_samples=FALSE)
#Combine the plots into one figure:
#grid.arrange(pch1, pch2, ncol = 2, widths = c(2,1.6))
pch1 <- plot_contribution_heatmap(nmf_res$contribution,
                                  sig_order = c("Signature A", "Signature B"))
pch1

#--------------------------------------------------
#COSMIC mutational signatures                     #
#--------------------------------------------------

#Download mutational signatures from the COSMIC website:
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
                    "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)

# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)

# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]

# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type

# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])

#Plot mutational profile of the first two COSMIC signatures:
plot_96_profile(cancer_signatures[,c(17,28)], condensed = TRUE, ymax = 0.3)

hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
plot(hclust_cosmic)

#Calculate pairwise cosine similarity between mutational profiles and COSMIC signatures:
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
# Plot heatmap with specified signature order
plot_cosine_heatmap(cos_sim_samples_signatures,
                        col_order = cosmic_order,
                        cluster_rows = TRUE)

#Fit mutation matrix to the COSMIC mutational signatures:
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)

#Plot relative contribution of the cancer signatures in each sample as a heatmap with sample
#clustering:
plot_contribution_heatmap(fit_res$contribution,
                              cluster_samples = TRUE,
                              method = "complete")

#Genomic distribution
#A rainfall plot visualizes mutation types and intermutation distance. Rainfall plots can be
#used to visualize the distribution of mutations along the genome or a subset of chromosomes.
#The y-axis corresponds to the distance of a mutation with the previous mutation and is log10
#transformed. Drop-downs from the plots indicate clusters or “hotspots” of mutations.
#Make rainfall plot of sample 1 over all autosomal chromosomes

# Define autosomal chromosomes
chromosomes <- seqnames(get(ref_genome))[1:22]
# Make a rainfall plot
palette <- c("pink", "orange", "blue", "lightblue", "green", "red")
plot_rainfall(vcfs[[1]], title = names(vcfs[1]), chromosomes = chromosomes, cex = 1.5, ylim = 1e+09, colors=palette)

library(biomaRt)
regulatory <- useEnsembl(biomart="regulation",
                         dataset="hsapiens_regulatory_feature",
                          GRCh = 37)

## Download the regulatory CTCF binding sites and convert them to
## a GRanges object.

CTCF <- getBM(attributes = c('chromosome_name',
  'chromosome_start',
  'chromosome_end',
  'feature_type_name'),
  filters = "regulatory_feature_type_name",
  values = "CTCF Binding Site",
  mart = regulatory)

CTCF_g <- reduce(GRanges(CTCF$chromosome_name,
  IRanges(CTCF$chromosome_start,
  CTCF$chromosome_end)))
  
CTCF_g <- readRDS(system.file("states/CTCF_g_data.rds",
                                package="MutationalPatterns"))

# Download the promoter regions and convert them to a GRanges object.

promoter = getBM(attributes = c('chromosome_name', 'chromosome_start',
  'chromosome_end', 'feature_type_name'),
  filters = "regulatory_feature_type_name",
  values = "Promoter",
  mart = regulatory)

promoter_g = reduce(GRanges(promoter$chromosome_name,
  IRanges(promoter$chromosome_start,
  promoter$chromosome_end)))
  
promoter_g <- readRDS(system.file("states/promoter_g_data.rds",
                                      package="MutationalPatterns"))

## Download the promoter flanking regions and convert them to a GRanges object.

flanking = getBM(attributes = c('chromosome_name',
  'chromosome_start',
  'chromosome_end',
  'feature_type_name'),
  filters = "regulatory_feature_type_name",
  values = "Promoter Flanking Region",
  mart = regulatory)

flanking_g = reduce(GRanges(
  flanking$chromosome_name,
  IRanges(flanking$chromosome_start,
  flanking$chromosome_end)))
  
flanking_g <- readRDS(system.file("states/promoter_flanking_g_data.rds",
                                      package="MutationalPatterns"))


#Combine all genomic regions (GRanges objects) in a named list:
regions <- GRangesList(promoter_g, flanking_g, CTCF_g)
names(regions) <- c("Promoter", "Promoter flanking", "CTCF")
#Use the same chromosome naming convention consistently:
seqlevelsStyle(regions) <- "UCSC"


## Get the filename with surveyed/callable regions
surveyed_file <- system.file("extdata/callableloci-sample.bed",
                                 package = "MutationalPatterns")
## Import the file using rtracklayer and use the UCSC naming standard
library(rtracklayer)
surveyed <- import(surveyed_file)
seqlevelsStyle(surveyed) <- "UCSC"
## For this example we use the same surveyed file for each sample.
surveyed_list <- rep(list(surveyed), 20)

## Calculate the number of observed and expected number of mutations in
## each genomic regions for each sample.
distr <- genomic_distribution(vcfs, surveyed_list, regions)

## Perform the enrichment/depletion test by tissue type.
distr_test <- enrichment_depletion_test(distr, by = tissue)
head(distr_test)
plot_enrichment_depletion(distr_test)
