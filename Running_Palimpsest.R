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
	"ggrepel", "stringr", "maftools", "Palimpsest")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Palimpsest")

#We first define the directory containing input data and the directory where result files should
#be exported.

#------------------------------------------------------------------
# directory containing input files 
datadir <- "input"
# directory to export results
resdir <- "Results";if(!file.exists(resdir)) dir.create(resdir)
#load input files 
annotation = fread(paste(datadir, "/annotation_data_palimpsest_input.txt", sep=""))            
cnas = fread(paste(datadir, "/copy_number_alteration_data_palimpsest_input.txt", sep=""))            
mut_data = readRDS(paste(datadir, "/SNV_input_Palimpsest.rds", sep=""))
mut_data$Gene_Name = NA
mut_data$Driver = NA
#------------------------------------------------------------------

#set up genome 

#------------------------------------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19)
ref_genome <- BSgenome.Hsapiens.UCSC.hg19
load("/cluster/home/kisaev/Palimpsest/data/ensgene_hg19.RData")
load("/cluster/home/kisaev/Palimpsest/data/cytoband_hg19.RData")
#------------------------------------------------------------------

#Preprocessing and annotating input data
#load("/cluster/home/kisaev/Palimpsest/RUNNING_PALIMPSEST_EXAMPLE/LiC1162/mut_data.RData")
#load("/cluster/home/kisaev/Palimpsest/RUNNING_PALIMPSEST_EXAMPLE/LiC1162/cna_data.RData")
#load("/cluster/home/kisaev/Palimpsest/RUNNING_PALIMPSEST_EXAMPLE/LiC1162/annot_data.RData")
#load("/cluster/home/kisaev/Palimpsest/RUNNING_PALIMPSEST_EXAMPLE/LiC1162/sv_data.RData")

#-------------------------------------------------------------------------------------------------
# 1] De novo mutational signature analysis
#-------------------------------------------------------------------------------------------------

mut_data = as.data.frame(mut_data)
vcf <- preprocessInput_snv(input_data = mut_data,ensgene=ensgene,reference_genome = ref_genome)

propMutsByCat <- palimpsestInput(vcf = vcf,type = "SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = TRUE)
denovo_signatures <- deconvolution_nmf(input_data = propMutsByCat,type = "SNV",range_of_sigs = 2:10,nrun = 20,method = "brunet",resdir = resdir)

# Compare with existing signatures from COSMIC database:
pdf(file.path(resdir, "Cosine_Similarity.pdf"))
cosine_similarities <- deconvolution_compare(denovo_signatures,COSMIC_Signatures)# missing lsa
dev.off()

# Define color codes for signatures
library(RColorBrewer);qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
mycol <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mycol<- mycol[sample.int(length(mycol),nrow(denovo_signatures))];names(mycol) <- rownames(denovo_signatures)

# Calculating contributions (exposures) of signatures in each sample and generate tumor-wise graphical outputs:
signatures_exp <- deconvolution_fit(vcf=vcf,type = "SNV",input_data = propMutsByCat,threshold = 5,input_signatures = denovo_signatures,sig_cols = mycol,plot = T,resdir = resdir)

# Plotting the exposures of signatures across the series:
pdf(file.path(resdir,"signature_content_plot.pdf"),width=10,height=7)
signature_content_plot <- deconvolution_exposure(signatures_exp$sig_nums,signatures_exp$sig_props,sig_cols = mycol)
print(signature_content_plot)
dev.off()


