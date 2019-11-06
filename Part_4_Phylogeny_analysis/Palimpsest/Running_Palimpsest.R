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

#-------------------------------------------------------------------------------------------------
# 2] Extract previously known signatures (COSMIC)
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Signatures_COSMIC");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# define list of colors for visualizing mutational signatures. Selecting default colors
#mycol <- c("darkgreen","deepskyblue4","grey","orangered1","darkred","goldenrod1","deeppink4","royalblue4","darkolivegreen3","purple4");names(mycol) <- liver_signature_names

# Calculating contributions (exposures) of signatures in each sample:
signatures_exp <- deconvolution_fit(vcf=vcf,type = "SNV",input_data = propMutsByCat,threshold = 6,input_signatures = COSMIC_Signatures,sig_cols = mycol,plot = T,resdir = resdir.)

# Plotting the exposures of signatures across the series:
pdf(file.path(resdir.,"signature_content_plot.pdf"),width=10,height=7)
mycol <- readRDS("/cluster/home/kisaev/data/palette_96_cols.rds")
mycol = mycol[1:nrow(COSMIC_Signatures)]
names(mycol) <- rownames(COSMIC_Signatures)
signature_content_plot <- deconvolution_exposure(signatures_exp$sig_nums,signatures_exp$sig_props, sig_cols = mycol)
print(signature_content_plot)
dev.off()

#-------------------------------------------------------------------------------------------------
# 4] Clonality analysis
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Clonality");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Calculate the Cancer Cell Fraction (CCF) of each mutation.
# This step is a bit long... Be patient!
vcf$Tumor_Varcount = as.numeric(vcf$Tumor_Varcount)
vcf$Tumor_Depth = as.numeric(vcf$Tumor_Depth)
vcf$Normal_Depth = as.numeric(vcf$Normal_Depth)

vcf <- cnaCCF_annot(vcf=vcf,annot_data = annotation,cna_data = cnas,CCF_boundary=0.95)

# Generate graphical representations of clonality analysis
cnaCCF_plots(vcf=vcf,resdir=resdir.)

#-------------------------------------------------------------------------------------------------
# 5] Compare mutational signatures between early clonal and late subclonal mutations in each tumor
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Signatures_early_vs_late");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Estimate the contribution of each signature to clonal and subclonal mutations in each tumor
vcf.clonal <- vcf[which(vcf$Clonality=="clonal"),]
propMutsByCat.clonal <- palimpsestInput(vcf = vcf.clonal,type="SNV",sample.col = "Sample", mutcat.col = "mutcat3", proportion = TRUE)
signatures_exp_clonal <- deconvolution_fit(vcf = vcf.clonal,type = "SNV",input_data = propMutsByCat.clonal,threshold = 6,input_signatures = COSMIC_Signatures,sig_cols = mycol,plot = F,resdir = resdir.)

vcf.subclonal <- vcf[which(vcf$Clonality=="subclonal"),]
propMutsByCat.subclonal <- palimpsestInput(vcf = vcf.subclonal,type="SNV",sample.col = "Sample",mutcat.col = "mutcat3",proportion = TRUE)
signatures_exp_subclonal <- deconvolution_fit(vcf = vcf.subclonal,type = "SNV",input_data = propMutsByCat.subclonal,threshold = 6,input_signatures = COSMIC_Signatures,sig_cols = mycol,plot = F,resdir = resdir.)

# Generate tumor-wise comparisons of clonal and subclonal mutations
palimpsest_DissectSigs(vcf=vcf, signatures_exp_clonal = signatures_exp_clonal, signatures_exp_subclonal = signatures_exp_subclonal,sig_cols = mycol,resdir=resdir.)

# Generate across the series comparisons of signature assigned to clonal and subclonal mutations
palimpsest_clonalitySigsCompare(clonsig = signatures_exp_clonal$sig_nums, subsig = signatures_exp_subclonal$sig_nums, msigcol = mycol, resdir = resdir.)

#-------------------------------------------------------------------------------------------------
# 7] Timing Chromosomal Gains
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"ChromosomeDups_timing");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

# Annotate vcf with chromomal gain timings
chrom_dup_time <- chrTime_annot(vcf=vcf,cna_data = cnas, cyto=cyto)
vcf <- chrom_dup_time$vcf; point.mut.time <- chrom_dup_time$point.mut.time; cna_data <- chrom_dup_time$cna_data

# Visualizing timing plots
chrTime_plot(vcf = vcf, point.mut.time = point.mut.time, resdir = resdir.)

#-------------------------------------------------------------------------------------------------
# 8] Visualize the natural history of tumor samples:
#-------------------------------------------------------------------------------------------------
resdir. <- file.path(resdir,"Natural_history");if(!file.exists(resdir.)){dir.create(resdir.)}# Defining the results directory

palimpsest_plotTumorHistories(vcf = vcf, sv.vcf = NULL, cna_data, point.mut.time, clonsig=signatures_exp_clonal$sig_props, subsig=signatures_exp_subclonal$sig_props, msigcol=mycol, resdir=resdir.)





