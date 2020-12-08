#----------------------------------------------------------------------
#karin isaev
#Nov 1st, 2019
#----------------------------------------------------------------------

#make sure R 4.0.0 is loaded
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

#-------------------------------------------------------------------------------------------------
# [1] set up data and directories
#-------------------------------------------------------------------------------------------------

#get patient ID
args = commandArgs(trailingOnly = TRUE) #patient ID
index = args[1]
print(index)
#index="LY_RAP_0001"
samp=index

#We first define the directory containing input data and the directory where result files should
#be exported.

# directory containing input files
datadir <- "input"
# directory to export results
resdir <- paste("Results", index, sep="_");if(!file.exists(resdir)) dir.create(resdir)
#load input files
annotation = fread(paste(datadir, "/annotation_data_palimpsest_input.txt", sep=""))
cnas = fread(paste(datadir, "/copy_number_alteration_data_palimpsest_input.txt", sep=""))
mut_data = readRDS(paste(datadir, "/SNV_input_Palimpsest.rds", sep=""))

#filter annotation data to only include it for patient being analyzed
annotation = annotation[which(str_detect(annotation$Sample, samp)),]

#filter cna data to only include it for patient being analyzed
cnas = cnas[which(str_detect(cnas$Sample, samp)),]

#filter mutation data to only include it for patient being analyzed
mut_data = mut_data[which(str_detect(mut_data$Sample, samp)),]

#-------------------------------------------------------------------------------------------------
# [2] set up genome
#-------------------------------------------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
ref_genome <- BSgenome.Hsapiens.UCSC.hg19
load("/cluster/home/kisaev/Palimpsest/data/ensgene_hg19.RData")
load("/cluster/home/kisaev/Palimpsest/data/cytoband_hg19.RData")

#-------------------------------------------------------------------------------------------------
# [3] prepare mutation data (SNVs only)
#-------------------------------------------------------------------------------------------------

mut_data = as.data.frame(mut_data)
#vcf <- preprocessInput_snv(input_data = mut_data,ensgene=ensgene,reference_genome = ref_genome)
vcf <- annotate_VCF(vcf = mut_data, ref_genome = BSgenome.Hsapiens.UCSC.hg19)
SBS_input <- palimpsest_input(vcf = vcf, Type = "SBS")

#-------------------------------------------------------------------------------------------------
# [4] Estimating the exposures of mutational signatures
#-------------------------------------------------------------------------------------------------

# Calculate and plot the exposure of the signatures across the series
SBS_col <- signature_colour_generator(rownames(SBS_cosmic))

SBS_signatures_exp = deconvolution_fit(input_matrices = SBS_input,
	input_signatures = SBS_cosmic, signature_colours = SBS_col, resdir = resdir)

deconvolution_exposure(signature_colours = SBS_col,
	signature_contribution = SBS_signatures_exp)

#-------------------------------------------------------------------------------------------------
# [5] Assigning the most likely signature at the origin of each mutation
#-------------------------------------------------------------------------------------------------

vcf <- signature_origins(input = vcf, Type = "SBS",
	input_signatures = SBS_cosmic, signature_contribution = SBS_signatures_exp)

vcf$Driver[vcf$Driver == ""] = NA
vcf.cod <- vcf[(!is.na(vcf$Driver) & vcf$Type=="SNV"),]
vcf.cod <- signature_origins(input = vcf.cod, Type = "SBS", input_signatures = SBS_cosmic,
                             signature_contribution = SBS_signatures_exp)

# Estimate and represent the cumulative contribution of signatures to each driver gene
drivers <- unique(vcf.cod$Driver)

matprob <- matrix(nrow=length(drivers),ncol=length(SBS_cosmic),dimnames=list(drivers, SBS_cosmic))
sig.cols <- paste0(rownames(SBS_cosmic),".prob")#grep("prob",colnames(vcf.cod))
for(i in 1:nrow(matprob)){
		g <- rownames(matprob)[i]
		ind <- which(vcf.cod$gene_name==g)
		matprob[i,] <- apply(vcf.cod[ind,sig.cols],2,sum,na.rm=T)
}
barplot(t(matprob),col = SBS_col, border = SBS_col, las=2)
legend("top",names(SBS_col)[names(SBS_col) %in% rownames(SBS_cosmic)],fill=SBS_col,ncol=5,
	cex=0.75,bty ="n",inset = c(0,-0.3),xpd = T)

#-------------------------------------------------------------------------------------------------
# [6] Clonality analysis
#-------------------------------------------------------------------------------------------------

#Copy number alterations and Cancer cell fraction (CCF)
vcf_cna <- cnaCCF_annot(vcf = vcf, annot_data = annotation, cna_data = cnas, CCF_boundary = 0.95)

#clonality plots
cnaCCF_plots(vcf= vcf_cna, resdir = resdir)

#temporal evolution of mutational signatures
# Estimate the contribution of each signature to clonal and
# subclonal mutations in each tumour
vcf.clonal <- vcf_cna[which(vcf_cna$Clonality=="clonal"),]

SBS_input_clonal <- palimpsest_input(vcf = vcf.clonal, Type = "SBS")

sig_exp_clonal <- deconvolution_fit(input_matrices = SBS_input_clonal, input_signatures = SBS_cosmic,
	resdir =  resdir, save_signatures_exp = F)

vcf.subclonal <- vcf_cna[which(vcf_cna$Clonality=="subclonal"),]

SBS_input_subclonal <- palimpsest_input(vcf = vcf.subclonal,Type = "SBS")

sig_exp_subclonal <- deconvolution_fit(input_matrices = SBS_input_subclonal,
	input_signatures = SBS_cosmic, resdir = resdir, save_signatures_exp = F)

# Generate per tumour comparisons of clonal and subclonal mutations
palimpsest_DissectSigs(vcf=vcf_cna, signatures_exp_clonal = sig_exp_clonal,
	signatures_exp_subclonal = sig_exp_subclonal, sig_cols = SBS_col, resdir=resdir)

# Generate across the series comparisons of signature assigned to clonal and
#subclonal mutations
palimpsest_clonalitySigsCompare(clonsig = sig_exp_clonal$sig_nums,
	subsig = sig_exp_subclonal$sig_nums,msigcol = SBS_col, resdir = resdir)

# Timing chromosomal gains

# Annotate vcf with chromomal gain timings
chrom_dup_time <- chrTime_annot(vcf = vcf_cna, cna_data = cnas,
	cyto = cytoband_hg19)
vcf_cna <- chrom_dup_time$vcf
point.mut.time <- chrom_dup_time$point.mut.time
cnas <- chrom_dup_time$cna_data
# Visualising timing plots
chrTime_plot(vcf = vcf_cna, point.mut.time = point.mut.time,
	resdir = resdir,cyto = cytoband_hg19)

#-------------------------------------------------------------------------------------------------
# 7] Visualise the natural history of tumour samples:
#-------------------------------------------------------------------------------------------------
resdir <- file.path(resdir,"Natural_history/");if(!file.exists(resdir)){dir.create(resdir)}

palimpsest_plotTumorHistories(vcf = vcf_cna, cna_data =  cnas,
	                              point.mut.time = point.mut.time,
	                              clonsig = sig_exp_clonal$sig_props,
	                              subsig = sig_exp_subclonal$sig_props,
	                              msigcol = SBS_col, resdir = resdir)
