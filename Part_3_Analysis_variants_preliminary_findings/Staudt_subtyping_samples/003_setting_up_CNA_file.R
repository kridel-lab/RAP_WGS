#----------------------------------------------------------------------
#001_setting_up_sample_annotation_file.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/TITAN_CNA/results/titan/hmm/optimalClusterSolution_files/titanCNA_ploidy2")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust",
              "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr",
              "ggExtra", "broom", "ggthemes", "GenomicRanges")

lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare sample annotation file to be able to run the tool by Staudt
#https://www.sciencedirect.com/science/article/pii/S1535610820301550

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#sample annotation
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#cnas

#output from TitanCNA, combine all samples into one dataframe
#use optimalClusterSolution.txt file to identify optimal cluster for each sample
#use sample to identifier conversion to get actual sample name

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

#read in data files
all_cnas = as.data.table(ldply(llply(files, function(x){fread(x)})))
colnames(all_cnas)[1] = "barcode"
all_cnas = merge(all_cnas, clusters, by = "barcode")
all_cnas$CHROM = paste("chr", all_cnas$Chromosome, sep="")

#entrez id conversion
#gene annotations
genes = unique(fread("/cluster/home/kisaev/data/annotables_grch37.txt"))
genes = as.data.table(filter(genes, biotype == "protein_coding"))
genes = as.data.table(filter(genes, !(is.na(entrez))))

#intersect CNA regions with genes
#prepare dataframe for conversion to granges
genes_gr = genes[,c("chr", "start", "end", "strand", "entrez", "symbol")]
genes_gr$strand = as.character(genes_gr$strand)
genes_gr$strand = "*"
genes_gr = makeGRangesFromDataFrame(genes_gr, keep.extra.columns=TRUE)

all_cnas_gr = all_cnas[,c("Chromosome", "Start", "End", "numClust",
"Sample", "TITAN_call")]
all_cnas_gr$numClust = as.character(all_cnas_gr$numClust)
all_cnas_gr$numClust = "*"
all_cnas_gr = makeGRangesFromDataFrame(all_cnas_gr, keep.extra.columns=TRUE)

hits <- findOverlaps(all_cnas_gr, genes_gr, ignore.strand=TRUE)
hits_overlap = cbind(all_cnas[queryHits(hits),], genes[subjectHits(hits),])
hits_overlap = unique(hits_overlap[,c("Sample", "entrez", "TITAN_call")])
colnames(hits_overlap)=c("Sample.ID", "ENTREZ.ID", "Type")

#format file needs to be in
#column1 - Sample.ID
#column2 - ENTREZ.ID
#column3 - Type (GAIN (single copy increase), AMP (> 2 or more),
#HETLOSS , HOMDEL)

hits_overlap$Type[hits_overlap$Type == "ASCNA"] = "AMP"
hits_overlap$Type[hits_overlap$Type == "BCNA"] = "AMP"
hits_overlap$Type[hits_overlap$Type == "UBCNA"] = "AMP"
hits_overlap$Type[hits_overlap$Type == "DLOH"] = "HETLOSS"
hits_overlap$Type[hits_overlap$Type == "HOMD"] = "HOMDEL"

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/LymphGen")

#make CNA list of genes
gene_list = as.data.table(hits_overlap$ENTREZ.ID)
colnames(gene_list) = "ENTREZ.ID"
gene_list = unique(gene_list)
write.table(gene_list, paste(date, "LymphGen_Sample_Flat_CNAs_gene_list.txt", sep="_"), quote=F, row.names=F, sep="\t")

hits_overlap = as.data.table(filter(hits_overlap, Type %in% c("AMP",
"HETLOSS", "GAIN", "HOMDEL")))
write.table(hits_overlap, paste(date, "LymphGen_Sample_Flat_CNAs.txt", sep="_"), quote=F, row.names=F, sep="\t")
