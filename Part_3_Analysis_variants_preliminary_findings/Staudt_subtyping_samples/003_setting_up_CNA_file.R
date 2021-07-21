#----------------------------------------------------------------------
#001_setting_up_sample_annotation_file.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS")

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

cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")
cnas$Patient = sapply(cnas$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
cnas = filter(cnas, Patient == "LY_RAP_0003")

#sample annotation
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")
colnames(samps)[2] = "Sample"

#cnas

#output from Sequenza, combine all samples into one dataframe

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

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
genes_gr$chr = paste("chr", genes_gr$chr, sep="")
genes_gr = makeGRangesFromDataFrame(genes_gr, keep.extra.columns=TRUE)

all_cnas_gr = cnas[,c("CHROM", "Start", "End", "Nmin",
"Sample", "Nmaj", "ntot")]
all_cnas_gr$Nmin = as.character(all_cnas_gr$Nmin)
all_cnas_gr$Nmin = "*"
all_cnas_gr = makeGRangesFromDataFrame(all_cnas_gr, keep.extra.columns=TRUE)

hits <- findOverlaps(all_cnas_gr, genes_gr, ignore.strand=TRUE)
hits_overlap = cbind(cnas[queryHits(hits),], genes[subjectHits(hits),])
upload_onedrive = merge(hits_overlap, samps, by="Sample")
write.table(upload_onedrive, file=paste(date,
  "all_CNAs_protein_coding_samples.txt", sep="_"), quote=F, row.names=F, sep="\t")

hits_overlap = unique(hits_overlap[,c("Sample", "entrez", "ntot", "Nmaj")])
colnames(hits_overlap)=c("Sample.ID", "ENTREZ.ID", "ntot", "Nmaj")

#format file needs to be in
#column1 - Sample.ID
#column2 - ENTREZ.ID
#column3 - Type (GAIN (single copy increase), AMP (> 2 or more),
#HETLOSS , HOMDEL)

hits_overlap$Type = ""

#LOH
z = which((hits_overlap$ntot ==  hits_overlap$Nmaj) & (hits_overlap$ntot == 2))
hits_overlap$Type[z] = "HETLOSS"

#Neutral
z = which((hits_overlap$ntot !=  hits_overlap$Nmaj) & (hits_overlap$ntot == 2))
hits_overlap$Type[z] = "Neutral"

#Deletions
z = which((hits_overlap$ntot ==  hits_overlap$Nmaj) & (hits_overlap$ntot == 0))
hits_overlap$Type[z] = "HOMDEL"

z = which((hits_overlap$ntot ==  hits_overlap$Nmaj) & (hits_overlap$ntot == 1))
hits_overlap$Type[z] = "HETLOSS"

#Amplications
z = which((hits_overlap$Type == "") & (hits_overlap$ntot ==3))
hits_overlap$Type[z] = "GAIN"

z = which((hits_overlap$Type == "") & (hits_overlap$ntot > 3))
hits_overlap$Type[z] = "AMP"

z = which(is.na(hits_overlap$ntot))
if(!(length(z) == 0)){
  hits_overlap=hits_overlap[-z,]
}

hits_overlap = filter(hits_overlap, !(Type == "Neutral"))
hits_overlap$ntot = NULL
hits_overlap$Nmaj = NULL

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/LymphGen")

#make CNA list of genes
gene_list = as.data.table(hits_overlap$ENTREZ.ID)
colnames(gene_list) = "ENTREZ.ID"
gene_list = unique(gene_list)
write.table(gene_list, paste(date, "LymphGen_Sample_Flat_CNAs_gene_list.txt", sep="_"), quote=F, row.names=F, sep="\t")

hits_overlap = as.data.table(filter(hits_overlap, Type %in% c("AMP",
"HETLOSS", "GAIN", "HOMDEL")))
write.table(hits_overlap, paste(date, "LymphGen_Sample_Flat_CNAs.txt", sep="_"), quote=F, row.names=F, sep="\t")
