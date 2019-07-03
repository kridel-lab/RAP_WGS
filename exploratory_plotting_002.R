#----------------------------------------------------------------------
#exploratory_plotting_002.R
#karin isaev
#last updated: June 24th 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("~/Documents/RAP_analysis")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarized snvs and cnvs from 21 sequencing folders 
#here, summarize number of mutations/sample/location
#which genes are mutated across all sites which are unique?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data 
muts = fread("2019-06-24_summary_dnaseq_mutations_RAP.txt")

#2. Summary CNV data 
cnvs = fread("2019-06-24_summary_dnaseq_CNVs_RAP.txt")

#3. Sample summary 
dna = fread("RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
dna = merge(dna, biops, by="barcode")
colnames(dna)[2] = "sample"
colnames(dna)[7] = "Tissue_Site"
colnames(dna)[8] = "Specimen_Type" 

dna = as.data.table(filter(dna, sample %in% muts$sample))

muts = merge(muts, dna, by="sample", all=TRUE)
cnvs = merge(cnvs, dna, by="sample", all=TRUE)

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

muts$Tissue_Site[muts$sample == "LY_RAP_0003_Ctl_FzG_01"] = "control"
cnvs$Tissue_Site[cnvs$sample == "LY_RAP_0003_Ctl_FzG_01"] = "control"

z = which(str_detect(muts$dbsnp_wind, "rs"))
muts = muts[-z,]
muts$Tissue_Site[is.na(muts$Tissue_Site)] = "FFPE"

#1. -------------------------------------------------------------------

#total mutations 
table(muts$sample)

#genes mutations
genes_muts = as.data.table(filter(as.data.table(table(muts$Tissue_Site, muts$gene_symbol)), N >0))
tots = as.data.table(table(muts$Tissue_Site))
colnames(tots)[2]="total"
genes_muts = merge(genes_muts, tots, by="V1")
genes_muts$fraction = genes_muts$N/genes_muts$total
genes_muts = genes_muts[order(-fraction)]
genes_muts$V2 = factor(genes_muts$V2, levels=unique(genes_muts$V2))
top50 = unique(genes_muts$V2)[1:30]

g = ggbarplot(filter(genes_muts, V2 %in% top50), x="V2", y="fraction", fill="V1")
ggpar(g, font.tickslab = c(10,"plain", "black"), xtickslab.rt=65, legend ="right")


