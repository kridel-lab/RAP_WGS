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
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats", "ggthemes")
lapply(packages, require, character.only = TRUE)

library(RColorBrewer)
display.brewer.all()
display.brewer.pal(9, "Set1")

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
muts = readRDS("final_somatic_mutations_RAP_WGS_20_samples.rds")

#2. Sample summary 
dna = fread("RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
dna = merge(dna, biops, by="barcode")
colnames(dna)[2] = "Indiv"
colnames(dna)[7] = "Tissue_Site"
colnames(dna)[8] = "Specimen_Type" 

dna = as.data.table(filter(dna, Indiv %in% muts$Indiv))
ffpe = as.data.table(matrix(ncol=ncol(dna), nrow=3))
colnames(ffpe) = colnames(dna)
ffpe = as.data.frame(ffpe)
ffpe$Indiv = unique(muts$Indiv[which(!(muts$Indiv %in% dna$Indiv))])
ffpe$Specimen_Type = "FFPE"
ffpe$DNA = "DNA"
ffpe$STUDY_PATIENT_ID = "LY_RAP_0003"
dna = rbind(dna, ffpe)

muts = merge(muts, dna, by="Indiv", all=TRUE)
#increase minimum depth to 60
muts = as.data.table(filter(muts, DP>=60))

#Mutation data from WGS Morin et al 2013
tables2 = fread("TableS2.csv")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

muts$Tissue_Site[is.na(muts$Tissue_Site)] = "FFPE"

#1. -------------------------------------------------------------------

#total mutations 
table(muts$Indiv)

#genes mutations
genes_muts = as.data.table(filter(as.data.table(table(muts$Tissue_Site, muts$hg19.ensemblToGeneName.value)), N >0))
tots = as.data.table(table(muts$hg19.ensemblToGeneName.value))
colnames(tots)[1:2]=c("V2", "total")
genes_muts = merge(genes_muts, tots, by="V2")
genes_muts$fraction = genes_muts$N/genes_muts$total
genes_muts = genes_muts[order(-total, -fraction)]
genes_muts$V2 = factor(genes_muts$V2, levels=unique(genes_muts$V2))
top50 = unique(genes_muts$V2)[1:80]

cols = colorRampPalette(solarized_pal()(3))(17)

g = ggbarplot(filter(genes_muts, V2 %in% top50), x="V2", y="fraction", fill="V1")
ggpar(g, font.tickslab = c(8,"plain", "black"), xtickslab.rt=90, legend ="right")+
  scale_fill_manual(values =cols)  

#genes mutations that occur in the less than 10 areas
t=as.data.table(table(genes_muts$V2))
t=as.data.table(filter(t, N <10))
t=t[order(-N)]
top50 = unique(t$V1[1:100])
g = ggbarplot(filter(genes_muts, V2 %in% top50), x="V2", y="fraction", fill="V1")
ggpar(g, font.tickslab = c(8,"plain", "black"), xtickslab.rt=90, legend ="right")+
  scale_fill_manual(values =cols)  

#genes mutations that are in previously published datatsets
cols = colorRampPalette(solarized_pal()(3))(17)

g = ggbarplot(filter(genes_muts, V2 %in% tables2$Gene), x="V2", y="fraction")
ggpar(g, font.tickslab = c(8,"plain", "black"), xtickslab.rt=90, legend ="right")+
  facet_wrap(~V1, nrow=6, ncol= 3, scales= "free_x") 

g = ggbarplot(filter(genes_muts, V2 %in% tables2$Gene), x="V2", y="fraction", fill="V1")
ggpar(g, font.tickslab = c(8,"plain", "black"), xtickslab.rt=90, legend ="right")+
  scale_fill_manual(values =cols)  

#types of mutations across tissues
genes_muts = as.data.table(filter(as.data.table(table(muts$Tissue_Site, muts$Func.ensGene)), N >0))
tots = as.data.table(table(muts$Func.ensGene))
colnames(tots)[1:2]=c("V2", "total")
genes_muts = merge(genes_muts, tots, by="V2")
genes_muts$fraction = genes_muts$N/genes_muts$total
genes_muts = genes_muts[order(-total, -fraction)]
genes_muts$V2 = factor(genes_muts$V2, levels=unique(genes_muts$V2))

cols = colorRampPalette(solarized_pal()(3))(17)

g = ggbarplot(genes_muts, x="V2", y="fraction")
ggpar(g, font.tickslab = c(12,"plain", "black"), xtickslab.rt=90, legend ="right")+
  facet_wrap(~V1, nrow=6, ncol= 3, scales= "free_x")  


