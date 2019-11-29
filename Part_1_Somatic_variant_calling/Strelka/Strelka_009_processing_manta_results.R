#----------------------------------------------------------------------
#Strelka_009_processing_manta_results.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
              "ggrepel", "stringr", "maftools", "VariantAnnotation", "ggpubr")
lapply(packages, require, character.only = TRUE)
library(GenomicRanges)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.lincRNAsTranscripts)
library(biomaRt)
library(openxlsx)
library(plotly)

date = Sys.Date()

print(date)
#args = commandArgs(trailingOnly = TRUE) #patient ID 
#index = args[1] 
#print(index) 

setwd("/Users/kisaev/Documents/RAP_ANALYSIS")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

svs = fread(list.files(pattern="all_SVs_samples.txt")[1])
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#distribution of types of structural variation
pats_svs = as.data.table(table(svs$SVTYPE, svs$pat)) ; pats_svs = pats_svs[order(-N)]

ggbarplot(pats_svs, x = "V2", y="N", fill="V1") +theme_bw() + 
  rotate_x_text(75) + ylab("Types of SVs") + xlab("Sample")

#genes most impacted by SVs 
genes_svs = as.data.table(table(svs$hgnc_symbol, svs$pat, svs$SVTYPE)) ; genes_svs = genes_svs[order(-N)]
genes_svs = as.data.table(filter(genes_svs, N >0, !(V1 == "")))
theme_update(text = element_text(size=7))

p = ggplot(filter(genes_svs, V1 %in% morin$Gene), aes(V2, V1)) +
  geom_tile(aes(fill = V3), colour = "grey50") + 
  rotate_x_text(90)

ggplotly(p)


