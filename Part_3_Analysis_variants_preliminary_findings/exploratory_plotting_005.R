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
library(openxlsx)
library(plotly)

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
muts = fread("2019-11-19_READ_ONLY_ALL_MERGED_MUTS.txt") 
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. -------------------------------------------------------------------

#total mutations 
sum_muts = as.data.table(table(muts$id)) ; sum_muts = sum_muts[order(-N)]
print(sum_muts)
ggbarplot(sum_muts, x = "V1", y="N", fill="grey") +theme_bw() + 
  rotate_x_text(65) + ylab("Number of SNVs") + xlab("Sample")
ggsave(paste(date, "num_muts_per_samples_final.pdf", sep="_"))  

#how many founder mutatiotns 
founds = filter(as.data.table(table(muts$mut_id, muts$Indiv)), N >0)
founds = as.data.table(filter(as.data.table(table(founds$V1)), N ==20)) #8268 unique variants present in everyone, 18902/52302 = 36%
founds_muts = as.data.table(filter(muts, mut_id %in% founds$V1))

#how many unique mutations per patient? 
unique = filter(as.data.table(table(muts$mut_id, muts$id)), N >0)
unique_muts = as.data.table(filter(as.data.table(table(unique$V1)), N ==1)) #9133 unique variants present, 14800/52302 = 28%
unique_muts = merge(unique_muts, unique, by = c("V1", "N"))
unique_muts_sum = as.data.table(table(unique_muts$V2)) ;  unique_muts_sum = unique_muts_sum[order(-N)]
print(unique_muts_sum)
ggbarplot(unique_muts_sum, x = "V1", y="N", fill="grey") +theme_bw() + 
  rotate_x_text(65) + ylab("Number of unique SNVs") + xlab("Sample")
ggsave(paste(date, "num_unique_muts_per_samples_final.pdf", sep="_"))  

#where are they?
unique_muts = (as.data.table(filter(muts, mut_id %in% unique_muts$V1)))
unique_muts = unique(unique_muts[,c("mut_id", "Func.ensGene", "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene", "cosmic68", "hg19.ensemblToGeneName.value")])
genes = as.data.table(table(unique_muts$hg19.ensemblToGeneName.value))

#Functional of genes
unique_muts_plot = unique(muts[,c("mut_id", "Func.ensGene")])
sum_muts = as.data.table(table(unique_muts_plot$Func.ensGene)) ; sum_muts = sum_muts[order(-N)]
print(sum_muts)
sum_muts$N = sum_muts$N
ggbarplot(sum_muts, x = "V1", y="N", fill="grey") +theme_bw() + 
  rotate_x_text(65) + ylab("Number of SNVs") + xlab("Gene region")
ggsave(paste(date, "num_muts_per_genomic_region.pdf", sep="_"))  

#Functional of exons
sum_muts = as.data.table(table(muts$ExonicFunc.ensGene)) ; sum_muts = sum_muts[order(-N)]
print(sum_muts)
sum_muts$N = sum_muts$N / 1000000
ggbarplot(sum_muts, x = "V1", y="N", fill="grey", label = TRUE) +theme_bw() + 
  rotate_x_text(65) + ylab("Number of SNVs/ Megabase") + xlab("Gene region")
ggsave(paste(date, "num_muts_per_exonic_region.pdf", sep="_"))  


#remove founder and unique variants since they are not informative for this task 
muts_matrix = as.data.frame(dcast(muts, mut_id ~ id, value.var = "gt_AF"))
muts_matrix = subset(muts_matrix, !(mut_id %in% founds$mut_id))
muts_matrix = subset(muts_matrix, !(mut_id %in% unique_muts$mut_id))

rownames(muts_matrix) = muts_matrix$mut_id
muts_matrix$mut_id = NULL

#replace all NAs with zeros for now
muts_matrix[is.na(muts_matrix)] = 0

require(vegan)
muts_matrix = t(muts_matrix)
dist.mat<-vegdist(muts_matrix,method="jaccard", na.rm = TRUE)
clust.res<-hclust(dist.mat)
plot(clust.res)
ggsave(paste(date, "hclust_prelim_tree.pdf", sep="_"))  

#3. -------------------------------------------------------------------
#which genes most mutated ? 
#in founders? unique samples? across the board?

#1. in founds 
genes_founds = as.data.table(table(founds$hg19.ensemblToGeneName.value)) ; genes_founds = genes_founds[order(-N)]

#2. in unique 
genes_unique = as.data.table(table(unique_muts$hg19.ensemblToGeneName.value)) ; genes_unique = genes_unique[order(-N)]

#3. across the board 
genes_all = as.data.table(table(muts$hg19.ensemblToGeneName.value)) ; genes_all = genes_all[order(-N)]

#4. -------------------------------------------------------------------
#just Y-RNAs? what's going on with them?
all_Y_RNA = as.data.table(filter(muts, hg19.ensemblToGeneName.value == "Y_RNA"))
ensgs = as.data.table(table(all_Y_RNA$Gene.ensGene))
ensgs = ensgs[order(-N)]
y_rna_muts = as.data.table(table(all_Y_RNA$mut_id))
y_rna_muts = y_rna_muts[order(-N)]

#5. diagnostic versus autopsy 
sample_gene = as.data.table(table(muts$Gene.ensGene, muts$Specimen_Type))
sample_gene = sample_gene[order(-N)]
sample_gene = as.data.table(filter(sample_gene, N >= 1))
genes_both = as.data.table(table(sample_gene$V1))
genes_unique = as.data.table(filter(genes_both, N == 1))
genes_both = as.data.table(filter(genes_both, N >=2))
conv = unique(muts[,c("Specimen_Type", "Gene.ensGene")])
colnames(genes_unique)[1] = "Gene.ensGene"
genes_unique = merge(genes_unique, conv, by = "Gene.ensGene")

