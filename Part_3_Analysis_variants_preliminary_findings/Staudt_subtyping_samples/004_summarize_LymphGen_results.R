#----------------------------------------------------------------------
#004_summarize_LymphGen_result.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)

date = Sys.Date()
options(scipen=999) #no scientific notation

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust",
              "data.table", "openxlsx",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr",
              "ggExtra", "broom", "ggthemes")

lapply(packages, require, character.only = TRUE)
library("readxl")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare sample annotation file to be able to run the tool by Staudt
#https://www.sciencedirect.com/science/article/pii/S1535610820301550

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#ABC GCB subtypes from Dave paper
genes_class = read.xlsx("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure2_data_ABC_GCB_gene_Dave.xlsx")

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))

#samples id names
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

read_only = filter(read_only, STUDY_PATIENT_ID == "LY_RAP_0003")
read_only$driver = ""
read_only$driver[read_only$symbol %in% reddy$Gene] = "yes"
read_only$driver[!(read_only$symbol %in% reddy$Gene)] = "no"

read_only$subtype = ""
read_only$subtype[read_only$symbol %in% genes_class$Genes[genes_class$Subtype.association == "ABC"]] = "ABC"
read_only$subtype[read_only$symbol %in% genes_class$Genes[genes_class$Subtype.association == "GCB"]] = "GCB"

table(read_only$driver, read_only$Func.ensGene)
table(read_only$driver, read_only$ExonicFunc.ensGene)
table(filter(read_only, biotype == "protein_coding")$driver, filter(read_only, biotype == "protein_coding")$Func.ensGene)

#summarize ABC GCB gene mutations
subtypes = as.data.table(filter(read_only, !(subtype=="")))
subtypes = subtypes[order(subtype, symbol)]
subtypes$symbol = factor(subtypes$symbol, levels=unique(subtypes$symbol))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/summary_COO_genes_across_samples.pdf", width=6, height=6)
p = ggplot(subtypes, aes(symbol, Sample)) +
  geom_tile(aes(fill = subtype), colour = "grey50") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
   xlab("Gene") + ylab("Sample name") #+
  #scale_fill_gradient2(low = "black", mid = "grey", midpoint = 0.5,
  #                     high = "red", na.value="transparent")
p
dev.off()

#only exonic and non synonymous mutations
subtypes = as.data.table(filter(read_only, !(subtype==""), ExonicFunc.ensGene=="nonsynonymous_SNV"))
subtypes = subtypes[order(subtype, symbol)]
subtypes$symbol = factor(subtypes$symbol, levels=unique(subtypes$symbol))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/summary_COO_genes_across_samples_nonsynonymous_only.pdf", width=6, height=6)
p = ggplot(subtypes, aes(symbol, Sample)) +
  geom_tile(aes(fill = subtype), colour = "grey50") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
   xlab("Gene") + ylab("Sample name") #+
  #scale_fill_gradient2(low = "black", mid = "grey", midpoint = 0.5,
  #                     high = "red", na.value="transparent")
p
dev.off()

#results files
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/LymphGen")

#result
res = fread("july2021_Result.txt")

#compare
comp = fread("july2021_Compare.txt")

colnames(res)[1] = "Indiv"
res = merge(res, samps, by="Indiv")

#summarize confience for each subtype
res_sum = melt(res, id.vars = "id", measure.vars=c("Confidence.BN2",
"Confidence.EZB", "Confidence.MCD", "Confidence.N1",
"Confidence.ST2", "Confidence.A53"))
pdf("001_LymphGen_summary.pdf", width=6, height=6)
# Basic barplot
p<-ggplot(data=res_sum, aes(x=id, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+theme_minimal()+ggtitle("Confidence for each subtype (LymphGen)")+
  xlab("Sample") + ylab("Confidence")+
    coord_flip()
print(p)
dev.off()

#summarize number of features associated with each subtype
res_sum = melt(res, id.vars = "id", measure.vars=c("BN2.Feature.Count",
"EZB.Feature.Count", "MCD.Feature.Count", "N1.Feature.Count",
"ST2.Feature.Count", "A53.Feature.Count"))
pdf("002_LymphGen_summary.pdf", width=6, height=6)
# Basic barplot
p<-ggplot(data=res_sum, aes(x=id, y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge())+theme_minimal()+ggtitle("Number of features in each subtype (LymphGen)")+
  xlab("Sample") + ylab("Number of features")+
    coord_flip()
print(p)
dev.off()
