#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

#GCB vs ABC mutations across samples

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

#generate intermutation plot

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data
muts = fread("2019-11-27_READ_ONLY_ALL_MERGED_MUTS.txt")
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")
genes_class = read.xlsx("Figure2_data_ABC_GCB_gene_Dave.xlsx")
colnames(genes_class)[1] = "hg19.ensemblToGeneName.value"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#for each sample extract mutations associated with ABC or GCB subtype genes

samples = unique(muts$id)

get_class = function(sample){
  samp_dat = as.data.table(filter(muts, id == sample))
  samp_dat = merge(samp_dat, genes_class, by = "hg19.ensemblToGeneName.value")
  return(samp_dat)
}

all_class_genes = as.data.table(ldply(llply(samples, get_class)))
all_class_genes = all_class_genes[order(-Subtype.association)]
all_class_genes$hg19.ensemblToGeneName.value = factor(all_class_genes$hg19.ensemblToGeneName.value, levels=unique(all_class_genes$hg19.ensemblToGeneName.value))

#make heatmap
#x-axis = sample
#y-axis = gene
#fill = VAF? copy number? "gt_AF"  "Nmaj"
#colour lettering based on ABC GCB

p = ggplot(all_class_genes, aes(id, hg19.ensemblToGeneName.value)) +
  geom_tile(aes(fill = Subtype.association), colour = "grey50") +
  rotate_x_text(90) + xlab("Sample") + ylab("Gene") #+
  #scale_fill_gradient2(low = "black", mid = "grey", midpoint = 0.5,
  #                     high = "red", na.value="transparent")
p
ggsave(paste(date, "summary_COO_genes_across_samples.pdf"))
