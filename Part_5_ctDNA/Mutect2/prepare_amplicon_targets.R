#----------------------------------------------------------------------
#prepare intervals for picard collect metrics tool
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/ctDNA/Capture panel/final_regions_052021")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr", "maftools")
library("readxl")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#Interval files provided/confirmed by Larry from OICR May 11 2021
#----------------------------------------------------------------------

#amplicons/probe or baits
amplicons = as.data.table(read_excel("19_09_03_TGL50_All_probes_Full_and_partial_genes_agena.xlsx"))
amplicons = amplicons[,c("Chromosome", "Start", "Stop")] #1,675 baits including 38 which correspond to the Agena sample identity probes (TargetNN)

#targets list for 46 genes
targets_pcg = fread("my-targets.bed")
targets_pcg = unique(targets_pcg[,c(1:3)]) #571 targets
targets_pcg$V1 = sapply(targets_pcg$V1, function(x){unlist(strsplit(x, "chr"))[2]})
colnames(targets_pcg) = c("Chr", "Start", "Stop")

#targets list for additional 152 regions
target_regs = as.data.table(read_excel("NGS-Targets.xlsx"))
target_regs = target_regs[,c("Chr", "Start", "Stop")] #1,675 baits including 38 which correspond to the Agena sample identity probes (TargetNN)

#combine targets
all_targets = rbind(targets_pcg, target_regs) #723 total target regions

colnames(amplicons) = c("chr", "start", "stop")
write.table(amplicons, file="picard_tools_amps_input.bed", quote=F,
row.names=F, col.names=F, sep="\t")

#collect all target intervals
colnames(all_targets) = c("chr", "start", "stop")
write.table(all_targets, file="picard_tools_targets_input.bed", quote=F,
row.names=F, col.names=F, sep="\t")
