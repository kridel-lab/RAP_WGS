#----------------------------------------------------------------------
#prepare intervals for picard collect metrics tool
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/ctDNA/Capture panel/46 genes_599 regions_129KB")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#Interval files provided by ?
#----------------------------------------------------------------------

#amplicons/probe or baits
amplicons = fread("designed-probe-coords.bed")
amplicons = unique(amplicons[,c(1:3)]) #1,337 baits
amplicons$V1 = sapply(amplicons$V1, function(x){unlist(strsplit(x, "chr"))[2]})

#targets list
targets = fread("my-targets.bed")
targets = unique(targets[,c(1:3)]) #571 targets
targets$V1 = sapply(targets$V1, function(x){unlist(strsplit(x, "chr"))[2]})

colnames(amplicons) = c("chr", "start", "stop")
write.table(amplicons, file="picard_tools_amps_input.bed", quote=F,
row.names=F, col.names=F, sep="\t")

#collect all target intervals
colnames(targets) = c("chr", "start", "stop")
write.table(targets, file="picard_tools_targets_input.bed", quote=F,
row.names=F, col.names=F, sep="\t")
