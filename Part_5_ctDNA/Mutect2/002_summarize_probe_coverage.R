#----------------------------------------------------------------------
#Summarize probe coverage from picard
#----------------------------------------------------------------------

#Karin Isaev
#Started August 5th 2020
#tested on R version 3.5.0

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#this script takes output from picard collect targeted pcr metrics
#merges it with sample information
#make summary csv file with coverage per gene per patient
#simple barplot summaries of coverage per gene/probe/region
#calculate differences in coverage between coding and non-coding regions

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
"plyr",
"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)
library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)

#----------------------------------------------------------------------
#Load sample info
#----------------------------------------------------------------------

all_targets = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/my-targets.bed")
all_probes = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/designed-probe-coords.bed")

colnames(all_targets)[1:4]=c("chr", "target_start", "stop", "target_gene")
all_targets$chr = sapply(all_targets$chr, function(x){unlist(strsplit(x, "chr"))[2]})

#----------------------------------------------------------------------
#load in output files from picard tools
#----------------------------------------------------------------------

#all files with coverage per target
all_res = list.files(pattern="per_target_coverage.txt")

#read in all files and append sample names
read_file = function(file_name){

  f = paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls/", file_name, sep="")
  f = fread(f)
  sample = unlist(strsplit(file_name, ".per"))[1]
  f$Library = sample
  return(f)
}

all_res = as.data.table(ldply(llply(all_res, read_file)))
all_res$id = paste(all_res$chrom, all_res$end, sep="_")
all_targets$id = paste(all_targets$chr, all_targets$stop, sep="_")

dim(all_res)
dim(all_targets)

all_res = merge(all_res, all_targets, by = "id")
dim(all_res)

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#save all_res
write.csv(all_res, file=paste(date,
  "picard_tools_coverage_summary_targets_DNA_sequencing.csv", sep="_"), quote=F, row.names=F)

#calculate mean coverage between coding and non-coding
target_regions = as.data.table(all_res %>% group_by(target_gene) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))

#to do ====================================

all_res$gene = sapply(all_res$target_gene, function(x){unlist(strsplit(x, "\\("))[2]})
all_res$gene = sapply(all_res$gene, function(x){unlist(strsplit(x, ")"))[1]})

#calculate mean coverage by gene
gene_cov = as.data.table(all_res %>% group_by(gene, Library) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))
gene_cov = gene_cov[order(-mean_cov)]
gene_cov$gene = factor(gene_cov$gene, levels=unique(gene_cov$gene))
write.csv(gene_cov, file=paste(date, "gene_based_coverage_summary.csv", sep="_"),
  quote=F, row.names=F)

#mean target based coverage
target_cov = as.data.table(all_res %>% group_by(target_gene, Library) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))
target_cov = target_cov[order(-mean_cov)]
target_cov$target_gene = factor(target_cov$target_gene, levels=unique(target_cov$target_gene))
write.csv(target_cov, file=paste(date, "target_based_coverage_summary.csv", sep="_"),
  quote=F, row.names=F)

#calculate mean coverage by sample
sample_cov = as.data.table(all_res %>% group_by(Library) %>% dplyr::summarize(mean_cov = mean(mean_coverage)))
sample_cov = sample_cov[order(-mean_cov)]
write.csv(sample_cov, file=paste(date, "sample_based_coverage_summary.csv", sep="_"),
  quote=F, row.names=F)

#plot barplot summary per gene coverage
pdf("ctDNA_coverage_detailed_summaries.pdf",
width=20)

#mean gene based coverage
ggplot(gene_cov, aes(x=gene, y=mean_cov, fill=Library))+
geom_bar(stat="identity", color="black", position=position_dodge())+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#mean target based coverage
ggplot(target_cov, aes(x=target_gene, y=mean_cov, fill=Library))+
geom_bar(stat="identity")+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4))

#mean sample based coverage
ggplot(sample_cov, aes(x=Library, y=mean_cov, fill=Library))+
geom_bar(stat="identity", color="black")+
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
