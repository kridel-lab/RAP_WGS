#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)
library("ggsci")
library(annotables)
hg19_genes = as.data.table(grch37)
#hg19_genes = filter(hg19_genes, biotype == "protein_coding")
hg19_genes = unique(hg19_genes[,c("chr", "start", "end", "symbol")])
z = which(hg19_genes$chr %in% c(1:22))
hg19_genes = hg19_genes[z,]#52955 genes

#----------------------------------------------------------------------
#overlap CNA calls with protein coding genes coordinates
#----------------------------------------------------------------------

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")
cnas$Patient = sapply(cnas$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
#cnas = filter(cnas, !(ntot==2))
cnas$CNA[cnas$ntot >2] = "Amplification"
cnas$CNA[cnas$ntot <2] = "Deletion"
cnas$cna_id = paste(cnas$CHROM, cnas$Start, cnas$End, sep="_")

cnas_gr = GRanges(
  seqnames = cnas$CHROM,
  ranges = IRanges(cnas$Start, end = cnas$End),
  strand = rep("*", length(cnas$Start)),
  depth_ratio = cnas$depth.ratio,
sample = cnas$Sample,
Nmin = cnas$Nmin,
Nmaj = cnas$Nmaj,
ntot = cnas$ntot,
cna_id = cnas$cna_id)

genes_test = unique(hg19_genes$symbol)

overlap_genes_cnas = function(gene){
  print(gene)
  gene_cords = as.data.table(filter(hg19_genes, symbol==gene))

  #make granges objects
  gene_cords$CHROM = paste("chr", gene_cords$chr, sep="")

  gene_cords_gr = GRanges(
    seqnames = gene_cords$CHROM,
    ranges = IRanges(gene_cords$start, end = gene_cords$end),
    strand = rep("*", length(gene_cords$start)),
    score = 1:length(gene_cords$start))

  #intersect them
  #Then subset the original objects with the negative indices of the overlaps:
  hits <- findOverlaps(gene_cords_gr, cnas_gr, ignore.strand=TRUE)
  hits_overlap = cbind(gene_cords[queryHits(hits),], cnas[subjectHits(hits),])
  #print(head(hits_overlap))
  return(hits_overlap)
}

#ONLY RUN ONCE

#all_genes_cnas_samples = as.data.table(ldply(llply(genes_test, overlap_genes_cnas, .progress="text")))

#keep only non-neutral CNAs
#all_genes_cnas_samples$Patient = sapply(all_genes_cnas_samples$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
#saveRDS(all_genes_cnas_samples, file="cnas_across_hg19_genomes_all_samples.rds")

all_genes_cnas_samples = readRDS("cnas_across_hg19_genomes_all_samples.rds")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#dir.create(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))
#setwd(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))

t=as.data.table(table(all_genes_cnas_samples$symbol, all_genes_cnas_samples$CNA, all_genes_cnas_samples$Sample, all_genes_cnas_samples$Patient))
t=filter(t, N >0)
tt=as.data.table(table(t$V1, t$V2, t$V4))
tt=filter(tt, N>0)
tt=tt[order(-N)]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

samples_per_mut = tt
colnames(samples_per_mut) = c("Gene", "CNA", "Patient", "num_of_samples_with_CNAinGene")

#data table for barplot
barplot = as.data.table(table(samples_per_mut$num_of_samples_with_CNAinGene, samples_per_mut$CNA, samples_per_mut$Patient))
barplot = as.data.table(filter(barplot, N >0))

colnames(barplot) = c("num_of_samples_with_mut", "CNA", "patient", "num_of_muts")
#barplot$num_of_samples_with_mut = factor(barplot$num_of_samples_with_mut, levels=unique(barplot$num_of_samples_with_mut))
barplot$patient[barplot$patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
barplot$patient[barplot$patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
barplot$patient[barplot$patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
barplot$patient = factor(barplot$patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

saveRDS(barplot, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure2_input_data_sample_dist_CNAs.rds")

#colnames(barplot)[2]="Patient"
#pdf("001_samples_per_mutation_lineplot_CNAs.pdf", width=5, height=5)
# Basic barplot
#p<-ggline(barplot, x="num_of_samples_with_mut", y="num_of_muts",
#palette = c("#00AFBB", "#E7B800", "#FC4E07"), color="Patient")+
#xlab("# of samples sharing mutation") + ylab("CNAs count")#+
#scale_y_continuous(breaks=seq(0,15000,1000))
#print(p)
#dev.off()
