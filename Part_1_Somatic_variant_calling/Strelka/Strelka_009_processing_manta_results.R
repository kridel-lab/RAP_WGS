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
library(ggpubr)

date = Sys.Date()

print(date)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Manta")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

svs = fread(list.files(pattern="all_SVs_samples.txt")[length(list.files(pattern="all_SVs_samples.txt"))])
colnames(svs)[18:20] = paste("GENE", colnames(svs)[18:20], sep="_")

#sample info
samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")
colnames(samps)[4] ="Indiv"
z = which(samps$Indiv %in% svs$pat)
samps = samps[z,]
samps[18,2] = "Kidney, NOS 2"

colnames(svs)[which(colnames(svs)=="pat")] = "Indiv"
svs = merge(svs, samps, by = "Indiv")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#genes most impacted by SVs
genes_svs = as.data.table(table(svs$gene, svs$Indiv, svs$SVTYPE)) ; genes_svs = genes_svs[order(-N)]
genes_svs = as.data.table(filter(genes_svs, N >0, !(V1 == "")))

#for structural variants only (those that have "mates") - restructure file such that
#main id and mate appear in same row --> easy look up and visualization

all_bnds = unique(svs$id[svs$SVTYPE == "BND"])

svs$chromosome_name = NULL
svs$start_position = NULL
svs$end_position = NULL
svs$geneid = NULL

get_row = function(bnd){
  print(bnd)
  dat = as.data.table(filter(svs, id == bnd))
  other_dat = as.data.table(filter(svs, id %in% dat$MATEID))
  colnames(other_dat) = paste("sv2", colnames(other_dat), sep="_")
  new_bnd_entry = cbind(dat, other_dat)
  return(new_bnd_entry)
}

just_bnds = as.data.table(ldply(llply(all_bnds, get_row)))
just_bnds$gene_pairs = paste(just_bnds$gene, just_bnds$sv2_gene, sep="_")
just_bnds$chr_pairs = paste(just_bnds$SV_CHR, just_bnds$sv2_SV_CHR, sep="_")

#want to have just one gene per breakpoint coordinates
#and one entry per BND

all_bnds = unique(just_bnds$id)
new_just_bnds = as.data.table(matrix(nrow=1, ncol=35))
colnames(new_just_bnds)=c("sv1_CHR", "sv1_start", "sv1_end",
                          "sv1_REF", "sv1_ALT", "sv1_SVTYPE", "sv1_SVLEN",
                          "sv1_BND_DEPTH", "sv1_MATE_BND_DEPTH",
                          "sv1_HOMLEN", "sv1_HOMSEQ", "sv1_gene_chr",
                          "sv1_gene_start", "sv1_gene_end", "sv1_gene_width",
                          "sv1_gene_name", "sample", "tissue_site", "patient",
                          "sv2_CHR", "sv2_start", "sv2_end",
                          "sv2_REF", "sv2_ALT", "sv2_SVTYPE", "sv2_SVLEN",
                          "sv2_BND_DEPTH", "sv2_MATE_BND_DEPTH",
                          "sv2_HOMLEN", "sv2_HOMSEQ", "sv2_gene_chr",
                          "sv2_gene_start", "sv2_gene_end", "sv2_gene_width",
                          "sv2_gene_name")

for(i in 1:length(all_bnds)){
  bnd = all_bnds[i]
  print(bnd)
  dat = as.data.table(filter(just_bnds, id == bnd | MATEID == bnd))
  dat=dat[order(SV_CHR, SV_start, -sv2_gene)]
  dat = dat[1,]
  dat = dat[,c(2:4,6:7, 9:10, 12:15, 18:21, 23, 1, 25, 24, 28:30, 32:33, 35:36, 38:41, 44:47, 49)]
  colnames(dat) = colnames(new_just_bnds)
  new_just_bnds = rbind(new_just_bnds, dat)
}

new_just_bnds = new_just_bnds[-1,]

z = which(new_just_bnds$sv1_gene_name == "")
if(!(length(z)==0)){
new_just_bnds$sv1_gene_name[z] = NA}

z = which(new_just_bnds$sv2_gene_name == "")
if(!(length(z)==0)){
new_just_bnds$sv2_gene_name[z] = NA}

new_just_bnds$gene_pairs = paste(new_just_bnds$sv1_gene_name, new_just_bnds$sv2_gene_name, sep="_")
new_just_bnds$sv_ID = paste(new_just_bnds$sv1_CHR, new_just_bnds$sv1_start, new_just_bnds$sv1_end,
                            new_just_bnds$sv2_CHR, new_just_bnds$sv2_start, new_just_bnds$sv2_end, sep="_")
new_just_bnds = unique(new_just_bnds)
write.csv(new_just_bnds, file=paste(date, "RAP_WGS_BNDs_ONLY_with_annotations.csv", sep="_"), quote=F, row.names=F)

#summarize all other SVs that aren't BNDs

all_nonbnds = as.data.table(filter(svs, SVTYPE != "BND"))

#for each unique SV (not including BNDs) - if multiple genes mapped - merge into one entry and
#keep gene coordinates of just the first gene

unique_nonbnds = unique(all_nonbnds$id)
clean = function(sv){
  dat=as.data.table(filter(all_nonbnds, id == sv))
  genes=paste(dat$gene, collapse=",")
  dat$gene = genes
  dat=dat[order(GENE_start)]
  dat$GENE_end = dat$GENE_end[nrow(dat)]
  dat=dat[1,]
  return(dat)
}

all_nonbnds = as.data.table(ldply(llply(unique_nonbnds, clean)))

#make final heatmap summary of genes with mutations
#clean up BND only matrix and nonBND matrix
#include only patient, type of SV and genes affected

bnds_sum = unique(new_just_bnds[,c("sample", "tissue_site", "patient", "sv1_SVTYPE", "gene_pairs", "sv_ID")])
nonbnds_sum = unique(all_nonbnds[,c("Indiv", "Tissue_Site", "STUDY_PATIENT_ID", "SVTYPE", "gene", "id")])
colnames(bnds_sum) = c("Indiv", "Tissue_Site", "STUDY_PATIENT_ID", "SVTYPE", "gene", "id")

all_svs_sum_plot = rbind(bnds_sum, nonbnds_sum)
z=which(str_detect(all_svs_sum_plot$gene, "MT"))
all_svs_sum_plot$gene[z] = "MT_genes"

t=as.data.table(table(all_svs_sum_plot$gene))
t=t[order(-N)]
all_svs_sum_plot$gene = factor(all_svs_sum_plot$gene, levels = t$V1)

#save summary file used to generate plot
saveRDS(all_svs_sum_plot, file=paste(date, "RAP_WGS_all_SVs_heatmap_plot.rds", sep="_"))

#save main file for translocations
write.csv(new_just_bnds, file=paste(date, "RAP_WGS_translocations_ONLY.csv", sep="_"), quote=F, row.names=F)

#save main file for all SVs other than translocations
write.csv(all_nonbnds, file=paste(date, "RAP_WGS_nontranslocations_ONLY.csv", sep="_"), quote=F, row.names=F)

#save SV masterfile that has all details regarding all SVs
write.csv(new_just_bnds, file=paste(date, "RAP_WGS_MANTA_masterlist.csv", sep="_"), quote=F, row.names=F)
