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

setwd("/Users/kisaev/OneDrive - UHN/RAP_ANALYSIS")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

svs = fread(list.files(pattern="all_SVs_samples.txt")[length(list.files(pattern="all_SVs_samples.txt"))])
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")
colnames(svs)[18:20] = paste("GENE", colnames(svs)[18:20], sep="_")

#sample info
samp_info = fread("RAP_samples_information.txt")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#distribution of types of structural variation
pats_svs = as.data.table(table(svs$SVTYPE, svs$pat)) ; pats_svs = pats_svs[order(-N)]
pats_tot = as.data.table(table(svs$pat)) ; pats_tot = pats_tot[order(-N)]
pats_svs$V2 = factor(pats_svs$V2, levels = pats_tot$V1)

ggbarplot(pats_svs, x = "V2", y="N", fill="V1") +theme_bw() + 
  rotate_x_text(75) + ylab("Number of SVs") + xlab("Sample")
ggsave(paste(date, "Manta_SVs_summary_patients.pdf"))

#genes most impacted by SVs 
genes_svs = as.data.table(table(svs$gene, svs$pat, svs$SVTYPE)) ; genes_svs = genes_svs[order(-N)]
genes_svs = as.data.table(filter(genes_svs, N >0, !(V1 == "")))
theme_update(text = element_text(size=6))

p = ggplot(filter(genes_svs, V1 %in% morin$Gene), aes(V2, V1)) +
  geom_tile(aes(fill = V3), colour = "grey50") + 
  rotate_x_text(90) + xlab("Sample") + ylab("Gene")

#ggplotly(p)
print(p)
ggsave(paste(date, "Manta_SVs_summary_Morin_genes.pdf"))

write.csv(svs, paste(date,"RAP_WGS_SVs_all_with_annotations.csv", sep="_"), quote=F, row.names=F)

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
  other_dat = as.data.table(filter(svs, id == dat$MATEID))
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
new_just_bnds = as.data.table(matrix(nrow=1, ncol=33)) ; colnames(new_just_bnds)=c("sv1_CHR", "sv1_start", "sv1_end",
                                                                                 "sv1_REF", "sv1_ALT", "sv1_SVTYPE", "sv1_SVLEN", 
                                                                                 "sv1_BND_DEPTH", "sv1_MATE_BND_DEPTH", 
                                                                                "sv1_HOMLEN", "sv1_HOMSEQ", "sv1_gene_chr", 
                                                                                "sv1_gene_start", "sv1_gene_end", "sv1_gene_width",
                                                                                "sv1_gene_name", "patient",
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
  dat = dat[,c(1:3,5:6, 8:9, 11:14, 17:20, 22:26, 28:29, 31:32, 34:37, 40:43, 45)]
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

bnds_sum = unique(new_just_bnds[,c("patient", "sv1_SVTYPE", "gene_pairs", "sv_ID")])
nonbnds_sum = unique(all_nonbnds[,c("pat", "SVTYPE", "gene", "id")])
colnames(bnds_sum) = c("pat", "SVTYPE", "gene", "id")
all_svs_sum_plot = rbind(bnds_sum, nonbnds_sum)
z=which(str_detect(all_svs_sum_plot$gene, "MT"))
all_svs_sum_plot$gene[z] = "MT_genes"

library(RColorBrewer)
t=as.data.table(table(all_svs_sum_plot$gene))
t=t[order(-N)]
all_svs_sum_plot$gene = factor(all_svs_sum_plot$gene, levels = t$V1)

theme_update(text = element_text(size=6))
mypalette <- brewer.pal(9, "YlOrRd") # 6 for geom_tile, 3 for geom_scatterpie
p <- ggplot(aes(x = gene , y = pat), data = all_svs_sum_plot )+
  geom_tile(aes(fill=SVTYPE), colour = "grey50") + 
  rotate_x_text(90)
p
ggsave(paste(date, "Manta_SVs_summary_Morin_genes.pdf"))
ggplotly(p)

#save summary file used to generate plot
write.csv(all_svs_sum_plot, file=paste(date, "RAP_WGS_all_SVs_heatmap_plot.csv", sep="_"), quote=F, row.names=F)

#save main file for translocations 
write.csv(new_just_bnds, file=paste(date, "RAP_WGS_translocations_ONLY.csv", sep="_"), quote=F, row.names=F)

#save main file for all SVs other than translocations 
write.csv(all_nonbnds, file=paste(date, "RAP_WGS_nontranslocations_ONLY.csv", sep="_"), quote=F, row.names=F)

#save SV masterfile that has all details regarding all SVs 
write.csv(new_just_bnds, file=paste(date, "RAP_WGS_MANTA_masterlist.csv", sep="_"), quote=F, row.names=F)




