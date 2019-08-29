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
library(cowplot)

library(RColorBrewer)
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
#muts = readRDS("final_somatic_mutations_RAP_WGS_20_samples.rds")
muts = readRDS("2019-08-28_all_soft_filtered_SNVs_overlapping_titan_cna_calls.rds")

#2. Sample summary 
dna = fread("RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
dna = merge(dna, biops, by="barcode")
colnames(dna)[2] = "Indiv"
colnames(dna)[7] = "Tissue_Site"
colnames(dna)[8] = "Specimen_Type" 
dna$Specimen_Type = "FT"

dna = as.data.table(filter(dna, Indiv %in% muts$vcf_sample))
ffpe = as.data.table(matrix(ncol=ncol(dna), nrow=3))
colnames(ffpe) = colnames(dna)
ffpe = as.data.frame(ffpe)

#ffpe$Indiv = unique(muts$vcf_sample[which(!(muts$vcf_sample %in% dna$Indiv))])
ffpe$Indiv = c("LY_RAP_0003_Dia_FoT_05", "LY_RAP_0003_Dia_FoT_01" ,"LY_RAP_0003_Dia_FoT_03")
ffpe$barcode =c("15:S12966E", "15:S12966A", "15:S12966C")
ffpe$Tissue_Site = c("left_breast", "right_neck_LN", "left_axilla_LN")
ffpe$Specimen_Type = "FFPE"
ffpe$DNA = "DNA"
ffpe$STUDY_PATIENT_ID = "LY_RAP_0003"
dna = rbind(dna, ffpe)

colnames(dna)[2] = "vcf_sample"
muts = merge(muts, dna, by="vcf_sample", all=TRUE)

#Mutation data from WGS Morin et al 2013
tables2 = fread("TableS2.csv")
muts$id = paste(muts$Specimen_Type, muts$Tissue_Site, muts$barcode, sep="_")
muts$Cellular_Prevalence = as.numeric(muts$Cellular_Prevalence)

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. -------------------------------------------------------------------

#total mutations 
t = as.data.table(table(muts$id)) ; t = t[order(-N)]
#barplot
ggbarplot(t, x="V1", y="N") + coord_flip() + ylab("Somatic Mutations (N)") + xlab("Site")
ggsave(paste(date, "001_summary_mutations.pdf", sep = "_"))

#total CNAs
t = as.data.table(table(muts$Corrected_Call, muts$id)) ; t = t[order(-N)] ; t = as.data.table(filter(t, N >0))
colnames(t) = c("CNA_status", "Site", "N")
#barplot
ggbarplot(t, x="Site", y="N", fill="CNA_status") +
 coord_flip() + ylab("CNA status of SNVs (N)")
ggsave(paste(date, "002_summary_CNAs_overlapping_mutations.pdf", sep="_"))

#CNAs across chrosmomes 
t = as.data.table(table(muts$chr_snv, muts$Corrected_Call)) 
colnames(t) = c("Chromosome", "CNA_status", "N")
#barplot
ggbarplot(t, x="Chromosome", y="N", fill="CNA_status") +
  coord_flip() + ylab("CNA status (N)")
ggsave(paste(date, "003_summary_CNAs_overlapping_mutations_across_chromosomes.pdf", sep="_"))

#what kind of genes are affected 
t = as.data.table(table(muts$region, muts$id)) ;  t = t[order(-N)] ; t = as.data.table(filter(t, N >0))
colnames(t) = c("Region", "Site", "N")
totals = as.data.table(table(muts$region)) ; colnames(totals) = c("Region", "total")
t = merge(t, totals, by="Region")
t$fraction = t$N/t$total

#barplot
ggbarplot(t, x="Site", y="fraction", fill="Region") +
  coord_flip() + ylab("Regions affected by SNVs (%)")
ggsave(paste(date, "004_summary_SNVS_types_of_regions.pdf", sep="_"))

#which genes are mutated across all samples versus only some?
t = as.data.table(table(muts$Symbol, muts$id)) ; t = t[order(-N)] ; t = as.data.table(filter(t, N >0))
#remove Y=RNA

#2. -------------------------------------------------------------------
ezh2 = as.data.table(filter(muts, Symbol=="EZH2")) ; ezh2 = ezh2[order(Specimen_Type)]

ggplot(ezh2, aes(x=id, y=mutation_id)) + geom_tile(aes(fill=Corrected_Call)) +
coord_flip() +xlab("Region")
ggsave(paste(date, "005_summary_SNVS_CNAs_in_EZH2.pdf", sep="_"))

#3. ----------
#explore these "likely driver gene mutations" 
drivers = c("GNAS", "HIST1H1E", "KEL", "MSH3", "MUC6", "NUP133", "PLCB4", "TCEB1", "TLR4", "ZFP36L1", "BCL2", "BRAF", "CNBD1", "ESR1", "EZH2", "FAT1", "HLA-A", "IL7R", "JAK2", "KMT2C", "KRAS", "MTOR", "MYC", "MYD88", "NOTCH2", "PIM1", "PTEN", "PTPRC", "PTPRD", "SF3B1", "SPOP", "TBX3", "ZFHX3")
drivers = as.data.table(filter(muts, Symbol %in% drivers))
drivers$gene_mut = paste(drivers$Symbol, drivers$start_snv)
orderlist = as.data.table(table(drivers$Symbol)) ; orderlist=orderlist[order(-N)]
barplot = ggbarplot(orderlist, x="V1", y="N") + theme_light()+ rremove("x.text") + rremove("x.ticks") + rremove("xlab") + ylab("total muts(N)") 

#histo 
drivers$Symbol = factor(drivers$Symbol, levels = orderlist$V1)
drivers$CNA_Cellular_Prevalence = drivers$Cellular_Prevalence

g = ggplot(drivers, aes(x=id, y=Symbol)) + geom_tile(aes(fill=Corrected_Call,width=0.75, height=0.75),size=0.55, colour = "black") + theme_light() +
   xlab("Region")
g= ggpar(g, x.text.angle = 90, legend ="bottom") + coord_flip() + scale_fill_brewer(palette="Set3")

plot_grid(plotlist=list(barplot, g), ncol=1, 
            rel_heights=c(1,4), align="v")

ggsave(paste(date, "006_summary_SNVS_CNAs_in_potential_drivers.pdf", sep="_"))

g = ggplot(drivers, aes(x=id, y=Symbol)) + geom_tile(aes(fill=CNA_Cellular_Prevalence,width=0.75, height=0.75),size=0.55, colour = "black") + theme_light() +
  xlab("Region")
g= ggpar(g, x.text.angle = 90, legend ="bottom") + coord_flip() + scale_fill_gradient(low = "grey", high = "purple", na.value = 'white')

plot_grid(plotlist=list(barplot, g), ncol=1, 
          rel_heights=c(1,4), align="v")

ggsave(paste(date, "007_summary_SNVS_CNAs_in_potential_drivers_CP_from_TITAN.pdf", sep="_"))

#make gene specific geom_tile plot
function(gene){
  gene_dat = as.data.table(filter(muts, Symbol == gene)) 
}
