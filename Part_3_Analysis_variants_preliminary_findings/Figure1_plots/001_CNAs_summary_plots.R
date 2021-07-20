#-------------------------------------------------------------------------------
#001_CNAs_summary_plots.R
#This script preps the CNA data for Figure 1 plots
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)
library("ggsci")
library(annotables)
hg19_genes = as.data.table(grch37)
#hg19_genes = filter(hg19_genes, biotype == "protein_coding")
hg19_genes = unique(hg19_genes[,c("chr", "start", "end", "symbol")]) #22810 genes
z = which(hg19_genes$chr %in% c(1:22))
hg19_genes = hg19_genes[z,]

#----------------------------------------------------------------------
#overlap CNA calls with protein coding genes coordinates
#----------------------------------------------------------------------

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")
cnas$Patient = sapply(cnas$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
cnas$cna_id = paste(cnas$CHROM, cnas$Start, cnas$End, sep="_")

#summarize number of mutations per sample
z = which((cnas$Tissue_Site == "Adrenal gland, NOS") & (cnas$Patient == "LY_RAP_0003"))
cnas$Tissue_Site[z] = "Adrenal gland"

z = which((cnas$Tissue_Site == "Aorta, ascending, not specified \n\n") & (cnas$Patient == "LY_RAP_0001"))
cnas$Tissue_Site[z] = "Aorta, ascending"

z = which((cnas$Tissue_Site == "Kidney, NOS") & (cnas$Sample == "LY_RAP_0003_Aut_FzT_12"))
cnas$Tissue_Site[z] = "Kidney, NOS 2"

tum_info_save=unique(cnas[,c("Sample", "Patient", "Tissue_Site", "Ploidy", "Purity")])
tum_info_save$Patient[tum_info_save$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
tum_info_save$Patient[tum_info_save$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
tum_info_save$Patient[tum_info_save$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
write.table(tum_info_save, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure1_MAIN_purity_ploidy.txt", quote=F, row.names=F, sep="\t")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

dir.create(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))
setwd(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))

#add categories for CNAs
cnas$CNA = ""

#LOH
z = which((cnas$ntot ==  cnas$Nmaj) & (cnas$ntot >= 2))
cnas$CNA[z] = "Somatic LOH"

#Neutral
z = which((cnas$ntot !=  cnas$Nmaj) & (cnas$ntot == 2))
cnas$CNA[z] = "Neutral"

#Deletions
z = which((cnas$ntot ==  cnas$Nmaj) & (cnas$ntot == 0))
cnas$CNA[z] = "Homozygous Del"

z = which((cnas$ntot ==  cnas$Nmaj) & (cnas$ntot == 1))
cnas$CNA[z] = "Hemizygous Del"

#Amplications
z = which((cnas$CNA == "") & (cnas$ntot > 5))
cnas$CNA[z] = ">5_N_Gain"

z = which(cnas$CNA == "")
cnas$CNA[z] = paste(cnas$ntot[z], "_N_", "Gain", sep="")

z = which(is.na(cnas$ntot))
if(!(length(z) == 0)){
  cnas=cnas[-z,]
}

#prep data for plot
cnas$width = cnas$End - cnas$Start

#get total number of basepairs explored in each sample
tot_widths = as.data.table(cnas %>% dplyr::group_by(Sample) %>% dplyr::summarize(sum_widths = sum(width)))

#get total number Base Pairs affected by each category of CNAs
cnas_by_widths = as.data.table(cnas %>% dplyr::group_by(Sample, CNA) %>% dplyr::summarize(tot_bps_cna = sum(width)))
cnas_by_widths = merge(cnas_by_widths, tot_widths, by="Sample")
cnas_by_widths$Patient = sapply(cnas_by_widths$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
cnas_by_widths$Patient[cnas_by_widths$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
cnas_by_widths$Patient[cnas_by_widths$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
cnas_by_widths$Patient[cnas_by_widths$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
cnas_by_widths$Patient = factor(cnas_by_widths$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

write.table(cnas_by_widths, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure1_MAIN_CNAs_ALL.txt", quote=F, row.names=F, sep="\t")
