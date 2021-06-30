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
hg19_genes = unique(hg19_genes[,c("chr", "start", "end", "symbol")]) #22810 genes
z = which(hg19_genes$chr %in% c(1:22, "X"))
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

muts_per_sample = as.data.table(table(cnas$Patient, cnas$Tissue_Site, cnas$CNA))
muts_per_sample = as.data.table(filter(muts_per_sample, N >0))
muts_per_sample = muts_per_sample[order(-N)]

colnames(muts_per_sample) = c("Patient", "Sample", "type_CNA", "num_of_muts")
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
muts_per_sample$Patient = factor(muts_per_sample$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

#get order of samples for barplot
t = as.data.table(table(cnas$Patient, cnas$Tissue_Site))
t = as.data.table(filter(t, N >0))
t = t[order(-N)]

muts_per_sample$Sample = factor(muts_per_sample$Sample, levels=unique(t$V2))
write.table(muts_per_sample, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure1_MAIN_CNAs_ALL.txt", quote=F, row.names=F, sep="\t")
