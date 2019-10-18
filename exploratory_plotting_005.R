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
muts = readRDS("2019-10-17_SSMs_went_into_phyloWGS.rds")

#2. Sample summary 
dna = fread("~/Documents/RAP_analysis/RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("~/Documents/RAP_analysis/RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
dna = merge(dna, biops, by="barcode")
colnames(dna)[2] = "Indiv"
colnames(dna)[7] = "Tissue_Site"
colnames(dna)[8] = "Specimen_Type" 
dna$Specimen_Type = "FT"

dna = as.data.table(filter(dna, Indiv %in% muts$Indiv))
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
dna$id = paste(dna$Specimen_Type, dna$Tissue_Site, dna$barcode, sep="_")

muts = merge(muts, dna, by="Indiv", all=TRUE)
#increase minimum depth to 60
muts = as.data.table(filter(muts, DP>=60, DP < 500))

#Mutation data from WGS Morin et al 2013
tables2 = fread("TableS2.csv")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

muts$Tissue_Site[is.na(muts$Tissue_Site)] = "FFPE"

#1. -------------------------------------------------------------------

#total mutations 
table(muts$Indiv)

heatmap = as.data.table(filter(muts, Gene.ensGene %in% tables2$`Ensembl ID`))

#make heatmap so summarize VAFs of mutations across the samples
heatmap = melt(heatmap, id.vars = c("Indiv", "Specimen_Type", "Tissue_Site", "mut_id", "hg19.ensemblToGeneName.value"),
     measure.vars = c("gt_AF"))

# Create the heatmap
p = ggplot(heatmap, aes(Tissue_Site, hg19.ensemblToGeneName.value, fill = value))+
  geom_tile(color="black") +
  scale_fill_gradient(low = "white", high = "steelblue", na.value="transparent") +
  theme_classic()
p + rotate_x_text(90) + rremove("y.ticks")#+
  #rremove("y.text")

