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
library(openxlsx)

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
muts = readRDS("2019-10-31_all_muts_merged_mut_calls.rds")
#blacklist mutations were already removed 

#additional filters 
#TLOD - min 10 - confidence score for variant being really somatic 
gghistogram(muts, "TLOD", bins=200)
muts = as.data.table(filter(muts, TLOD > 10)) 
print(length(unique(muts$mut_id)))
gghistogram(muts, "TLOD", bins=200)

#min counts for alternative reads 
gghistogram(muts, "alt_counts", bins=200)
muts = as.data.table(filter(muts, alt_counts > 10))
length(unique(muts$mut_id))
print(length(unique(muts$mut_id)))
gghistogram(muts, "alt_counts", bins=200)

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

#Mutation data from WGS Morin et al 2013
#Table S3. All somatic SNVs identified from 40 genome pairs and 13 cell lines (XLSX, 918 KB)
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")

#Copy number data 
cnas = readRDS()

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

muts$Tissue_Site[is.na(muts$Tissue_Site)] = "FFPE"

#1. -------------------------------------------------------------------

#total mutations 
sum_muts = as.data.table(table(muts$id)) ; sum_muts = sum_muts[order(-N)]
print(sum_muts)
ggbarplot(sum_muts, x = "V1", y="N") + rotate_x_text(90)

#how many founder mutatiotns 
founds = filter(as.data.table(table(muts$mut_id, muts$Indiv)), N >0)
founds = as.data.table(filter(as.data.table(table(founds$V1)), N ==20)) #8268 unique variants present in everyone, 18902/52302 = 36%

#where are they?
founds = (as.data.table(filter(muts, mut_id %in% founds$V1)))
founds = unique(founds[,c("mut_id", "Func.ensGene", "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene", "cosmic68", "hg19.ensemblToGeneName.value")])

#any of the genes from Morin paper in founders?
morin$chr = sapply(morin$Chromosome, function(x){unlist(strsplit(x, "chr"))[2]})
morin$mut_id = paste(morin$chr, morin$Position, sep="_")

#how many unique mutations per patient? 
unique = filter(as.data.table(table(muts$mut_id, muts$id)), N >0)
unique_muts = as.data.table(filter(as.data.table(table(unique$V1)), N ==1)) #9133 unique variants present, 14800/52302 = 28%
unique_muts = merge(unique_muts, unique, by = c("V1", "N"))
unique_muts_sum = as.data.table(table(unique_muts$V2)) ;  unique_muts_sum = unique_muts_sum[order(-N)]
print(unique_muts_sum)
ggbarplot(unique_muts_sum, x = "V1", y="N") + rotate_x_text(90)

#where are they?
unique_muts = (as.data.table(filter(muts, mut_id %in% unique_muts$V1)))
unique_muts = unique(unique_muts[,c("mut_id", "Func.ensGene", "Gene.ensGene", "GeneDetail.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene", "cosmic68", "hg19.ensemblToGeneName.value")])
genes = as.data.table(table(unique_muts$hg19.ensemblToGeneName.value))

#save final list of mutations to filter out VCF files 
muts_keep = unique(muts[,c("CHROM", "POS")])
muts_keep = muts_keep[order(CHROM, POS)]
muts_keep$pos_end = muts_keep$POS
write.table(muts_keep, paste(date, "final_SNVs_include_in_VCFs.bed", sep="_"), quote=F, row.names=F, sep="\t", col.names = F)
saveRDS(muts, file=paste(date, file="final_list_of_mutations_input_palimpsest.rds", sep="_"))

#2. -------------------------------------------------------------------

#cluster samples based on mutation profiles, see which samples are more related to each other than others 
#generate binary matrix with 0s and 1s or VAFs? 

#remove founder and unique variants since they are not informative for this task 
muts_matrix = as.data.frame(dcast(muts, mut_id ~ id, value.var = "gt_AF"))
muts_matrix = subset(muts_matrix, !(mut_id %in% founds$mut_id))
muts_matrix = subset(muts_matrix, !(mut_id %in% unique_muts$mut_id))

rownames(muts_matrix) = muts_matrix$mut_id
muts_matrix$mut_id = NULL

#replace all NAs with zeros for now
muts_matrix[is.na(muts_matrix)] = 0

require(vegan)
muts_matrix = t(muts_matrix)
dist.mat<-vegdist(muts_matrix,method="jaccard", na.rm = TRUE)
clust.res<-hclust(dist.mat)
plot(clust.res)







