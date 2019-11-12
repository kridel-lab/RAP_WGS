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
library(plotly)


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
table(muts$Indiv)

#min counts for alternative reads 
gghistogram(muts, "alt_counts", bins=200)
muts = as.data.table(filter(muts, alt_counts > 10))
length(unique(muts$mut_id))
print(length(unique(muts$mut_id)))
gghistogram(muts, "alt_counts", bins=200)
table(muts$Indiv)

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
cnas = fread("copy_number_alteration_data_palimpsest_input.txt")

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
#write.table(muts_keep, paste(date, "final_SNVs_include_in_VCFs.bed", sep="_"), quote=F, row.names=F, sep="\t", col.names = F)
#saveRDS(muts, file=paste(date, file="final_list_of_mutations_input_palimpsest.rds", sep="_"))

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

#3. -------------------------------------------------------------------
#which genes most mutated ? 
#in founders? unique samples? across the board?

#1. in founds 
genes_founds = as.data.table(table(founds$hg19.ensemblToGeneName.value)) ; genes_founds = genes_founds[order(-N)]

#2. in unique 
genes_unique = as.data.table(table(unique_muts$hg19.ensemblToGeneName.value)) ; genes_unique = genes_unique[order(-N)]

#3. across the board 
genes_all = as.data.table(table(muts$hg19.ensemblToGeneName.value)) ; genes_all = genes_all[order(-N)]

#4. -------------------------------------------------------------------
#just Y-RNAs? what's going on with them?
all_Y_RNA = as.data.table(filter(muts, hg19.ensemblToGeneName.value == "Y_RNA"))
ensgs = as.data.table(table(all_Y_RNA$Gene.ensGene))
ensgs = ensgs[order(-N)]
y_rna_muts = as.data.table(table(all_Y_RNA$mut_id))
y_rna_muts = y_rna_muts[order(-N)]

#5. -------------------------------------------------------------------
#remove founder mutations for phyloWGS input
phylo_input_current = fread("ssm_data.txt")
phylo_input_current = as.data.table(filter(phylo_input_current, !(gene %in% founds$mut_id)))

#overlap SNVs with CNAs 
library(GenomicRanges)

overlap_snvs_cnas = function(sample){
  print(sample)
  snvs = as.data.table(filter(muts, Indiv==sample))
  cnas_pat = as.data.table(filter(cnas, Sample == sample))
  
  #make granges objects
  snvs$CHROM = paste("chr", snvs$CHROM, sep="")
  snvs_gr = GRanges(
    seqnames = snvs$CHROM,
    ranges = IRanges(snvs$POS, end = snvs$POS),
    strand = rep("*", length(snvs$POS)),
    score = 1:length(snvs$POS))
  
  cnas_gr = GRanges(
    seqnames = cnas_pat$CHROM,
    ranges = IRanges(cnas_pat$POS_START, end = cnas_pat$POS_END),
    strand = rep("*", length(cnas_pat$POS_START)),
    logR = cnas_pat$LogR)
    
  #intersect them
  #Then subset the original objects with the negative indices of the overlaps:
  
  hits <- findOverlaps(snvs_gr, cnas_gr, ignore.strand=TRUE)
  hits_overlap = cbind(snvs[queryHits(hits),], cnas_pat[subjectHits(hits),])
  print(head(hits_overlap))
  return(hits_overlap)
}

muts_wCNAs = as.data.table(ldply(llply(unique(muts$Indiv), overlap_snvs_cnas, .progress = "text")))
muts_wCNAs = muts_wCNAs[order(CHROM, POS)]

y_rnas_wcnas = as.data.table(filter(muts_wCNAs, mut_id %in% y_rna_muts$V1))

#plot just founders 
founds_wCNAs = as.data.table(filter(muts_wCNAs, mut_id %in% founds$mut_id))

pdf("Exploratory_Analysis/all_founder_mutations_wCNAs.pdf", width=25)
g = ggplot(founds_wCNAs, aes(x=mut_id, y=id)) + geom_tile(aes(fill=ntot)) +
  xlab("Region")
g= ggpar(g, x.text.angle = 90, legend ="bottom") + scale_fill_gradient2(low = "blue", mid = "white",
                                                                        high = "red", midpoint = 2) +
  rremove("ticks")+
  rremove("x.text")
print(g)
#ggplotly(g)
dev.off()

#plot just unique muts 
unique_wCNAs = as.data.table(filter(muts_wCNAs, mut_id %in% unique_muts$mut_id))

pdf("Exploratory_Analysis/all_unique_mutations_wCNAs.pdf", width=25)
g = ggplot(unique_wCNAs, aes(x=mut_id, y=id)) + geom_tile(aes(fill=ntot)) +
  xlab("Region")
g= ggpar(g, x.text.angle = 90, legend ="bottom") + scale_fill_gradient2(low = "blue", mid = "white",
                                                                        high = "red", midpoint = 2) +
  rremove("ticks")+
  rremove("x.text")
print(g)
#ggplotly(g)
dev.off()

#CNAs 
get_plot = function(sample){
  cnas_sam = as.data.table(filter(cnas, Sample == sample))
  cnas_sam$mut_id = paste(cnas_sam$CHROM, cnas_sam$POS_START, sep="_")
  
  g = ggplot(cnas_sam, aes(x=mut_id, y=ntot)) + geom_tile(aes(fill=LogR)) +
    xlab("Region")
  
  g= g + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 2) +
    rremove("ticks")+
    rremove("x.text") + ggtitle(sample)
  print(g)
  #ggplotly(g)
}

#llply(unique(muts$Indiv), get_plot)

