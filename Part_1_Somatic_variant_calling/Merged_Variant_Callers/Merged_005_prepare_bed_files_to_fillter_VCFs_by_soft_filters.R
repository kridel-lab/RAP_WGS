#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

#purpose: complete final soft filters 
#mainly intersect with copy number file 
#generate filters required for Treeomics, PhyloWGS and Pyclone 
#so that VCFs can be filtered accordingly 

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats", "ggthemes")
lapply(packages, require, character.only = TRUE)
library(RColorBrewer)
library(openxlsx)
library(plotly)

#display.brewer.all()
#display.brewer.pal(9, "Set1")

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
muts = readRDS(list.files(pattern="mut_calls.rds")[length(list.files(pattern="mut_calls.rds"))]) #get most recent mutation file

#blacklist mutations were already removed 
#TLOD filter and minimum # of alternative alleles was also filtered

table(muts$Indiv)

#2. Sample summary 
dna = fread("RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
biops = fread("RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
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
table(muts_wCNAs$ntot)

#----------------------------------------------------------------------
#SEPERATE FILES INTO FINAL SETS BASED ON TOOL INPUT 
#----------------------------------------------------------------------

#1. READ-ONLY FILE = FINAL MUTATIONS = KEEP THIS WAY UNLESS MAJOR CHANGE NEEDED 
#FINAL soft filter for this is based on CNAs - keep only variants with less than 5 total copy number 
#AND ALSO THOSE WITHOUT 0 COPIES BECAUSE THAT MEANS DE

#2. TREEOMICS INPUT VCF FILES = COPY NEUTRAL MUTATIONS ONLY 

#3. PHYLOWGS INPUT = REMOVE FOUNDER MUTATIONS AND UNIQUE MUTATIONS 

#4. PYCLONE = REMOVE FOUNDER MUTATIONS AND UNIQUE MUTATIONS
#ONLY COPY NEUTRAL MUTATIONS 
#OTHER FILTERS? DEPTH IS BIASED BECAUSE MIGHT NOT REACH THREHOLD IN DIAGNOSTIC BUT STILL BE 
#PRESENT AT LOWER DEPTH = MUTATIONA ACTUALLY STILL THERE 
#INCREASE THRESHOLD FOR AUTOPSY SAMPLES BUT KEEP THE SAME FOR DIAGNOSTIC? 











