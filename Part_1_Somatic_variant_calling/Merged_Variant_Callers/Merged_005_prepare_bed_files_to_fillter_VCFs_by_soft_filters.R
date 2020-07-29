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
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom")
lapply(packages, require, character.only = TRUE)
library(RColorBrewer)
library(openxlsx)
library(plotly)
library(readxl)
library(GenomicRanges)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare input mutation files for pyclone, treeomics and phylowgs

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
#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))
treeomics_drivers = reddy
colnames(treeomics_drivers)[1] = "Gene_Symbol"
write.csv(treeomics_drivers, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/reddy_drivers.csv",
quote=F, row.names=F)

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_TITAN.rds")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#overlap SNVs with CNAs

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
    ranges = IRanges(cnas_pat$Start, end = cnas_pat$End),
    strand = rep("*", length(cnas_pat$Start)),
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
table(muts_wCNAs$Copy_Number)

#----------------------------------------------------------------------
#SEPERATE FILES INTO FINAL SETS BASED ON TOOL INPUT
#----------------------------------------------------------------------

#1. READ-ONLY FILE = FINAL MUTATIONS = KEEP THIS WAY UNLESS MAJOR CHANGE NEEDED
read_only = as.data.table(muts_wCNAs)
write.table(read_only, file=paste(date, "READ_ONLY_ALL_MERGED_MUTS.txt", sep="_"), quote=F, row.names=F, sep="\t")

#2. PHYLOWGS INPUT = REMOVE UNIQUE MUTATIONS
founds = filter(as.data.table(table(read_only$mut_id)), N == 20)
unique = filter(as.data.table(table(read_only$mut_id)), N == 1)
phylowgs_input = as.data.table(filter(read_only, !(mut_id %in% unique$V1),
!(Func.ensGene %in% c("ncRNA_intronic", "intronic", "intergenic"))))
phylowgs_input = unique(phylowgs_input[,c("CHROM", "POS")])
#remove CHR
phylowgs_input$CHROM = sapply(phylowgs_input$CHROM, function(x){unlist(strsplit(x, "chr"))[2]})

#set.seed(100)
#z = sample(dim(phylowgs_input)[1], 10000)
#phylowgs_input = phylowgs_input[z,]
write.table(phylowgs_input, file=paste(date, "PHYLOWGS_INPUT_MUTS.bed", sep="_"), quote=F, row.names=F, sep="\t")

#3. PYCLONE = REMOVE UNIQUE MUTATIONS
#REMOVE MUTATION WITH NMAJ OF 0
#KEEP WES GENE MUTATIONS TO SIMPLIFY

#run one version of pyclone with all mutations except for unique ones
#and those with major copy number greater than 0
pyclone_full = as.data.table(filter(read_only, !(mut_id %in% unique$V1),
MajorCN > 0))

#also remove noncoding mutations here, will make analysis easier
pyclone_input = as.data.table(filter(read_only, !(mut_id %in% unique$V1),
!(Func.ensGene %in% c("ncRNA_intronic", "intronic", "intergenic")), MajorCN > 0))

#diagnostic = as.data.table(filter(pyclone_input, Specimen_Type == "FFPE", alt_counts >=20)) #median alt count
#autopsy = as.data.table(filter(pyclone_input, Specimen_Type == "FT", alt_counts >=32)) #median alt count
#pyclone_input = rbind(diagnostic, autopsy) #1364 unique mutations...

length(unique(pyclone_full$mut_id)) #38605 mutations
length(unique(pyclone_input$mut_id)) #1466 mutations

#for mutations that are not present in all samples need to generate an entry for them
#ideally need to get count of reads mapping there but for now just gonna put in zeros

#get input for bamreadcount
get_bam_read = function(dat, name_analysis){

t = filter(as.data.table(table(dat$mut_id)), (N==20))
muts_all = as.data.table(filter(dat, mut_id %in% t$V1))
muts_some = as.data.table(filter(dat, !(mut_id %in% t$V1)))

#save muts_some so that can run bam readcount and extract counts in those mutations
#across all samples
#need chr, start, end, ref and alt save as bed file, no colnames, tab sep
muts_some_bam_readcount = unique(muts_some[,c("CHROM", "POS", "mut_id", "REF", "ALT")])
muts_some_bam_readcount$mut_id = muts_some_bam_readcount$POS
#remove chr from CHROM
muts_some_bam_readcount$CHROM = sapply(muts_some_bam_readcount$CHR, function(x){
unlist(strsplit(x, "chr"))[2]
})
print(dim(muts_some_bam_readcount))
file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/", name_analysis, "_pyclone_bam_readcount_input.bed", sep="")
write.table(muts_some_bam_readcount,
  file,
  col.names=F, quote=F, row.names=F, sep="\t")
}

get_bam_read(pyclone_full, "all_muts") #14317 mutations
get_bam_read(pyclone_input, "small_subset") #555 mutations 
