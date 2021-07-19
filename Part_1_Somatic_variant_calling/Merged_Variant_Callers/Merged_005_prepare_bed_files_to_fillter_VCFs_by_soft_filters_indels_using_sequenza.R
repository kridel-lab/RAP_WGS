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
options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr", "openxlsx", "readxl", "GenomicRanges",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom")
lapply(packages, require, character.only = TRUE)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare input mutation files for pyclone, treeomics and phylowgs

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data
muts = readRDS(list.files(pattern="merged_mut_calls_indels.rds")[length(list.files(pattern="merged_mut_calls_indels.rds"))]) #get most recent mutation file

#blacklist mutations were already removed
#TLOD filter and minimum # of alternative alleles was also filtered
table(muts$Indiv)
unique(muts$Indiv)

#2. Sample summary
samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")
colnames(samps)[4] ="Indiv"
z = which(samps$Indiv %in% muts$Indiv)
samps = samps[z,]

muts = merge(muts, samps, by="Indiv", all=TRUE)
unique(muts$Indiv)

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
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")

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
table(muts_wCNAs$ntot)

#----------------------------------------------------------------------
#SEPERATE FILES INTO FINAL SETS BASED ON TOOL INPUT
#----------------------------------------------------------------------

#1. READ-ONLY FILE = FINAL MUTATIONS = KEEP THIS WAY UNLESS MAJOR CHANGE NEEDED
read_only = as.data.table(muts_wCNAs)
read_only$tot_counts = read_only$Ref_counts+read_only$alt_counts

#2. PYCLONE = REMOVE UNIQUE MUTATIONS
saveRDS(read_only, file=paste(date, "READ_ONLY_ALL_MERGED_MUTS_INDELS.rds", sep="_"))
patients = unique(read_only$STUDY_PATIENT_ID)
