#----------------------------------------------------------------------
#variants_003_read_in_VCFs_into_matrix.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
  "ggrepel", "stringr", "maftools", "VariantAnnotation")
lapply(packages, require, character.only = TRUE)
library(GenomicRanges)

date = Sys.Date()

print(date)
args = commandArgs(trailingOnly = TRUE) #patient ID 
index = args[1] 
print(index) 

setwd(paste("MANTA_RESULTS/", "MANTA_WORKDIR_", index, "_files/", "gatk/", index, ".sorted.dup.recal.cram.bam/results/variants", sep=""))

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#how did we get here?

#----------------------------------------------------------------------
#data 
#----------------------------------------------------------------------

vcf = "somaticSV.vcf.gz" #normalized annotated vcf files (n=131)

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(vcf){

  pat = index
  if(!(dim(fread(vcf))[1] == 3)){ #to prevent empty VCF files from throwing errors

  #read in VCF file 
  vcf_dat <- readVcf(vcf, "hg19")

  #extract genotype info into data table 
  vcf_dat_ranges = as.data.frame(rowRanges(vcf_dat))
  vcf_dat_ranges$id = rownames(vcf_dat_ranges)
  vcf_dat_info = info(vcf_dat)  
  vcf_dat_info$id = rownames(vcf_dat_info)
  vcf_dat = as.data.table(merge(vcf_dat_ranges, vcf_dat_info, by = "id"))

  #add patient ID 
  vcf_dat$patient = pat

  #filter to include SVs that passed 
  vcf_dat = as.data.table(filter(vcf_dat, FILTER == "PASS"))

  #turn SV coordinates into Granges object and intersect with genes 
  vcf_dat$MATEID = unlist(vcf_dat$MATEID)
  vcf_dat_coords = unique(vcf_dat[,c("seqnames", "start", "end", "MATEID", "REF", "ALT", "id", "SVTYPE", "IMPRECISE", "BND_DEPTH", "MATE_BND_DEPTH")])
  vcf_dat_coords$seqnames = paste("chr", vcf_dat_coords$seqnames, sep="")
  vcf_dat_coords = makeGRangesFromDataFrame(vcf_dat_coords)

  #overlap with gene IDs
  vcf_dat_coords = as.data.table(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), vcf_dat_coords))
  
  #get gene symbols
  gene_symb = as.data.table(org.Hs.egSYMBOL)

  vcf_dat_coords = merge(vcf_dat_coords, gene_symb, by = "gene_id")

  #retrun
  return(vcf_dat)
}
}

clean_up_001(vcf)
print("done")
print(index)
saveRDS(all_vcfs_text, file=paste("/cluster/projects/kridelgroup/FLOMICS/variant_analysis_folder/SNP_matrices_all_algortihms/", date, index, ".rds", sep="_"))
