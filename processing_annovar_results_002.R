#----------------------------------------------------------------------
#processing_annovar_results.R
#karin isaev
#July 11th, 2019
#----------------------------------------------------------------------

date = Sys.Date()

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/vcfs_annovar_annotated/vcf_summary_text")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 was run on paired mode compaing cns to diagnostic tumour
#now it's time to:
#summarize cns specific mutations
#but first should still filter out false positives (note, these are unfilitered variants)
#see how many appear in multiple comparisons (n=5 total)

#note these vcf files have been normalized and fed through annovar 
#for annotations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

paired = list.files(pattern="vcf_file_filtered.bed")

#gene annotations
genes = unique(fread("/cluster/projects/kridelgroup/paired_cns/ucsc_table_browser_gene_IDs.txt"))
colnames(genes)[2] = "Gene.ensGene"

#sample data 
samp_dat = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#SSM input data for PHYLOWGS 
ssm_dat = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/ssm_data.txt")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#combine files into one matrix to make input file for PhyloWGS 
get_dat = function(file){
  f = fread(file)
  f = f[,c("Indiv", "mut_id")]
  return(f)
}

t = as.data.table(ldply(llply(paired, get_dat)))

#subset SSM input file for phylowgs to only include these mutatoins 

z = which(ssm_dat$gene %in% t$mut_id)
ssm_dat = ssm_dat[z,]
ssm_dat$id = paste("s", 0:(nrow(ssm_dat)-1), sep="")

write.table(ssm_dat, file="ssm_dat_35512_variants_filtered.txt", sep="\t", row.names=F, quote=F)

#-----------------------------------------------------------------------
#all mutation data summary for downstream analysis 
#-----------------------------------------------------------------------

get_dat = function(file){
  f = fread(file)
  return(f)
}

all_muts = as.data.table(ldply(llply(paired, get_dat)))
saveRDS(all_muts, file = paste(date, "SSMs_went_into_phyloWGS.rds", sep="_")) 
















