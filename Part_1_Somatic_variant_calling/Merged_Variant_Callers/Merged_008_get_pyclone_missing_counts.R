#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

date = Sys.Date()
print(date)
options(scipen=999) #no scientific notation

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
source("/cluster/home/kisaev/scripts/bam_readcount_parseline.R")

#load libraries
packages <- c("dplyr", "ggplot2", "tidyr", "data.table", "plyr",
	"stringr")
lapply(packages, require, character.only = TRUE)
library("readxl")
library(GenomicRanges)
library(params)
library(readr)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#used bam readcount to extract counts from all samples
#overlapping mutations that were only found in some samples and not all
#Pyclone input requires these values

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#Our mutations
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone")
#sample info
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#bam readcount values for mutations not found in all patients
bamreadcount=list.files(pattern="missing")

#file with missing variants
miss_vars = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_bam_readcount_input.bed")
colnames(miss_vars) = c("chr", "start", "end",
																		 "ref_allele", "alt_allele")
write.table(miss_vars, "pyclone_bam_readcount_input.bed", quote=F, row.names=F, sep="\t")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#first combine all results from bamreadcount into one matrix
get_bam_readcount = function(file_res){

	#samplename
	samplename = unlist(strsplit(file_res, "_missing_muts"))[1]

	#readcount output
	x = file_res

	#get output
	output_t1reads_t2muts = as.data.table(bam_readcount.parse(x, samplename = samplename, "pyclone_bam_readcount_input.bed"))

	#save output
	return(output_t1reads_t2muts)
}

missing_mutations = as.data.table(ldply(llply(bamreadcount, get_bam_readcount)))

#unique mutations
unique = filter(as.data.table(table(read_only$mut_id)), N == 1)

pyclone_input = as.data.table(filter(read_only, !(mut_id %in% unique$V1),
!(Func.ensGene %in% c("ncRNA_intronic", "intronic", "intergenic"))))

diagnostic = as.data.table(filter(pyclone_input, Specimen_Type == "FFPE", alt_counts >=20)) #median alt count
autopsy = as.data.table(filter(pyclone_input, Specimen_Type == "FT", alt_counts >=32)) #median alt count

pyclone_input = rbind(diagnostic, autopsy) #1364 unique mutations...
t=filter(as.data.table(table(pyclone_input$mut_id)), !(N==1)) #remove unique mutations
z = which(pyclone_input$mut_id %in% t$V1)
pyclone_input = pyclone_input[z,] #1313 unique mutations

#for mutations that are not present in all samples need to generate an entry for them
#ideally need to get count of reads mapping there but for now just gonna put in zeros
t = filter(as.data.table(table(pyclone_input$mut_id)), (N==20))
muts_all = as.data.table(filter(pyclone_input, mut_id %in% t$V1))
muts_some = as.data.table(filter(pyclone_input, !(mut_id %in% t$V1)))

#save muts_some so that can run bam readcount and extract counts in those mutations
#across all samples
#need chr, start, end, ref and alt save as bed file, no colnames, tab sep
muts_some_bam_readcount = unique(muts_some[,c("CHROM", "POS", "mut_id", "REF", "ALT")])
muts_some_bam_readcount$mut_id = muts_some_bam_readcount$POS
#remove chr from CHROM
muts_some_bam_readcount$CHROM = sapply(muts_some_bam_readcount$CHR, function(x){
unlist(strsplit(x, "chr"))[2]
})

#for muts in muts_some ... need to generate a record for samples that dont have mutation
get_record = function(mutation){
  print(mutation)
  mut_dat = as.data.table(filter(muts_some, mut_id == mutation))

  #columns that we need to edit
  #c("mut_id", "Ref_counts", "alt_counts", "normal_cn", "Nmin", "Nmaj", "hg19.ensemblToGeneName.value", "Func.ensGene", "id")
  #which patients don't have
  pats = unique(read_only$id)[which(!(unique(read_only$id) %in% mut_dat$id))]

    make_pat = function(pat){
      pat_dat = mut_dat[1,]
      pat_dat$id = pat
      pat_dat$MajorCN=2
      pat_dat$MinorCN=0
      pat_dat$Ref_counts=0
      pat_dat$alt_counts=0
      return(pat_dat)
    }

  missing_pats = as.data.table(ldply(llply(pats, make_pat, .progress="text")))
  mut_dat = rbind(mut_dat, missing_pats)
  return(mut_dat)

}

all_muts = unique(muts_some$mut_id)
missing_records = as.data.table(ldply(llply(all_muts, get_record, .progress="text")))
all_records = rbind(muts_all, missing_records)

#check that now all mutations appear in all 20 samples
t = as.data.table(table(all_records$mut_id))
t=t[order(-N)]

write.table(all_records, file=paste(date, "PYCLONE_INPUT_MUTS.txt", sep="_"), quote=F, row.names=F, sep="\t")
