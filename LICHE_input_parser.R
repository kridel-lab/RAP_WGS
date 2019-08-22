#----------------------------------------------------------------------
#phylowgs_input_parser.R
#karin isaev
#July 11th, 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr", 
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final")

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

#remove X and Y chromsomes and indels

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

paired = list.files(pattern="final_vcf_file_filtered.bed")

get_mut = function(file){
  f = fread(file)
  return(f)
}

all_muts = (llply(paired, get_mut))

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#filter to only keep SNVs (no indels) and DP > 60

#input for VCF data for phylowgs 
#ID (S0,S1,S2....), gene(chr_start) , a (#of reference alleles), d (#of total reads), mu_r (0.999), mu_v (0.499)

get_input = function(f){
  colnames(f) =  c("Gene.ensGene" ,"ChromKey" ,                   
 "POS"  ,                        "CHROM",                       
 "mut_id",                       "Indiv" ,                      
 "gt_AD"  ,                      "gt_AF"  ,                     
 "gt_DP"   ,                     "gt_F1R2" ,                    
"gt_F2R1"   ,                   "gt_GQ"    ,                   
 "gt_GT"     ,                   "gt_MBQ"    ,                  
 "gt_MFRL"    ,                  "gt_MMQ" ,                     
 "gt_MPOS"     ,                 "gt_PGT" ,                     
 "gt_PID"       ,                "gt_PL"  ,                     
 "gt_SA_MAP_AF"  ,               "gt_SA_POST_PROB",             
 "gt_GT_alleles"  ,              "ID"  ,                        
 "REF"             ,             "ALT"  ,                       
 "QUAL"             ,            "FILTER" ,                     
 "DP"                ,           "ECNT" ,                       
 "IN_PON"             ,          "NLOD" ,                       
 "N_ART_LOD"           ,         "POP_AF" ,                     
 "P_CONTAM"             ,        "P_GERMLINE",                  
 "RPA"                   ,       "RU" ,                         
 "STR"                    ,      "TLOD" ,                       
 "ANNOVAR_DATE"            ,     "Func.ensGene" ,               
 "GeneDetail.ensGene"       ,    "ExonicFunc.ensGene",          
 "AAChange.ensGene"          ,   "AF",                          
 "AF_popmax"                  ,  "AF_male",                    
 "AF_female"                   , "AF_raw",                      
 "AF_afr"                       ,"AF_sas",                     
 "AF_amr"                    ,   "AF_eas",                     
 "AF_nfe"                     ,  "AF_fin",                      
 "AF_asj"                     ,  "AF_oth",                     
 "non_topmed_AF_popmax"       ,  "non_neuro_AF_popmax" ,        
 "non_cancer_AF_popmax"       ,  "controls_AF_popmax",          
 "cosmic68"                   ,  "avsnp142",              
 "ALLELE_END"                 ,  "#hg19.ensGene.chrom" ,        
 "hg19.ensemblToGeneName.value")
  f = as.data.table(f)
  f= as.data.table(filter(f, DP >=100))
  #keep only SNVs 
  f = as.data.table(filter(f, ((ALT %in% c("A", "C", "G", "T" ) )& (REF %in% c("A", "C", "G", "T")))))
  #get reads 
  f = f %>% separate(gt_AD, c("ref", "alt"))
  final = f[,c("POS", "mut_id", "ref", "alt", "REF", "ALT", "CHROM", "POS")]
  colnames(final) = c("ID", "gene", "a", "d", "REF", "ALT", "CHROM", "POS")
  final$mu_r = 0.999
  final$mu_v = 0.499
  final$d = as.numeric(final$a) + as.numeric(final$d)
  final$patient = f$Indiv[1]
  final$ID = paste("S", 0:(nrow(final)-1), sep="")
  return(final)
} 

ssm_files = as.data.table(ldply(llply(all_muts, get_input)))

#remove mutations that appear "twice" in same patient
t = as.data.table(table(ssm_files$gene, ssm_files$patient))
rm = filter(t, N >=2)$V1
ssm_files = as.data.table(filter(ssm_files, !(gene %in% rm)))
dim(ssm_files)
length(unique(ssm_files$gene))

#actual SSM file generated by parser using VCFs - read in and keep only variants obtained here (more in-depth filtering)
#then write new SSM file and use for phylowgs run 
real_ssm = fread("ssm_data.txt")
z = which(real_ssm$gene %in% ssm_files$gene)
unique(real_ssm$gene[z])
real_ssm = real_ssm[z,]
#rename IDs
real_ssm$id = paste("s", 1:nrow(real_ssm), sep="")
write.table(real_ssm, file="ssm_data_additional_soft_filters.txt", quote=F, row.names=F, sep="\t") #unique 33120
ssm_files = as.data.table(filter(ssm_files, gene %in% real_ssm$gene))

#find mutations that are not present in all samples
#need to fill in the blank --> retrieve reads from BAM files
t = as.data.table(table(ssm_files$gene))
need_to_retrieve = as.data.table(filter(t, N <20))

#get chr, start, end, ref, alt for all the mutations that need retrival 
ssm_files$combo = paste(ssm_files$gene, ssm_files$patient, sep="_")
retrieve = as.data.table(filter(ssm_files, gene %in% need_to_retrieve$V1))
#retrieve = unique(retrieve[,c("CHROM", "POS", "POS", "REF", "ALT")])
#write.table(retrieve, file="mutations_to_retrieve_from_bam_files.bed", quote=F, row.names=F, sep="\t")

#remove mutations that aren't present in everyone from original file
#first add tag to seperate patients that had mutation versus those that didn't 
#ssm_files = as.data.table(filter(ssm_files, !(gene %in% need_to_retrieve$V1)))
ssm_files$tag = ""
z = which(ssm_files$combo %in% retrieve$combo)
ssm_files$tag[z] = "somatic_this_patient_recount_others"
ssm_files$tag[-z] = "somatic_across_all"

dim(ssm_files)

##results
#find -L . -name "*missing_counts.bed" > bam_readcount_results #20
files = (fread("bam_readcount_results", header=F))
files = files$V1

#in R 
library(params)
library(readr)
source("/cluster/home/kisaev/scripts/bam_readcount_parseline.R")

get_dat = function(file){
  #readcount output 
  x = file 
  #mutation data 
  bed_f = "mutations_to_retrieve_from_bam_files.bed"
  bed = fread(bed_f)
  colnames(bed) = c("chr", "start", "end", 
                                     "ref_allele", "alt_allele")

  write.table(bed, bed_f, quote=F, row.names=F, sep="\t")
  bed = bed_f
  samplename = unlist(strsplit(file, "/"))[6]
  #get output 
  output_reads_for_muts = bam_readcount.parse(x, samplename = samplename, bed)
  output_reads_for_muts = as.data.table(output_reads_for_muts)

  return(output_reads_for_muts)
}

missing_variants_data1 = llply(files, get_dat, .progress="text")
missing_variants_data = ldply(missing_variants_data1)
missing_variants_data$samplename = sapply(missing_variants_data$samplename, function(x){paste(unlist(strsplit(x, "_"))[1:6], collapse="_")})
missing_variants_data = as.data.table(missing_variants_data)

#now need to convert this into PhyloWGS format 
head(ssm_files)
head(missing_variants_data) 

#which of these patients actually had a somatic mutation
#colnames needed = chr, position, description = gene, profile = list of 0/1, VAF for each patient 
#combine somatic with bam recounts 
ssm_files$vaf = (as.numeric(ssm_files$d)-as.numeric(ssm_files$a)) / as.numeric(ssm_files$d)
ssm_files = ssm_files[,c("CHROM", "POS", "gene", "vaf", "patient", "tag")]

missing_variants_data$tag = "bam_readcounts"
missing_variants_data$vaf = missing_variants_data$alt_count / (missing_variants_data$alt_count + missing_variants_data$ref_count)
missing_variants_data$gene = paste(missing_variants_data$chr, missing_variants_data$start, sep="_")
#keep only variants in ssm_files
missing_variants_data = as.data.table(filter(missing_variants_data, gene %in% ssm_files$gene))
missing_variants_data = missing_variants_data[,c("chr", "start", "gene", "vaf", "samplename", "tag")]
colnames(missing_variants_data) = colnames(ssm_files)

all_muts = rbind(missing_variants_data, ssm_files)
unique_muts = unique(all_muts$gene) #33120 unique muts 

make_input = function(mut){
  mut_dat = as.data.table(filter(all_muts, gene == mut))
  dup=mut_dat[which(duplicated(mut_dat$patient))]$patient
  if(!(length(dup)==0)){
    z = which((mut_dat$patient %in% dup) & (mut_dat$tag == "bam_readcounts"))
    mut_dat = mut_dat[-z,]}
  if(dim(mut_dat)[1] == 20){
  mut_dat = mut_dat[order(patient)]
  mut_dat$binary[mut_dat$tag == "bam_readcounts"] = 0
  mut_dat$binary[mut_dat$tag == "somatic_this_patient_recount_others"] = 1
  description = paste(mut_dat$binary, collapse="")
  mut_dat$profile = description
  mut_dat$Normal = 0
  mut_dat_f = as.data.table(dcast(mut_dat, CHROM+ POS+gene + profile + Normal ~ patient, value.var = "vaf"))
  return(mut_dat_f)
  }
}

liche_input = as.data.table(ldply(llply(unique_muts, make_input, .progress="text")))
colnames(liche_input)[1:4] = c("#chr" , "position" , "description" , "profile") 

#some were not available through bam readcount, either remove or label as 0











