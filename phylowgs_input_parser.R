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
  f= as.data.table(filter(f, DP >=60))
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

#find mutations that are not present in all samples
#need to fill in the blank --> retrieve reads from BAM files

t = as.data.table(table(ssm_files$gene))
need_to_retrieve = as.data.table(filter(t, N <20))

#get chr, start, end, ref, alt for all the mutations that need retrival 
retrieve = as.data.table(filter(ssm_files, gene %in% need_to_retrieve$V1))
retrieve = unique(retrieve[,c("CHROM", "POS", "POS", "REF", "ALT")])
write.table(retrieve, file="mutations_to_retrieve_from_bam_files.bed", quote=F, row.names=F, sep="\t")






