#----------------------------------------------------------------------
#processing_annovar_results.R
#karin isaev
#July 11th, 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final")

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

paired = list.files(pattern="final_vcf_file_filtered.bed")

get_mut = function(file){
	f = fread(file)
	return(f)
}

all_muts = ldply(llply(paired, get_mut))
colnames(all_muts) =  c("Gene.ensGene" ,"ChromKey" ,                   
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
all_muts = as.data.table(all_muts)
saveRDS(all_muts, file="final_somatic_mutations_RAP_WGS_20_samples.rds")





