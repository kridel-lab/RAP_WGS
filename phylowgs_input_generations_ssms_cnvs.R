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

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution")

date = Sys.Date()

#what output should look like ----> below

#c0(id) 66023,50883,62757,36056,58777 (a in order) 126755,100469,121941,71263,115417 (d in order) s2,1,2;s4,0,1 (ssms in the region, minor, major copy numnber)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#tried using parser provided by PhyloWGS... did not work 
#with our VCF files, eventually worth debugging it and comparing results
#will generate CNV inout file using this script instead to match
#required input

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#-----------#
#parsed cnvs#
#-----------#

#parsed cnv files
parsed = list.files(pattern="segs.txt.parsed.txt")

#get full parsed cnv info - combine
get_cnv_file=function(file){
  f = fread(file)
  samp = paste(unlist(strsplit(file, "_"))[1:3], collapse="_")
  f$sample = samp
  return(f)
}
all_parsed_cnvs = as.data.table(ldply(llply(parsed, get_cnv_file)))

#minor modifs to parsed file to intersect with SSM 
parsed_cords = all_parsed_cnvs[,c("chromosome", "start", "end", "sample", "cellular_prevalence", "copy_number", "minor_cn", "major_cn")]
write.table(parsed_cords, file="cnv_cords_to_intersect_with_ssms.txt", sep="\t", quote=F, row.names=F, col.names=F)

#-----------#
#ssms info. #
#-----------#

#ssm file
ssm = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/2019-07-16_ssm_data_RAP_WGS_input.txt")
order_pats = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/2019-07-16_ssm_data_RAP_WGS_input_order.txt")

#minor modifs to ssm file
ssm$id = paste("s", 0:(nrow(ssm)-1), sep="_")
ssm$chr = sapply(ssm$gene, function(x){unlist(strsplit(x, "_"))[1]})
ssm$start = sapply(ssm$gene, function(x){unlist(strsplit(x, "_"))[2]})
ssm_cords = ssm[,c("chr", "start", "start", "id")]
colnames(ssm_cords)[3] = "end"
write.table(ssm_cords, file="ssm_cords_to_intersect_with_cnvs.txt", sep="\t", quote=F, row.names=F, col.names=F)

#-----------#
#sample ids #
#-----------#

#sample to patient id conversion 
conv = fread("2019-08-06-optimalClusterSolution.txt")
colnames(conv)[12] = "sample_id"
conv = conv[,c("barcode", "sample_id")]
colnames(conv)[1] = "sample"

##########################
#in bedtools - terminal. # 
##########################

#module load bedtools
#bedtools intersect -a cnv_cords_to_intersect_with_ssms.txt -b ssm_cords_to_intersect_with_cnvs.txt -wa -wb > ssms_in_cnvs.txt

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

ssms_in_cnvs = fread("ssms_in_cnvs.txt")
colnames(ssms_in_cnvs) = c("cnv_chr", "cnv_start", "cnv_end", "sample", "cellular_prevalence", "copy_number", "minor_cn", "major_cn", "ssm_chr", "ssm_start", "ssm_end", "ssm_id")

#add copy number info and depth info and sample id 
ssms_in_cnvs = merge(ssms_in_cnvs, conv, by="sample")

#for each one get a and d values
#average read depth = 60 (estimate for now)

ssms_in_cnvs$a = ""
ssms_in_cnvs$d = ""

get_a_d = function(cnv_end, cnv_start, cellular_prevalence){
  length_cna = (cnv_end-cnv_start)
  d = 60 * min(3000, ((7 * length_cna)/(10^4)))
  a = d * (1-(cellular_prevalence/2))
  print(c(a,d))
  return(c(a,d))
}

test = mapply(get_a_d, ssms_in_cnvs$cnv_end, ssms_in_cnvs$cnv_start, ssms_in_cnvs$cellular_prevalence)
test = t(test)
test = as.data.table(test)
colnames(test) = c("a", "d")

ssms_in_cnvs$a = test$a
ssms_in_cnvs$d = test$d
ssms_in_cnvs$cna = paste(ssms_in_cnvs$cnv_chr, ssms_in_cnvs$cnv_start, ssms_in_cnvs$cnv_end, sep="_")
length(unique(ssms_in_cnvs$cna))

#function that merged ssms per cna 
unique_cnas = unique(ssms_in_cnvs$cna)

get_list_ssms = function(ssm){
  res = paste(ssm, filter(dat, ssm_id == ssm)$major_cn, filter(dat, ssm_id == ssm)$minor_cn, sep=",")
  return(res)
}

get_cna_ssms = function(cna_id){
    dat = as.data.table(filter(ssms_in_cnvs, cna == cna_id))
    num_pats = length(unique(dat$sample_id))

    #cnv: identifier for each CNV. Identifiers must start at c0 and increment, so the first data row will have c0, the second row c1, and so forth.
    #a: number of reference reads covering the CNV.
    #d: total number of reads covering the CNV. This will be affected by factors such as total copy number at the locus, sequencing depth, and the size of the chromosomal region spanned by the CNV.
    #ssms: SSMs that overlap with this CNV. Each entry is a comma-separated triplet consisting of SSM ID, maternal copy number, and paternal copy number. These triplets are separated by semicolons.

    #multiple patients have cna
    if(num_pats == 1){
      cnv = cna
      a = dat$a[1]
      d = dat$d[1]
      ssms=paste(sapply(dat$ssm_id, get_list_ssms), collapse=";")
    }

    if((num_pats >1) & (num_pats < 20)){
      cna_status = unique(dat$major_cn)
      #are they the same copy number?
      if(length(unique(cna_status)==1)){

      }
      
      #are they diff copy number?
      if(length(unique(cna_status)==1)){
        #if diff need to split into multiple lines 
        
      }

    }

    #only one patient has cna 
    if(num_pats ==20){

    }

}







