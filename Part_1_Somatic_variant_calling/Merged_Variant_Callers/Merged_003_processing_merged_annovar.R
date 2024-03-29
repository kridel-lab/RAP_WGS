#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

date = Sys.Date()
print(date)

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)
library(VariantAnnotation)
#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#library(BSgenome.Hsapiens.UCSC.hg19)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcfs_annovar_annotated")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

paired = list.files(pattern="no_info_AF")
z = which(str_detect(paired, "indel"))
paired = paired[-z]
print(length(paired))

#gene annotations
genes = unique(fread("/cluster/home/kisaev/data/annotables_grch37.txt"))
genes = unique(genes[,c("ensgene", "symbol", "chr", "start", "end",
"strand", "biotype")])

clean_up_001 = function(paired_vcf){

  mutations_T1 =read.vcfR(paired_vcf)
  mutations_T1 = vcfR2tidy(mutations_T1)
  meta = as.data.table(mutations_T1$meta)

  vcf = as.data.table(mutations_T1$fix)
  chr_conv = unique(vcf[,c("ChromKey", "CHROM")])

  gt = as.data.table(mutations_T1$gt)
  gt = merge(gt, chr_conv, by="ChromKey")

  gt$mut_id = paste(gt$CHROM, gt$POS, sep="_")
  vcf$mut_id = paste(vcf$CHROM, vcf$POS, sep="_")

  #1. keep only the ones that passed default hard filter
  vcf = as.data.table(filter(vcf, FILTER=="PASS"))
  print(paste("number of variants that passed filtering=", dim(vcf)[1]))

  #2. filter by coverage
  #try 60X for all tumour samples except for ffpe diagnostic ones
	if(str_detect(paired_vcf, "Dia_")){
		vcf = as.data.table(filter(vcf, DP >=30))
	  print(paste("number of variants that passed coverage=", dim(vcf)[1]))
	}

	if(!(str_detect(paired_vcf, "Dia_"))){
		vcf = as.data.table(filter(vcf, DP >=60))
	 print(paste("number of variants that passed coverage=", dim(vcf)[1]))
	}

  #3. combine vcf and gt info
  cols = colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt = merge(gt, vcf, by= cols)

  #split file into two - tumour and control
  split <- split(gt , f = gt$Indiv)

  #keep only tum sample
  z = which(str_detect(gt$Indiv, "Ctl"))
  gt = gt[-z,]

  #6. keep only those with population allele frequency < 0.001
  gt$controls_AF_popmax = as.numeric(gt$controls_AF_popmax)
  gt = as.data.table(filter(gt, (controls_AF_popmax < 0.001 | is.na(controls_AF_popmax))))
  print(paste("number of variants that passed controls_AF_popmax=", dim(gt)[1]))

	#4. remove potential snps annotated by gnomad - these could include
	#known disease causing variants
  z = which((str_detect(gt$avsnp142, "rs")) & (gt$cosmic68 == "."))
  if(!(length(z)==0)){
  gt = gt[-z,]}
  print(paste("number of variants that passed avsnp142=", dim(gt)[1]))

  #7. keep only mutations where t2 VAF > 0.1
  gt$gt_AF = as.numeric(gt$gt_AF)
  gt = as.data.table(filter(gt, gt_AF >= 0.1))
  print(paste("number of variants that passed vaf >= 0.1=", dim(gt)[1]))

  #8. remove variants from chromosome X and Y
  gt = as.data.table(filter(gt, !(CHROM %in% c("X", "Y"))))
  print(paste("number of variants that passed X Y=", dim(gt)[1]))

  #9. add gene name
  #if multiple genes mapped to variant (usually if in between genes or upstream of genes)
  #keep id of first gene
  gt$ensgene = sapply(gt$Gene.ensGene, function(t){unlist(strsplit(t, '\\x', fixed=TRUE))[1]})
  gt = merge(gt, genes, by="ensgene")

  #11. generate bed file - summary of mutation and coordinates to intersect with cnvkit output
  pat = gt$Indiv[1]
  file = paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text/", pat, "_final_vcf_file_filtered.bed", collapse="_", sep="")
  write.table(gt, file, quote=F, row.names=F, sep="\t", col.names=T)
  print("done")
}

llply(paired, clean_up_001)
