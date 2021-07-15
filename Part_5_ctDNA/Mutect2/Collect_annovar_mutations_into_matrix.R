#----------------------------------------------------------------------
#karin isaev
#post mutect2 and annovar soft filtering and matrix prep
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
date=Sys.Date()

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls/Mutect2_VCF_output/annovar")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", "plyr",
	"ggrepel", "stringr", "maftools")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#mutect2 run on tumour only mode
#annotated by annovar

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

genes = fread("/cluster/home/kisaev/data/annotables_grch37.txt")
colnames(genes)[1] = "Gene.ensGene"

files = list.files(pattern="vcf.gz.hg19_multianno.vcf")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#1. clean up individual paired vcf files

clean_up_001 = function(paired_vcf){

  mutations_T1 =read.vcfR(paired_vcf)
  mutations_T1 = vcfR2tidy(mutations_T1)
  meta = as.data.table(mutations_T1$meta)

  vcf = as.data.table(mutations_T1$fix)
  chr_conv = unique(vcf[,c("ChromKey", "CHROM")])

  gt = as.data.table(mutations_T1$gt)
  gt = merge(gt, chr_conv, by="ChromKey")

	z = which(str_detect(gt$Indiv, "Ctl"))
	if(!(length(z)==0)){
		gt=gt[-z,]
	}

  gt$mut_id = paste(gt$CHROM, gt$POS, sep="_")
  vcf$mut_id = paste(vcf$CHROM, vcf$POS, sep="_")

	#1. combine vcf and gt info
  cols = colnames(gt)[which(colnames(gt) %in% colnames(vcf))]
  gt = merge(gt, vcf, by= cols)

	#2. get hugo gene names
  gt = merge(gt, genes, by= "Gene.ensGene")
  print(paste("number of variants that passed gene merge", dim(gt)[1]))

  #3. generate bed file - summary of mutation and coordinates to intersect with cnvkit output
  pat_name = unlist(strsplit(paired_vcf, "_annovar"))[1]
	pat = unlist(strsplit(pat_name, "\\."))[1]

	gt$sample=pat
	gt$correction_type = paste(unlist(strsplit(pat_name, "\\."))[2:3], collapse="_")

	colnames(gt)[which(colnames(gt)=="symbol")] = "Hugo_Symbol"
	gt$End_Position = gt$POS
	gt$Start_Position = gt$POS
	gt$Chromosome = gt$CHROM
	gt$Reference_Allele = gt$REF
	gt$Tumor_Seq_Allele2 = gt$ALT
	gt$Variant_Classification = paste(gt$Func.ensGene, gt$ExonicFunc.ensGene)
	gt$Variant_Type = "SNP"
	gt$Tumor_Sample_Barcode = gt$sample
	gt$Var_Freq = gt$gt_AF

	gt = gt[,c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
	"Tumor_Seq_Allele2", "avsnp142", "cosmic68", "AAChange.ensGene",
	"Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode", "Var_Freq", "FILTER", "correction_type", "DP")]

	return(gt)

  print("done")
}

all_muts = as.data.table(ldply(llply(files, clean_up_001, .progress="text")))

write.table(all_muts, file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls/",
date, "_Mutect2_annovar_mutations_all.txt", sep=""), quote=F, row.names=F, sep=";")

saveRDS(all_muts, file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls/",
date, "_Mutect2_annovar_mutations_all.rds", sep=""))
