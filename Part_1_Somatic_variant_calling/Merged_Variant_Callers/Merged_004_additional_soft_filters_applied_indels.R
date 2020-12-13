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
library("readxl")
library(GenomicRanges)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#here:
#1. remove blacklisted region mutations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

files = list.files(pattern="vcf_file_filtered_indels.bed")

#sample info
samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")

#dlbcl panel
panel = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/CCG.Lymphoma.DNA.targets.GRCh37.v1.0.annotation.xlsx"))

#++++++++++++++++++++
#oicr panel
#++++++++++++++++++++

#-->full genes
oicr_full = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Capture_Panel_ctDNA_Lymphoma_entire_genes.csv")

#-->partially overlapping genes
oicr_partial = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Capture_Panel_ctDNA_Lymphoma_partial_genes.xlsx"))

#read all files and rbind
all_muts = as.data.table(ldply(llply(files, function(x){fread(x)})))
all_muts_cords = unique(all_muts[,c("CHROM", "POS", "mut_id")])
all_muts_cords$POSend = all_muts_cords$POS
all_muts_cords$CHROM = paste("chr", all_muts_cords$CHROM, sep="")
all_muts_cords = all_muts_cords[,c("CHROM", "POS", "POSend", "mut_id")]
all_muts_cords=all_muts_cords[order(CHROM, POS)]
write.table(all_muts_cords, file="all_muts_pre_blacklist_filter_indels.bed", quote=F, row.names=F, sep="\t", col.names=F)

#DO NOT RUN
#module load bedtools
#bedtools intersect -wa -wb \
#    -a all_muts_pre_blacklist_filter_indels.bed \
#    -b hg19-blacklist.v2.bed -f 1.0 > all_muts_post_blacklist_filter_indels.bed \

blacklist_muts = fread("all_muts_post_blacklist_filter_indels.bed")
colnames(blacklist_muts)[1:4] = c("CHROM", "POS", "POS_end", "mut_id")
z = which(blacklist_muts$mut_id %in% all_muts$mut_id) #1185 mutations
z = which(all_muts$mut_id %in% blacklist_muts$mut_id) #3227 entries mutations
all_muts = all_muts[-z,] #remove mutations in blacklisted regions

#other filters to apply to these variants?
#split into num ref and alt allele counts
all_muts = all_muts %>% separate(gt_AD, c("Ref_counts", "alt_counts"))
all_muts$Ref_counts = as.numeric(all_muts$Ref_counts)
all_muts$alt_counts = as.numeric(all_muts$alt_counts)
#remove empty columns
all_muts$gt_DP = NULL
all_muts$gt_GQ = NULL
all_muts$gt_PGT = NULL
all_muts$gt_PID = NULL
all_muts$gt_PL = NULL
all_muts$ID = NULL
all_muts$QUAL = NULL
all_muts$AC = NULL
all_muts$AF = NULL
all_muts$AN = NULL
all_muts$IN_PON = NULL
all_muts$NLOD = NULL
all_muts$N_ART_LOD = NULL
all_muts$RPA = NULL
all_muts$RU = NULL
all_muts$STR = NULL
all_muts$ALLELE_END = NULL
length(unique(all_muts$mut_id)) #146,690

#TLOD - min 10 - confidence score for variant being really somatic
all_muts = as.data.table(filter(all_muts, TLOD > 10))
length(unique(all_muts$mut_id)) #145,617 across three patients

#min counts for alternative reads
all_muts = as.data.table(filter(all_muts, alt_counts > 5))
length(unique(all_muts$mut_id)) #145,499

#save file so far --> download locally and plot summary stats for various variable
saveRDS(all_muts, file=paste(date, "all_muts_merged_mut_calls_indels.rds", sep="_")) #this file is used for downstream exploratory analysis

colnames(samps)[which(colnames(samps)=="LY_RAP_ID")] = "Indiv"

#save clean version just columns that are important
clean = all_muts[,c("Indiv", "CHROM", "POS", "mut_id",
"REF", "ALT", "Ref_counts", "alt_counts", "DP", "ensgene",
"Func.ensGene", "ExonicFunc.ensGene", "AAChange.ensGene",
"cosmic68", "symbol", "chr", "start", "end", "biotype")]

#add sample info on where it came from
clean = merge(clean, samps, by = "Indiv")
write.table(clean, file=paste(date, "all_indels_samples.txt", sep="_"), quote=F,
row.names=F, sep="\t")

#check which genes from lympohoma panel are mutated
#remove noncoding mutations - intronic, intergenic,
clean = as.data.table(filter(clean, !(Func.ensGene %in% c("intronic", "intergenic"))))
panel$rap_wgs_muts = ""
z = which(panel$gene %in% clean$symbol)
genes_w_mut = panel$gene[z]
panel$rap_wgs_muts[z] = "at_least_one_sample_wmut_in_gene"
#write.table(panel, file=paste(date, "CGC_lymphoma_targets_genes_wmuts_in_RAP_WGS.txt", sep="_"), quote=F, row.names=F, sep=";")

#check which genes from OICR lymphoma panel
#full genes vs partially overlapping genes
#which genes mutated that are captured in full?
oicr_full$rap_wgs_muts = ""
z = which(oicr_full$GENE %in% clean$symbol)
genes_w_mut = oicr_full$GENE[z]
oicr_full$rap_wgs_muts[z] = "at_least_one_sample_wmut_in_gene"
oicr_full$coverage = "full"

#get intersection of mutations with oicr partial gene panel
clean$START = clean$POS
clean$END = clean$POS
clean$CHR = paste("chr", clean$CHROM, sep = "")

all_muts_clean = clean[,c("CHR", "START", "END", "DNA", "symbol")]
all_muts_clean$DNA = as.character(all_muts_clean$DNA)
all_muts_clean$DNA = "*"
all_muts_clean_gr = makeGRangesFromDataFrame(all_muts_clean, keep.extra.columns=TRUE)

#partial oicr panel make granges object
oicr_partial$STRAND[which(oicr_partial$STRAND == "NA")] = "*"
oicr_p = oicr_partial[,c("CHR", "START", "END", "STRAND", "GENE")]
oicr_p_gr = makeGRangesFromDataFrame(oicr_p, keep.extra.columns=TRUE)

hits <- findOverlaps(oicr_p_gr, all_muts_clean_gr, ignore.strand=TRUE)
hits_overlap = cbind(oicr_p[queryHits(hits),], all_muts_clean[subjectHits(hits),])
colnames(hits_overlap) = c("panel_CHR", "panel_START", "panel_END", "panel_STRAND",
"panel_GENE", "RAP_mut_CHR", "RAP_mut_START", "RAP_mut_END", "DNA", "RAP_mut_gene")
hits_overlap$DNA = NULL
hits_overlap = unique(hits_overlap)
hits_overlap$REGION = paste(hits_overlap$panel_CHR, hits_overlap$panel_START,
hits_overlap$panel_END, hits_overlap$panel_GENE, sep="_")
oicr_partial$REGION = paste(oicr_partial$CHR, oicr_partial$START,
oicr_partial$END, oicr_partial$GENE, sep="_")

#summarize results
oicr_full$REGION = "FULL_GENE"
oicr_part = unique(oicr_partial[,c("GENE", "REGION", "COMMENT")])
z = which(oicr_part$REGION %in% hits_overlap$REGION)
oicr_part$rap_wgs_muts=""
oicr_part$rap_wgs_muts[z] = "at_least_one_sample_wmut_in_region"
oicr_full$COMMENT = ""
oicr_full = oicr_full[,c("GENE", "REGION", "rap_wgs_muts", "COMMENT")]
oicr_part = oicr_part[,c("GENE", "REGION", "rap_wgs_muts", "COMMENT")]
all_oicr = rbind(oicr_part, oicr_full)
#write.table(all_oicr, file=paste(date, "OICR_lymphoma_panel_summary_RAP_0003_patient.txt", sep="_"),
#quote=F, row.names=F, sep=";")
