#compare ctDNA mutations to bulk DNA sequencing from tumour samples

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
#library(clonevol)
library("gplots")
library(threadr)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls")

#Data++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ctDNA = list.files(pattern="_Mutect2_annovar_mutations_all.rds")[length(list.files(pattern="_Mutect2_annovar_mutations_all.rds"))]
ctDNA = readRDS(ctDNA)

ctDNA$STUDY_PATIENT_ID = ""
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0001"] = "LY_RAP_0001"
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0002"] = "LY_RAP_0002"
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0003"] = "LY_RAP_0003"

ctDNA_muts = unique(ctDNA[,c("STUDY_PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Var_Freq", "FILTER")])
ctDNA_muts$Chromosome = as.numeric(ctDNA_muts$Chromosome)

ctDNA_muts=filter(ctDNA_muts, FILTER=="PASS")

all_muts = unique(read_only[,c("STUDY_PATIENT_ID", "Indiv", "symbol", "chr", "POS", "REF", "ALT")])
colnames(all_muts) = c("STUDY_PATIENT_ID", "Indiv", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")

merged_muts = merge(ctDNA_muts, all_muts, by = c("STUDY_PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
#merged_muts = filter(merged_muts, STUDY_PATIENT_ID == "LY_RAP_0002")
merged_muts$mut_id = paste(merged_muts$Chromosome, merged_muts$Start_Position, sep="_")
write.csv(merged_muts, "All_samples_mutations_found_in_ctDNA_and_autopsy_samples.csv", quote=F, row.names=F)
