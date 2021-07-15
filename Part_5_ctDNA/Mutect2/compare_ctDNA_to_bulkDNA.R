#compare ctDNA mutations to bulk DNA sequencing from tumour samples

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Load libraries and data from RAP samples
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
#library(clonevol)
library("gplots")
library(threadr)
library(GenomicRanges)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ctDNA data load
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ctDNA = list.files(pattern="_Mutect2_annovar_mutations_all.rds")[length(list.files(pattern="_Mutect2_annovar_mutations_all.rds"))]
ctDNA = readRDS(ctDNA)

ctDNA$STUDY_PATIENT_ID = ""
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0001"] = "LY_RAP_0001"
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0002"] = "LY_RAP_0002"
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0003"] = "LY_RAP_0003"

ctDNA_muts = unique(ctDNA[,c("STUDY_PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position",
"Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Var_Freq",
"FILTER", "Tumor_Sample_Barcode", "correction_type", "DP")])
ctDNA_muts = filter(ctDNA_muts, !(Chromosome == "X"))
ctDNA_muts$Chromosome = as.numeric(ctDNA_muts$Chromosome)
ctDNA_muts=filter(ctDNA_muts, FILTER=="PASS")
ctDNA_muts$mut_id = paste(ctDNA_muts$Chromosome, ctDNA_muts$Start_Position, sep="_")

#targets list for 46 genes
targets_pcg = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/my-targets.bed")
targets_pcg$V1 = sapply(targets_pcg$V1, function(x){unlist(strsplit(x, "chr"))[2]})
targets_pcg = unique(targets_pcg[,c(1:4)]) #571 targets
colnames(targets_pcg) = c("Chr", "Start", "Stop", "Target")
targets_pcg$type = "46_genes"

#targets list for additional 152 regions
target_regs = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/NGS-Targets.xlsx"))
target_regs = target_regs[,c("Chr", "Start", "Stop", "Target")] #1,675 baits including 38 which correspond to the Agena sample identity probes (TargetNN)
target_regs$type = "other_regions"
all_targets = rbind(targets_pcg, target_regs)

colnames(all_targets)[1:5]=c("chr", "target_start", "target_stop", "target_gene", "type")

#RAP mutations
all_muts = unique(read_only[,c("STUDY_PATIENT_ID", "Indiv", "symbol", "chr", "POS", "REF", "ALT", "gt_AF", "DP")])
colnames(all_muts) = c("STUDY_PATIENT_ID", "Indiv", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "RAP_AF", "RAP_DP")

#Merge mutations
merged_muts = merge(ctDNA_muts, all_muts, by = c("STUDY_PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
write.csv(merged_muts, paste(date, "All_samples_mutations_found_in_ctDNA_and_autopsy_samples.csv", sep="_"), quote=F, row.names=F)

#save targeted regions coverage summary (only those > 100)
#write.csv(cov_sum, paste(date, "Final_regions_evaluated_mean_targeted_coverage.csv", sep="_"), quote=F, row.names=F)

#save coordinates of regions evaluated in the end
#all_targets = filter(all_targets, target_gene %in% cov_sum$target_gene)
write.csv(all_targets, paste(date, "Final_regions_evaluated_mean_targeted_coordinates.csv", sep="_"), quote=F, row.names=F)
