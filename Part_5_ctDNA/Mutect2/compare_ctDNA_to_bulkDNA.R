#compare ctDNA mutations to bulk DNA sequencing from tumour samples

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

#Data++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ctDNA = list.files(pattern="_Mutect2_annovar_mutations_all.rds")[length(list.files(pattern="_Mutect2_annovar_mutations_all.rds"))]
ctDNA = readRDS(ctDNA)

ctDNA$STUDY_PATIENT_ID = ""
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0001"] = "LY_RAP_0001"
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0002"] = "LY_RAP_0002"
ctDNA$STUDY_PATIENT_ID[ctDNA$Tumor_Sample_Barcode == "LY_0003"] = "LY_RAP_0003"

ctDNA_muts = unique(ctDNA[,c("STUDY_PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Var_Freq", "FILTER", "Tumor_Sample_Barcode")])
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

cov_sum = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls/2021-07-12_target_based_coverage_summary.csv")
cov_sum = filter(cov_sum, mean_cov > 100) #1400/1905 records

#for each patient figure out regions that passed and intersect with mutation calls
get_muts_targets = function(patient){

  cov_sum_pat = filter(cov_sum, Library == patient)

  #keep only ctDNA mutations in targeted regions with enough coverage
  pat_targs = filter(all_targets, target_gene %in% cov_sum_pat$target_gene)
  pat_targs_gr = makeGRangesFromDataFrame(pat_targs, keep.extra.columns=TRUE)

  #get mutations coordinates
  ctdna_mut_pat = filter(ctDNA_muts, Tumor_Sample_Barcode == patient)
  ctdna_mut_pat = ctdna_mut_pat[,c("Chromosome", "Start_Position", "Reference_Allele", "mut_id")]
  ctdna_mut_pat$Reference_Allele = ctdna_mut_pat$Start_Position
  colnames(ctdna_mut_pat) = c("chr", "start", "end", "mut_id")
  ctdna_mut_pat_gr = makeGRangesFromDataFrame(ctdna_mut_pat, keep.extra.columns=TRUE)
  print(paste("num mutations originally = ", dim(ctdna_mut_pat)[1]))

  hits <- findOverlaps(ctdna_mut_pat_gr, pat_targs_gr, ignore.strand=TRUE)
  hits_overlap = cbind(ctdna_mut_pat[queryHits(hits),], pat_targs[subjectHits(hits),])

  mutations_keep = hits_overlap$mut_id
  patient_mutations_keep = filter(ctDNA_muts, Tumor_Sample_Barcode == patient, mut_id %in% mutations_keep)
  print(paste("num mutations after taget filter = ", dim(patient_mutations_keep)[1]))

  return(patient_mutations_keep)
}

#get ctDNA mutations in regions with at least 100x coverage
patients = c("LY_0001", "LY_0002", "LY_0003")
ctDNA_muts_pass = as.data.table(ldply(llply(patients, get_muts_targets)))

#RAP mutations
all_muts = unique(read_only[,c("STUDY_PATIENT_ID", "Indiv", "symbol", "chr", "POS", "REF", "ALT", "gt_AF", "DP")])
colnames(all_muts) = c("STUDY_PATIENT_ID", "Indiv", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "RAP_AF", "RAP_DP")

#Merge mutations
merged_muts = merge(ctDNA_muts_pass, all_muts, by = c("STUDY_PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
write.csv(merged_muts, paste(date, "All_samples_mutations_found_in_ctDNA_and_autopsy_samples.csv", sep="_"), quote=F, row.names=F)

#save targeted regions coverage summary (only those > 100)
write.csv(cov_sum, paste(date, "Final_regions_evaluated_mean_targeted_coverage.csv", sep="_"), quote=F, row.names=F)

#save coordinates of regions evaluated in the end
all_targets = filter(all_targets, target_gene %in% cov_sum$target_gene)
write.csv(all_targets, paste(date, "Final_regions_evaluated_mean_targeted_coordinates.csv", sep="_"), quote=F, row.names=F)
