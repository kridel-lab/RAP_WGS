#compare ctDNA mutations to bulk DNA sequencing from tumour samples

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Load libraries and data from RAP samples
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R") #load all mutations
library("gplots")
library(threadr)
library(GenomicRanges)
library(VennDiagram)
library(UpSetR)

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
ctDNA_muts=filter(ctDNA_muts, FILTER=="PASS", DP > 100)
ctDNA_muts$mut_id = paste(ctDNA_muts$Chromosome, ctDNA_muts$Start_Position, sep="_") #949 entries

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Find overlap of mutations across Consensus Cruncher corrections
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#sscs_NA is just sscs and dcs_NA is just dcs

#make venn diagram
patients = c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")

make_venn = function(patient, mut_data, type_analysis){

  print(type_analysis)

  #make folder to store plots for patient
  mainDir <- getwd()
  subDir <- paste("mutations_overlap_corrections", date, patient, sep="_")

  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))

  sscs = unique(filter(mut_data, STUDY_PATIENT_ID == patient, correction_type == "sscs_NA")$mut_id)
  sscs_sc = unique(filter(mut_data, STUDY_PATIENT_ID == patient, correction_type == "sscs_sc")$mut_id)
  dcs = unique(filter(mut_data, STUDY_PATIENT_ID == patient, correction_type == "dcs_NA")$mut_id)
  dcs_sc = unique(filter(mut_data, STUDY_PATIENT_ID == patient, correction_type == "dcs_sc")$mut_id)

  #myCol
  mycol = c("#B3E2CD" ,"#FDCDAC" ,"#CBD5E8" ,"#F4CAE4")

  if(!(length(dcs_sc)==0)){

  # Chart
  venn.diagram(
  x = list(sscs, sscs_sc, dcs, dcs_sc),
  category.names = c("sscs" , "sscs_sc", "dcs", "dcs_sc"),
  filename = paste(patient, type_analysis, "venn_diagram_consensus_cruncher_mutations.png", sep="_"),
  output=TRUE,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycol)

  #make barplot summary of venn diagram
  listInput <- list(sscs = sscs, sscs_sc = sscs_sc,
    dcs = dcs, dcs_sc=dcs_sc)

  pdf(paste(patient, type_analysis, "venn_diagram_barplot_consensus_cruncher_mutations.pdf", sep="_"))
  print(upset(fromList(listInput), order.by = "freq"))
  dev.off()

  print("done plots")}

  setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls")
}

llply(patients, make_venn, ctDNA_muts, "ctDNA_only")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#RAP mutations
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

all_muts = unique(read_only[,c("STUDY_PATIENT_ID", "Indiv", "symbol", "chr", "POS", "REF", "ALT", "gt_AF", "DP")])
colnames(all_muts) = c("STUDY_PATIENT_ID", "Indiv", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "RAP_AF", "RAP_DP")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Merge mutations
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

merged_muts = merge(ctDNA_muts, all_muts, by = c("STUDY_PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"))
write.csv(merged_muts, paste(date, "All_samples_mutations_found_in_ctDNA_and_autopsy_samples.csv", sep="_"), quote=F, row.names=F)

#get venn diagram for how many mutations of those merged appeared in each correction
llply(patients, make_venn, merged_muts, "merged_with_RAP")

#save coordinates of regions evaluated in the end
write.csv(all_targets, paste(date, "Final_regions_evaluated_mean_targeted_coordinates.csv", sep="_"), quote=F, row.names=F)
