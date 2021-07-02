#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)
library("ggsci")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#dir.create(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))
#setwd(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))

svs = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Manta/2021-04-20_RAP_WGS_all_SVs_heatmap_plot.rds")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#average number of SVs per sample/patient
t=as.data.table(table(svs$STUDY_PATIENT_ID, svs$Indiv))
t=filter(t, N >0)
t %>% group_by(V1) %>% dplyr::summarize(mean = mean(N))
t %>% group_by(V1) %>% dplyr::summarize(sd = sd(N))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

get_SVs = function(patient){

  #1. SNPs keep
  print(paste(patient, "started"))

  #keep only ancestral ones

  if(patient == "LY_RAP_0001"){
    patient_name="MCL blastoid stage IV"
  }

  if(patient == "LY_RAP_0002"){
    patient_name = "PMBCL stage IV bulky B symptoms"
  }

  if(patient == "LY_RAP_0003"){
    patient_name = "DLCBL double hit stage IV"
  }

  svs_keep = filter(svs, STUDY_PATIENT_ID == patient)

  bnds = filter(svs_keep, SVTYPE == "BND")
  not_bnds = filter(svs_keep, !(SVTYPE == "BND"))

  bnds$chr1 = sapply(bnds$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[2], "_"))[1]})
  bnds$chr2 = sapply(bnds$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[3], "_"))[1]})

  z = which(bnds$chr1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
"13", "14", "15", "16", "17", "18", "19", "20", "X", "Y"))
  bnds = bnds[z,]
  z = which(bnds$chr2 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
"13", "14", "15", "16", "17", "18", "19", "20", "X", "Y"))
  bnds = bnds[z,]

  bnds$gene1 = sapply(bnds$gene, function(x){unlist(strsplit(as.character(x), "_"))[1]})
  bnds$gene2 = sapply(bnds$gene, function(x){unlist(strsplit(as.character(x), "_"))[2]})
  bnds$new_id = paste(bnds$chr1, bnds$chr2, bnds$gene2, sep="_")

  z = which(bnds$gene2 == "NA")
  bnds$new_id[z] = paste(bnds$chr1[z], bnds$chr2[z], bnds$gene1[z], sep="_")

  not_bnds$new_id = paste(not_bnds$gene, not_bnds$SVTYPE, sep="_")

  bnds = bnds %>% select("Indiv", "Tissue_Site", "STUDY_PATIENT_ID", "SVTYPE", "gene", "id", "new_id")
  svs_keep = rbind(bnds, not_bnds)
  svs_keep = unique(svs_keep[,c("Indiv", "new_id")])
  t=as.data.table(table(svs_keep$new_id))
  colnames(t) = c("SV", "N")
  t = as.data.table(table(t$N))
  colnames(t) = c("num_of_samples_with_mut", "num_of_muts")
  t$Patient = patient_name
  t = t %>% select("num_of_samples_with_mut", "Patient", "num_of_muts")
  return(t)
}

patients = c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")
allSVs = as.data.table(ldply(llply(patients, get_SVs)))

saveRDS(allSVs, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Figure2_input_data_sample_dist_SVs.rds")

#pdf("001_samples_per_mutation_lineplot_SVs.pdf", width=4, height=5)
# Basic barplot
#p<-ggline(barplot, x="num_of_samples_with_mut", y="num_of_muts",
#palette = c("#00AFBB", "#E7B800", "#FC4E07"), color="Patient")+
#xlab("# of samples sharing mutation") + ylab("Structural variants count")+
#scale_y_continuous(breaks=seq(0,1000,100))
#print(p)
#dev.off()
