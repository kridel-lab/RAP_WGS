#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)
library("ggsci")
require(gridExtra)
library(cowplot)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#All SNVs--------------------------------------------------------------

#summarize number of mutations per sample
z = which((read_only$Tissue_Site == "Adrenal gland, NOS") & (read_only$STUDY_PATIENT_ID == "LY_RAP_0003"))
read_only$Tissue_Site[z] = "Adrenal gland"
z = which((read_only$Tissue_Site == "Aorta, ascending, not specified \n\n") & (read_only$STUDY_PATIENT_ID == "LY_RAP_0001"))
read_only$Tissue_Site[z] = "Aorta, ascending"

samples_per_mut = unique(samples_per_mut %>% select("mut_id", "patient", "phylogeny"))
colnames(samples_per_mut)[2] = "STUDY_PATIENT_ID"
read_only = merge(read_only, samples_per_mut, by=c("mut_id", "STUDY_PATIENT_ID"))

patient = "LY_RAP_0002"

make_geom_heatmap = function(patient){
  muts_pat = filter(read_only, STUDY_PATIENT_ID == patient)
  muts_pat = unique(muts_pat %>% select("mut_id", "Sample", "phylogeny", "ExonicFunc.ensGene", "cosmic68", "symbol"))

  if(patient == "LY_RAP_0001"){
    drivers = all_drivers$Gene[all_drivers$type == "MCL"]
  }

  if(patient == "LY_RAP_0002"){
    drivers = all_drivers$Gene[all_drivers$type == "PMBCL"]
  }

  if(patient == "LY_RAP_0003"){
    drivers = all_drivers$Gene[all_drivers$type == "DLBCL"]
  }

  muts_pat$driver_gene = ""
  z = which(muts_pat$symbol %in% drivers)
  if(!(length(z) == 0)){
    muts_pat$driver_gene[z] = "driver"
  }

  muts_pat$cosmic_label = ""
  muts_pat$cosmic_label[which(!(muts_pat$cosmic68 == "."))] = "cosmic"

  order_make = as.data.table(table(muts_pat$mut_id))
  order_make = order_make[order(-N)]
  muts_pat$mut_id = factor(muts_pat$mut_id, levels=order_make$V1)

  order_samples = as.data.table(table(muts_pat$Sample))
  order_samples = order_samples[order(-N)]
  muts_pat$Sample = factor(muts_pat$Sample, levels=order_samples$V1)
  muts_pat$driver_coding = paste(muts_pat$driver_gene, muts_pat$cosmic_label)
  muts_pat$driver_coding[!(muts_pat$driver_coding == "driver cosmic")] = ""
  muts_pat$driver_coding = factor(muts_pat$driver_coding, levels=c("", "driver cosmic"))

#  muts_pat = muts_pat[1:1000,]

  p = ggplot(muts_pat, aes(mut_id, Sample)) +
#    geom_tile(aes(color=driver_coding)) +
    geom_tile() +
     xlab("Mutation") + ylab("Sample")+
     theme(axis.title.x=element_blank(),
       axis.text.x=element_blank(),
       axis.ticks.x=element_blank(), legend.position="right")+
       scale_color_manual(values=c("#999696", "#ba0707"))

  print(p)

}

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_mutations_P1.pdf", width=15, height=3)
make_geom_heatmap("LY_RAP_0001")
dev.off()

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_mutations_P2.pdf", width=15, height=3)
make_geom_heatmap("LY_RAP_0002")
dev.off()

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_mutations_P3.pdf", width=15, height=5)
make_geom_heatmap("LY_RAP_0003")
dev.off()
