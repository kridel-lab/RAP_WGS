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

#data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")
cnas = unique(cnas %>% select("Sample", "Ploidy", "Purity", "Tissue_Site"))

private_muts = filter(samples_per_mut, phylogeny == "private")

private_muts_full = unique(filter(read_only, mut_id %in% private_muts$mut_id) %>%
              mutate(patient = STUDY_PATIENT_ID) %>%
              select("patient", "Tissue_Site", "Sample", "mut_id"))

private_muts = merge(private_muts, private_muts_full, by=c("patient", "mut_id"))

#summary number of private mutations per sample per patient
summary_private = as.data.table(table(private_muts$patient, private_muts$Sample)) %>% filter(N >0)
colnames(summary_private) = c("Patient", "Sample", "N")

summary_private = merge(summary_private, cnas, by="Sample")

#get total number of mutations per sample
t = as.data.table(table(read_only$Sample)) ; colnames(t) = c("Sample", "total_mutations")
summary_private = merge(summary_private, t, by="Sample")
summary_private$private_percent = summary_private$N/summary_private$total_mutations * 100

#plot

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/001_private_mutations_vs_purity.pdf",
width=5, height=5)
p = ggscatter(summary_private, x = "Purity", y = "private_percent",
          palette = "jco",           # Color by groups "cyl"
          color = "Patient" , shape=21,                             # Change point shape by groups "cyl"
        )+xlim(0.4, 1) + ylim(0, 100)+
  stat_cor(method = "spearman")+ylab("% private mutations")
print(p)
dev.off()
