#summarize number of mutations explained by hatchet CNA clones

library(data.table)
library(ggpubr)
library(dplyr)
library(plyr)


setwd("/Users/kisaev/Documents/Hatchet_analysis")

#1. all muts
muts = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files/SNVs/Mutect2_Strelka_merged_mutations_wCNA_status/2021-04-06_Mutect2_Strelka_merged_mutations_wCNA_status.rds")

#2. mutations vs cnas from hatchet
p1=fread("p001_SNVs_explained_by_hatchet_output.txt")
p1$patient = "LY_RAP_0001"

p2=fread("p002_SNVs_explained_by_hatchet_output.txt")
p2$patient = "LY_RAP_0002"

#p3=fread("p003_SNVs_explained_by_hatchet_output.txt")

#all_hatchet = rbind(p1, p2, p3)
all_hatchet = rbind(p1, p2)

snvs_explained = as.data.table(table(all_hatchet$Explained, all_hatchet$'PATIENT-SAMPLE', all_hatchet$patient)) %>% filter(N >0)
colnames(snvs_explained) = c("Explained", "Sample", "Patient", "Num_mutations")

#get percentage
tot_mutations = as.data.table(table(all_hatchet$'PATIENT-SAMPLE', all_hatchet$patient)) %>% filter(N >0)
colnames(tot_mutations) = c("Sample", "Patient", "Total_mutations")

snvs_explained = merge(snvs_explained, tot_mutations)
snvs_explained$perc_explained = snvs_explained$Num_mutations / snvs_explained$Total_mutations
snvs_explained$sample = sapply(snvs_explained$Sample, function(x){unlist(strsplit(x, "-"))[2]})

#plot

p1 = ggbarplot(filter(snvs_explained, Patient == "LY_RAP_0001"), x="sample", y="perc_explained", fill="Explained", palette=c("red", "blue")) +
 coord_flip() + ylab("% of mutations")
p1
ggsave("P001_percentage_mutations_explained_by_allele_specific_copy_number_events.pdf", width=4, height=3)

p2 = ggbarplot(filter(snvs_explained, Patient == "LY_RAP_0002"), x="sample", y="perc_explained", fill="Explained", palette=c("red", "blue")) +
 coord_flip() + ylab("% of mutations")
p2
ggsave("P002_percentage_mutations_explained_by_allele_specific_copy_number_events.pdf", width=4, height=3)
