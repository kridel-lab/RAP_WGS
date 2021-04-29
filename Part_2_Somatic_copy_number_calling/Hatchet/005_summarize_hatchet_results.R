#summarize number of mutations explained by hatchet CNA clones

library(data.table)
library(ggpubr)
library(dplyr)

setwd("/Users/kisaev/Documents/Hatchet_analysis")

#1. all muts
muts = readRDS("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files/SNVs/Mutect2_Strelka_merged_mutations_wCNA_status/2021-04-06_Mutect2_Strelka_merged_mutations_wCNA_status.rds")

#2. mutations vs cnas from hatchet
p1=fread("p001_SNVs_explained_by_hatchet_output.txt")
p2=fread("p002_SNVs_explained_by_hatchet_output.txt")
p3=fread("p003_SNVs_explained_by_hatchet_output.txt")

all_hatchet = rbind(p1, p2, p3)
