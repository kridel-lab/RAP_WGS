#-------------------------------------------------------------------------------
#Make main figure 1 combining summaries of all dimensions (SNVs, CNAs, SVs...)
#Karin Isaev
#-------------------------------------------------------------------------------

#load packages and data
library(dplyr)
library(ggpubr)
library(data.table)
library(BioCircos)
library(plyr)
require(gridExtra)
library(cowplot)
library(ggforce)
library(wesanderson)

options(scipen=999)
date=Sys.Date()

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files/Figure2_MAIN")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

snvs_all = readRDS("Figure2_input_data_sample_dist_SNVindels.rds")
svs = fread("")
cnas = fread("")

snvs_all$Patient = factor(snvs_all$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
svs$Patient = factor(svs$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
cnas$Patient = factor(cnas$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

dir.create(file.path(date))
setwd(file.path(date))
