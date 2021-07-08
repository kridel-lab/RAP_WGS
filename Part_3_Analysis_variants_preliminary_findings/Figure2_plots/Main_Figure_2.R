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
snvs_all$type = ""
snvs_all$mut = "SNV/Indel"

svs = readRDS("Figure2_input_data_sample_dist_SVs.rds")
svs$type = ""
svs$mut = "SV"
z = which(!(filter(snvs_all, Patient == "DLCBL double hit stage IV")$num_of_samples_with_mut %in% filter(svs, Patient == "DLCBL double hit stage IV")$num_of_samples_with_mut))
missing = as.character(filter(snvs_all, Patient == "DLCBL double hit stage IV")$num_of_samples_with_mut[z])
add_missing = as.data.frame(matrix(ncol=ncol(svs), nrow=length(missing)))
colnames(add_missing) = colnames(svs)
add_missing$num_of_samples_with_mut = as.character(missing)
add_missing$Patient = "DLCBL double hit stage IV"
add_missing$num_of_muts = 0
add_missing$type=""
add_missing$mut="SV"
svs = rbind(svs, add_missing)

cnas = readRDS("Figure2_input_data_sample_dist_CNAs.rds")
cnas$type = cnas$CNA
cnas$CNA = NULL
cnas$mut = "CNA"
colnames(cnas)[2] = "Patient"

all_data = rbind(snvs_all, svs, cnas)

all_data$Patient = factor(all_data$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

all_data$mut = factor(all_data$mut, levels=c("SNV/Indel", "CNA", "SV"))
all_data$type = factor(all_data$type, levels=c("", "Amplification", "Deletion"))

dir.create(file.path(date))
setwd(file.path(date))

#----------------------------------------------------------------------
#make plot
#----------------------------------------------------------------------

pdf("Figure2_main_right_p1.pdf", width=6, height=4)
p<-ggline(filter(all_data, Patient == "MCL blastoid stage IV"),
x="num_of_samples_with_mut", y="num_of_muts", color="type")+theme_bw()+
xlab("# of samples sharing mutation") + ylab("Count")+
ggforce::facet_row(. ~ mut, scales="free", space='free')+
scale_color_manual(values=c("black", "red", "blue"))+
theme(axis.text.x = element_text(color="black"),
axis.text.y = element_text(color="black"))
p1=ggpar(p,legend="none")
print(p1)
dev.off()

pdf("Figure2_main_right_p2.pdf", width=6, height=4)
p<-ggline(filter(all_data, Patient == "PMBCL stage IV bulky B symptoms"),
x="num_of_samples_with_mut", y="num_of_muts", color="type")+theme_bw()+
xlab("# of samples sharing mutation") + ylab("Count")+
ggforce::facet_row(. ~ mut, scales="free", space='free')+
scale_color_manual(values=c("black", "red", "blue"))+
theme(axis.text.x = element_text(color="black"),
axis.text.y = element_text(color="black"))
p2=ggpar(p,legend="none")
print(p2)
dev.off()

pdf("Figure2_main_right_p3.pdf", width=6, height=4)
p<-ggline(filter(all_data, Patient == "DLCBL double hit stage IV"),
x="num_of_samples_with_mut", y="num_of_muts", color="type")+theme_bw()+
xlab("# of samples sharing mutation") + ylab("Count")+
ggforce::facet_row(. ~ mut, scales="free", space='free')+
scale_color_manual(values=c("black", "red", "blue"))+
theme(axis.text.x = element_text(color="black", size=3),
axis.text.y = element_text(color="black"))
p3=ggpar(p,legend="none")
print(p3)
dev.off()

#Combine everything ---------------------------------------------------
all_plots = plot_grid(
  p1,
  p2,
  p3, align = "v", ncol = 1)

pdf("Figure2_main_all_patients.pdf", width=5, height=10)
print(all_plots)
dev.off()

#get one plot with legend
pdf("legend_for_figure2.pdf", width=6, height=4)
p<-ggline(filter(all_data, Patient == "DLCBL double hit stage IV"),
x="num_of_samples_with_mut", y="num_of_muts", color="type")+theme_bw()+
xlab("# of samples sharing mutation") + ylab("Count")+
ggforce::facet_row(. ~ mut, scales="free", space='free')+
scale_color_manual(values=c("black", "red", "blue"))+
theme(axis.text.x = element_text(color="black", size=3),
axis.text.y = element_text(color="black"))
p3=ggpar(p,legend="bottom")
print(p3)
dev.off()
