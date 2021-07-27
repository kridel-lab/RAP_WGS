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

#setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files/Figure1_MAIN")
setwd("/Users/kerenisaev/Documents/UHN/kridel-lab - RAP_WGS/Data-Files/Figure1_MAIN")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

snvs_all = fread("Figure1_MAIN_SNVs_Indels_ALL.txt")
snvs_coding = fread("Figure1_MAIN_SNVs_Indels_CODING.txt")
svs = fread("Figure1_MAIN_SVs.txt")
cnas = fread("Figure1_MAIN_CNAs_ALL.txt")
samples_order = fread("Figure1_MAIN_sample_order.txt", sep="")
ploidy = fread("Figure1_MAIN_purity_ploidy.txt")

#add tissue side to cnas
cnas = merge(cnas, ploidy, by=c("Sample", "Patient"))
cnas$Sample=NULL
cnas$Sample = cnas$Tissue_Site
cnas$Ploidy = NULL
cnas$Purity = NULL
cnas$perc_genome = cnas$tot_bps_cna / cnas$sum_widths

ploidy$Sample = NULL
ploidy$Sample = ploidy$Tissue_Site

snvs_all$Sample = factor(snvs_all$Sample, levels=samples_order$x)
snvs_all$Patient = factor(snvs_all$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

snvs_coding$Sample = factor(snvs_coding$Sample, levels=samples_order$x)
snvs_coding$Patient = factor(snvs_coding$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

svs$Sample = factor(svs$Sample, levels=samples_order$x)
svs$Patient = factor(svs$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

cnas$Sample = factor(cnas$Sample, levels=samples_order$x)
cnas$Patient = factor(cnas$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

ploidy$Sample = factor(ploidy$Sample, levels=samples_order$x)
ploidy$Patient = factor(ploidy$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#Main barplot summarizing all somatic events across each sample

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#All genomewide SNVs + Indels -----------------------------------------

full_snvs_indels = ggplot(data=snvs_all, aes(x=Sample, y=num_of_muts, fill=Patient)) +
  geom_bar(stat="identity")+theme_bw()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=5),
  axis.text.y = element_text(size=4),
  legend.position = "right") +
	ylab("Mutation count")+xlab("")+
  ggforce::facet_row(. ~ Patient, scales="free", space='free')+
	theme(legend.title=element_text(size=2),
    axis.title=element_text(size=6),
    legend.key.size = unit(0.3, "cm"),
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
	legend.text = element_text(size=6))+
	#scale_y_continuous(breaks=seq(0, 400000, by = 50000))+
  scale_fill_manual(values=c("snow4", "steelblue4", "rosybrown4"))


#Coding only SNVs + Indels ------------------------------------------

coding_snvs_indels = ggplot(data=snvs_coding, aes(x=Sample, y=num_of_muts, fill=mut_type)) +
  geom_bar(stat="identity")+theme_bw()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=5),
  axis.text.y = element_text(size=4),
  legend.position = "right") +
	ylab("Coding Mutation count")+xlab("")+
	ggforce::facet_row(. ~ Patient, scales="free", space='free')+
	theme(legend.title=element_text(size=2),
    axis.title=element_text(size=6),
  legend.key.size = unit(0.3, "cm"),
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
	legend.text = element_text(size=6))+
	#scale_y_continuous(breaks=seq(0, 3000, by = 300))+
  scale_fill_manual(values = wes_palette("Rushmore1", n = 5))

#CNAs ----------------------------------------------------------------

cnas$CNA = factor(cnas$CNA, levels=c(">5_N_Gain", "5_N_Gain",
"4_N_Gain", "3_N_Gain", "Neutral", "Somatic LOH", "Hemizygous Del", "Homozygous Del"))

cnas_plot = ggplot(data=cnas, aes(x=Sample, y=perc_genome, fill=CNA)) +
    geom_bar(stat="identity")+theme_bw()+#+ggtitle("Number of mutations per sample")+
  	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=4),
    axis.text.y = element_text(size=4),
    legend.position = "right") +
  	ylab("CNA % of genome")+xlab("")+
    ggforce::facet_row(. ~ Patient, scales="free", space='free')+
  	theme(legend.title=element_text(size=2),
      axis.title=element_text(size=6),
    legend.key.size = unit(0.3, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
  	legend.text = element_text(size=6))+
#  	scale_y_continuous(breaks=seq(0, 4000, by = 200))+
scale_fill_brewer(palette="RdYlBu")

#SVs -----------------------------------------------------------------

svs_plot = ggplot(data=svs, aes(x=Sample, y=num_of_muts, fill=type_SV)) +
  geom_bar(stat="identity")+theme_bw()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=4),
  axis.text.y = element_text(size=4),
  legend.position = "right") + xlab("Sample")+
	ylab("Structural variants count")+xlab("")+
	ggforce::facet_row(. ~ Patient, scales="free", space='free')+
	theme(legend.title=element_text(size=2),
    axis.title=element_text(size=4),
    legend.key.size = unit(0.3, "cm"),
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
	legend.text = element_text(size=6))+
	#scale_y_continuous(breaks=seq(0, 200, by = 20))+
  scale_fill_manual(values = wes_palette("Moonrise2", n = 4))


#Ploidy and Purity ----------------------------------------------------

ploidy = melt(ploidy)

purity_plot = ggplot(filter(ploidy, variable == "Purity"), aes(Sample, variable)) +
 geom_raster(aes(fill = value))+theme_bw()+
 theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=4),
 axis.text.y = element_text(size=4),
  legend.position = "right") + xlab("Sample")+
 ylab("Purity")+xlab("")+
 ggforce::facet_row(. ~ Patient, scales="free", space='free')+
 theme(
    axis.title=element_text(size=4),
    legend.key.size = unit(0.2, "cm"),
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_text(size=2),
 legend.text = element_text(size=2))+
  scale_fill_gradient(low = "yellow", high = "red")

ploidy_plot = ggplot(filter(ploidy, variable == "Ploidy"), aes(Sample, variable)) +
geom_raster(aes(fill = value))+theme_bw()+
theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=4),
axis.text.y = element_text(size=4),
legend.position = "right") + xlab("Sample")+
 ylab("Ploidy")+
 ggforce::facet_row(. ~ Patient, scales="free", space='free')+
 theme(
      axis.title=element_text(size=4),
      legend.key.size = unit(0.2, "cm"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.title=element_text(size=2),
   legend.text = element_text(size=2))+
    scale_fill_gradient(low = "blue", high = "red")

#Combine everything ---------------------------------------------------
all_plots = plot_grid(
  full_snvs_indels + theme(legend.justification = c(0,1)),
  coding_snvs_indels + theme(legend.justification = c(0,1)),
  cnas_plot + theme(legend.justification = c(0,1)),
  svs_plot + theme(legend.justification = c(0,1)),
  purity_plot + theme(legend.justification = c(0,1)),
  ploidy_plot + theme(legend.justification = c(0,1)),
  align = "v", ncol = 1, rel_heights=c(3, 2, 2, 2, 0.7, 2))

dir.create(file.path(date))
setwd(file.path(date))

pdf("001_MAIN_Figure.pdf")
print(all_plots)
dev.off()
