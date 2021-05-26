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

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files/Figure1_MAIN")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

snvs_all = fread()
snvs_coding = fread()
samples_order = fread()

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#All genomewide SNVs + Indels -----------------------------------------

full_snvs_indels = ggplot(data=snvs_all, aes(x=Sample, y=num_of_muts, fill=Patient)) +
  geom_bar(stat="identity")+theme_classic()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=5),
  legend.position = "right") + xlab("Sample")+
	ylab("Mutation count")+
  facet_grid(. ~ Patient, scales="free", space='free')+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
	legend.text = element_text(size=6))+
	#scale_y_continuous(breaks=seq(0, 400000, by = 50000))+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))


#Coding only SNVs + Indels ------------------------------------------

coding_snvs_indels = ggplot(data=snvs_coding, aes(x=Sample, y=num_of_muts, fill=mut_type)) +
  geom_bar(stat="identity")+theme_classic()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=5),
  legend.position = "right") + xlab("Sample")+
	ylab("Mutation count")+
	facet_grid(. ~ Patient, scales="free", space='free')+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
	legend.text = element_text(size=6))+
	#scale_y_continuous(breaks=seq(0, 3000, by = 300))+
  scale_fill_brewer(palette="Dark2")


pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/SNVs_panel.pdf", width=8, height=6)
plot_grid(full_snvs_indels, coding_snvs_indels, labels = "AUTO",
align = "v",
ncol = 1)
dev.off()
