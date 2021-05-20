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

dir.create(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))
setwd(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))

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

samples_per_mut = as.data.table(table(svs$id, svs$STUDY_PATIENT_ID))

mut_gene = unique(svs[,c("gene", "id", "SVTYPE")])

colnames(samples_per_mut) = c("id", "patient", "num_of_samples_with_mut")
samples_per_mut = merge(samples_per_mut, mut_gene, by="id")
samples_per_mut = samples_per_mut[order(-num_of_samples_with_mut)]
samples_per_mut = as.data.table(filter(samples_per_mut, num_of_samples_with_mut >0))

samples_per_mut$phylogeny = ""
z = which(samples_per_mut$patient == "LY_RAP_0001")
LY_RAP_0001 = samples_per_mut[z,]
LY_RAP_0001$phylogeny[LY_RAP_0001$num_of_samples_with_mut == 3] = "ancestor"
LY_RAP_0001$phylogeny[LY_RAP_0001$num_of_samples_with_mut == 1] = "private"
LY_RAP_0001$phylogeny[LY_RAP_0001$phylogeny == ""] = "shared"

z = which(samples_per_mut$patient == "LY_RAP_0002")
LY_RAP_0002 = samples_per_mut[z,]
LY_RAP_0002$phylogeny[LY_RAP_0002$num_of_samples_with_mut == 4] = "ancestor"
LY_RAP_0002$phylogeny[LY_RAP_0002$num_of_samples_with_mut == 1] = "private"
LY_RAP_0002$phylogeny[LY_RAP_0002$phylogeny == ""] = "shared"

z = which(samples_per_mut$patient == "LY_RAP_0003")
LY_RAP_0003 = samples_per_mut[z,]
LY_RAP_0003$phylogeny[LY_RAP_0003$num_of_samples_with_mut == 20] = "ancestor"
LY_RAP_0003$phylogeny[LY_RAP_0003$num_of_samples_with_mut == 1] = "private"
LY_RAP_0003$phylogeny[LY_RAP_0003$phylogeny == ""] = "shared"

samples_per_mut = rbind(LY_RAP_0001, LY_RAP_0002, LY_RAP_0003)

#data table for barplot
barplot = as.data.table(table(samples_per_mut$num_of_samples_with_mut, samples_per_mut$patient))
barplot = as.data.table(filter(barplot, N >0))
colnames(barplot) = c("num_of_samples_with_mut", "patient", "num_of_muts")
barplot$num_of_samples_with_mut = factor(barplot$num_of_samples_with_mut, levels=unique(barplot$num_of_samples_with_mut))
barplot$patient[barplot$patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
barplot$patient[barplot$patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
barplot$patient[barplot$patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
barplot$patient = factor(barplot$patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

colnames(barplot)[2]="Patient"
pdf("001_samples_per_mutation_lineplot_SVs.pdf", width=4, height=5)
# Basic barplot
p<-ggline(barplot, x="num_of_samples_with_mut", y="num_of_muts",
palette = c("#00AFBB", "#E7B800", "#FC4E07"), color="Patient")+
xlab("# of samples sharing mutation") + ylab("Structural variants count")+
scale_y_continuous(breaks=seq(0,1000,100))
print(p)
dev.off()

#summarize number of mutations per sample
z = which((svs$Tissue_Site == "Adrenal gland, NOS") & (svs$STUDY_PATIENT_ID == "LY_RAP_0003"))
svs$Tissue_Site[z] = "Adrenal gland"
z = which((svs$Tissue_Site == "Aorta, ascending, not specified \n\n") & (svs$STUDY_PATIENT_ID == "LY_RAP_0001"))
svs$Tissue_Site[z] = "Aorta, ascending"

muts_per_sample = as.data.table(table(svs$STUDY_PATIENT_ID, svs$Tissue_Site, svs$SVTYPE))
muts_per_sample = as.data.table(filter(muts_per_sample, N >0))
muts_per_sample = muts_per_sample[order(-N)]

colnames(muts_per_sample) = c("Patient", "Sample", "type_SV", "num_of_muts")
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
muts_per_sample$Patient = factor(muts_per_sample$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

#get order of samples for barplot
t = as.data.table(table(svs$STUDY_PATIENT_ID,svs$Tissue_Site))
t = as.data.table(filter(t, N >0))
t = t[order(-N)]

muts_per_sample$Sample = factor(muts_per_sample$Sample, levels=unique(t$V2))

pdf("002_mutations_per_sample_SVs.pdf",width=4, height=5)
# Basic barplot
p<-ggplot(data=muts_per_sample, aes(x=Sample, y=num_of_muts, fill=type_SV)) +
  geom_bar(stat="identity")+theme_classic()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=4),
  legend.position = "bottom") + xlab("Sample")+
	ylab("Structural variants count")+
	facet_grid(. ~ Patient, scales="free", space='free')+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
	legend.text = element_text(size=6))+
	scale_y_continuous(breaks=seq(0, 200, by = 20))+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07", "grey"))
print(p)
dev.off()
