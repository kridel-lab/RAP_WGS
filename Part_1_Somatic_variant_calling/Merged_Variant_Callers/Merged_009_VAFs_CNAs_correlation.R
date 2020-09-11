#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

date = Sys.Date()
print(date)
options(scipen=999) #no scientific notation

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)
source("/cluster/home/kisaev/scripts/bam_readcount_parseline.R")

#load libraries
packages <- c("dplyr", "ggplot2", "tidyr", "data.table", "plyr",
	"stringr", "readxl", "GenomicRanges", "params", "readr", "RColorBrewer")
lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#used bam readcount to extract counts from all samples
#overlapping mutations that were only found in some samples and not all
#Pyclone input requires these values

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#Our mutations
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone")
#sample info
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#want to get an idea of how VAFs correspond to the predicted copy number of the
#region they are found in

read_only$isdriver=""
read_only$isdriver[which(read_only$symbol %in% reddy$Gene)] = "yes"

#save plots in data folder
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/data")

#make density VAF plot colored by sample
# Add mean lines
pdf("read_only_muts_AF_summary_by_patient.pdf")
p<-ggplot(read_only, aes(x=gt_AF, color=id)) +
  geom_density()
p
dev.off()

pdf("read_only_muts_AF_summary_by_CNA.pdf")
p<-ggplot(read_only, aes(x=gt_AF, color=Corrected_Call)) +
  geom_density()
p
dev.off()

drivers = as.data.table(filter(read_only, isdriver == "yes", ExonicFunc.ensGene=="nonsynonymous_SNV"))
medians = as.data.table(drivers %>% group_by(symbol) %>% dplyr::summarize(median=median(gt_AF)))
medians = medians[order(median)]
drivers$symbol = factor(drivers$symbol, levels=medians$symbol)

pdf("read_only_muts_AF_summary_by_driver_gene.pdf", width=10)
ggplot(drivers, aes(x=symbol, y=gt_AF, color=Corrected_Call)) +
  geom_point(aes(colour = Corrected_Call),alpha=0.6, size=2) + theme_minimal()+
	theme(axis.text.x = element_text(angle = 90,size=8))+
	 scale_color_brewer(palette="Dark2")
dev.off()
