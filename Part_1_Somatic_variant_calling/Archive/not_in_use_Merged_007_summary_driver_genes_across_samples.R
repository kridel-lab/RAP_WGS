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

#load libraries
packages <- c("dplyr", "ggplot2", "tidyr", "data.table", "plyr",
	"stringr")
lapply(packages, require, character.only = TRUE)
library("readxl")
library(GenomicRanges)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#CNAs annotated by gene that they overlap (no SNVs)
#setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/TITAN_CNA/results/titan/hmm/optimalClusterSolution_files/titanCNA_ploidy2")
#cnas_all = fread(list.files(pattern="all_CNAs_protein_coding_samples.txt")[length(list.files(pattern="all_CNAs_protein_coding_samples.txt"))])

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#sample info
#samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")
samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#Our mutations
#muts = fread(list.files(pattern="all_SNVs_samples.txt")[length(list.files(pattern="all_SNVs_samples.txt"))])
read_only = fread(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.txt"))])

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#mut-gene summary table for downstream use
mut_gene = unique(read_only[,c("mut_id", "symbol", "biotype", "Func.ensGene", "ExonicFunc.ensGene",
"AAChange.ensGene", "cosmic68")])

#1. how many patient samples is each mutation found in?
samples_per_mut = as.data.table(table(read_only$mut_id))
colnames(samples_per_mut) = c("mut_id", "num_of_samples_with_mut")
samples_per_mut = merge(samples_per_mut, mut_gene, by="mut_id")
samples_per_mut = samples_per_mut[order(-num_of_samples_with_mut)]
samples_per_mut$phylogeny = ""
samples_per_mut$phylogeny[samples_per_mut$num_of_samples_with_mut == 20] = "ancestor"
samples_per_mut$phylogeny[samples_per_mut$num_of_samples_with_mut == 1] = "private"
samples_per_mut$phylogeny[samples_per_mut$phylogeny == ""] = "shared"

#data table for barplot
barplot = as.data.table(table(samples_per_mut$num_of_samples_with_mut))
colnames(barplot) = c("num_of_samples_with_mut", "num_of_muts")
barplot$num_of_samples_with_mut = factor(barplot$num_of_samples_with_mut, levels=barplot$num_of_samples_with_mut)

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/001_samples_per_mutation.pdf")
# Basic barplot
p<-ggplot(data=barplot, aes(x=num_of_samples_with_mut, y=num_of_muts)) +
  geom_bar(stat="identity")+theme_minimal()+ggtitle("Number of samples with a given mutation")
print(p)
dev.off()

#summarize number of mutations per sample
muts_per_sample = as.data.table(table(read_only$id))
muts_per_sample = muts_per_sample[order(-N)]
muts_per_sample$V1 = factor(muts_per_sample$V1, levels=muts_per_sample$V1)
colnames(muts_per_sample) = c("Sample", "num_of_muts")

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/002_mutations_per_sample.pdf")
# Basic barplot
p<-ggplot(data=muts_per_sample, aes(x=Sample, y=num_of_muts)) +
  geom_bar(stat="identity")+theme_minimal()+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Sample")+
	ylab("Number of mutations")
print(p)
dev.off()

#summarize biotypes - types of genes that are mutated
biotypes = as.data.table(table(read_only$biotype))
colnames(biotypes) = c("gene_type", "num_mutations_found_in_gene_type")
biotypes = biotypes[order(-num_mutations_found_in_gene_type)]
biotypes$gene_type = factor(biotypes$gene_type, levels=biotypes$gene_type)

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/003_all_mutation_overview_gene_types.pdf")
# Basic barplot
p<-ggplot(data=biotypes, aes(x=gene_type, y=num_mutations_found_in_gene_type)) +
  geom_bar(stat="identity")+theme_minimal()
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("Number of mutations found in each type of gene")
dev.off()

#biotypes but only unique mutations
biotypes = as.data.table(table(samples_per_mut$biotype))
colnames(biotypes) = c("gene_type", "num_mutations_found_in_gene_type")
biotypes = biotypes[order(-num_mutations_found_in_gene_type)]
biotypes$gene_type = factor(biotypes$gene_type, levels=biotypes$gene_type)

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/004_unique_mutation_overview_gene_types.pdf")
# Basic barplot
p<-ggplot(data=biotypes, aes(x=gene_type, y=num_mutations_found_in_gene_type)) +
  geom_bar(stat="identity")+theme_minimal()
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle("Number of mutations found in each type of gene")
dev.off()

#2. how many mutations found in all samples fall into 'driver genes' vs
#non driver genes
samples_per_mut$driver = ""
z = which(samples_per_mut$symbol %in% reddy$Gene)
samples_per_mut$driver[z] = "driver"
samples_per_mut$driver[-z] = "non_driver"

paste(length(unique(samples_per_mut$mut_id)), "unique mutations")
paste(length(unique(samples_per_mut$mut_id[samples_per_mut$driver == "driver"])), "unique mutations in driver genes including in non-exon regions")

samples_per_mut$morin = ""
z = which(samples_per_mut$symbol %in% genes_sum$Gene)
samples_per_mut$morin[z] = "morin"
samples_per_mut$morin[-z] = "not_in_morin"

write.table(samples_per_mut, file=paste(date, "mutation_summary_occurence_phylogeny_driver_gene_status.txt", sep="_"),
quote=F, row.names=F, sep=";")

#keep only functional protein coding gene mutations
samples_per_mut = as.data.table(filter(samples_per_mut, biotype == "protein_coding",
Func.ensGene %in% c("exonic", "splicing"), !(ExonicFunc.ensGene == "synonymous_SNV")))

drivers = as.data.table(table(samples_per_mut$driver, samples_per_mut$phylogeny))
colnames(drivers) = c("Driver", "Phylogeny", "Num_muts")
morin = as.data.table(table(samples_per_mut$morin, samples_per_mut$phylogeny))
colnames(morin) = c("Morin", "Phylogeny", "Num_muts")

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/005_drivers_vs_phylogeny_exonic_splicing_nonsynon_only.pdf")
# Basic barplot
p<-ggplot(data=drivers, aes(x=Driver, y=Num_muts, fill=Phylogeny)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal()
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
scale_fill_brewer(palette="Accent")+ggtitle("protein-coding genes mutations drivers vs non-drivers")
dev.off()

#which driver genes are found in both more than one phylogeny
convergent=as.data.table(table(samples_per_mut$symbol, samples_per_mut$phylogeny))
convergent = as.data.table(filter(convergent, N >0))
convergent=as.data.table(table(convergent$V1))
convergent=convergent[order(-N)]
convergent_cands = filter(convergent, N >=2)$V1

#check status of convergent genes
convergent_muts = as.data.table(filter(samples_per_mut, symbol %in% convergent_cands))
convergent_muts$phylogeny = factor(convergent_muts$phylogeny, levels=c("ancestor",
"shared", "private"))
convergent_muts$symbol = factor(convergent_muts$symbol, levels=convergent_cands)

#summarize "convergent genes" (found in more than 2 phylogeny types)
pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/006_convergent_genes_cands_summary_exonic_splicing_nonsynon_only.pdf")
p = ggplot(convergent_muts, aes(phylogeny, symbol)) +
  geom_tile(aes(fill = driver), colour = "grey50") +
  xlab("Phylogeny") + ylab("Gene") +
	ggtitle("Summary of protein-coding genes mutated in more than 1 phylogeny")
print(p)
dev.off()

#look at shared genes only
#are there examples of genes mutated in multiple ways across shared samples?
drivers = as.data.table(filter(samples_per_mut, driver=="driver",
!(ExonicFunc.ensGene == "synonymous_SNV")))

drivers$gene_mut = paste(drivers$symbol, drivers$mut_id)
pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/007_driver_genes_muts_summary_exonic_splicing_nonsynon_only.pdf")
p = ggplot(drivers, aes(phylogeny, gene_mut)) +
  geom_tile(aes(fill = morin), colour = "grey50") +
  xlab("Phylogeny") + ylab("Gene") +
	ggtitle("Summary of mutations in driver genes")
print(p)
dev.off()

#summarize copy number status of driver genes
cnas_driver = as.data.table(filter(cnas_all, symbol %in% drivers$symbol))
colnames(cnas_driver)[51] = "Sample_Name"
cnas_driver$TITAN_call = factor(cnas_driver$TITAN_call, levels=c("DLOH", "NLOH",
"HET", "ALOH", "ASCNA", "GAIN", "BCNA", "UBCNA"))
pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/008_driver_genes_CNA_status.pdf")
p = ggplot(cnas_driver, aes(Sample, symbol)) +
  geom_tile(aes(fill = TITAN_call), colour = "grey50") +
  xlab("Sample") + ylab("Gene") +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
	scale_fill_manual(values=c("dodgerblue1", "grey", "darkolivegreen3","salmon",
	"orangered1", "orangered3", "red4", "violetred4"))+
	ggtitle("Summary of driver genes CNAs")
print(p)
dev.off()
