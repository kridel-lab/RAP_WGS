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

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

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

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/001_samples_per_mutation.pdf")
# Basic barplot
p<-ggplot(data=barplot, aes(x=num_of_samples_with_mut, y=num_of_muts, fill=patient)) +
  geom_bar(stat="identity")+theme_minimal() + #+ggtitle("Number of samples with a given mutation") +
	xlab("Number of samples with mutation") + ylab("Number of mutations")+
	facet_grid(. ~ patient, scales="free", space='free')+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(), legend.position="bottom",
	legend.text = element_text(size=6))+
	scale_y_continuous(breaks=seq(0, 300000, by = 25000))

print(p)
dev.off()

colnames(barplot)[2]="Patient"

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/001_samples_per_mutation_lineplot.pdf",
width=5, height=5)
# Basic barplot
p<-ggline(barplot, x="num_of_samples_with_mut", y="num_of_muts",
palette = c("#00AFBB", "#E7B800", "#FC4E07"), color="Patient")+
xlab("# of samples sharing mutation") + ylab("Mutation count")+
scale_y_continuous(breaks=seq(0,300000,20000))
print(p)
dev.off()

#summarize number of mutations per sample
muts_per_sample = as.data.table(table(read_only$STUDY_PATIENT_ID,read_only$Tissue_Site))
muts_per_sample = as.data.table(filter(muts_per_sample, N >0))
muts_per_sample = muts_per_sample[order(-N)]
colnames(muts_per_sample) = c("Patient", "Sample", "num_of_muts")
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
muts_per_sample$Patient = factor(muts_per_sample$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
muts_per_sample$Sample = factor(muts_per_sample$Sample, levels=unique(muts_per_sample$Sample))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/002_mutations_per_sample.pdf",
width=4, height=5)
# Basic barplot
p<-ggplot(data=muts_per_sample, aes(x=Sample, y=num_of_muts, fill=Patient)) +
  geom_bar(stat="identity")+theme_classic()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1),
  legend.position = "none") + xlab("Sample")+
	ylab("Mutation count")+
	facet_grid(. ~ Patient, scales="free", space='free')+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
	legend.text = element_text(size=6))+
	scale_y_continuous(breaks=seq(0, 400000, by = 50000))+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))
print(p)
dev.off()

#summarize biotypes - types of genes that are mutated
biotypes = as.data.table(table(read_only$Func.ensGene, read_only$symbol, read_only$STUDY_PATIENT_ID))
colnames(biotypes) = c("gene_type", "gene", "patient", "num_muts_in_gene_type")
biotypes = as.data.table(filter(biotypes, num_muts_in_gene_type > 0))
biotypes = biotypes[order(-num_muts_in_gene_type)]
biotypes$gene_type = factor(biotypes$gene_type, levels=unique(biotypes$gene_type))
biotypes$driver = ""
z = which((biotypes$gene %in% all_drivers$Gene[all_drivers$type=="MCL"]) &(biotypes$patient == "LY_RAP_0001"))
biotypes$driver[z] = "driver"
z = which((biotypes$gene %in% all_drivers$Gene[all_drivers$type=="PMBCL"]) &(biotypes$patient == "LY_RAP_0002"))
biotypes$driver[z] = "driver"
z = which((biotypes$gene %in% all_drivers$Gene[all_drivers$type=="DLBCL"]) &(biotypes$patient == "LY_RAP_0003"))
biotypes$driver[z] = "driver"
biotypes$driver[biotypes$driver == ""] = "non_driver"
biotypes$patient[biotypes$patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
biotypes$patient[biotypes$patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
biotypes$patient[biotypes$patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
biotypes$patient = factor(biotypes$patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

introns = filter(biotypes, gene_type=="intronic")
introns$N = log1p(introns$num_muts_in_gene_type)
introns$driver = factor(introns$driver, levels=c("driver", "non-driver"))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/intronic_mutations_drivers_vs_nondrivers.pdf",
height=4, width=4)
p<-ggboxplot(data=introns, x="driver", y="N", fill="patient", facet.by ="patient") +
  theme_classic() +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))+
	theme(strip.background = element_blank(),
  strip.text.x = element_blank(), legend.position="none")#+
	#scale_y_continuous(breaks=seq(0, 400000, by = 25000))
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 stat_compare_means(label = "p.format") + ylab("log1p(# of mutations per gene)")
dev.off()

biotypes$N = log1p(biotypes$num_muts_in_gene_type)
biotypes$driver = factor(biotypes$driver, levels=c("driver", "non-driver"))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_muts_mutations_drivers_vs_nondrivers.pdf",
height=4, width=4)
p<-ggboxplot(data=biotypes, x="driver", y="N", fill="patient", facet.by ="patient") +
  theme_classic() +
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))+
	theme(strip.background = element_blank(),
  strip.text.x = element_blank(), legend.position="none")#+
	#scale_y_continuous(breaks=seq(0, 400000, by = 25000))
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 stat_compare_means(label = "p.format") + ylab("log1p(# of mutations per gene)")
dev.off()

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/003_all_mutation_overview_gene_types.pdf")
# Basic barplot
p<-ggplot(data=biotypes, aes(x=gene_type, y=num_muts_in_gene_type, fill=patient)) +
  geom_bar(stat="identity")+theme_minimal() +
	facet_grid(. ~ patient, scales="free", space='free')+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(), legend.position="bottom",
	legend.text = element_text(size=6))#+
	#scale_y_continuous(breaks=seq(0, 400000, by = 25000))
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ggtitle("Number of mutations found in each type of gene")
dev.off()

#2. how many mutations found in all samples fall into 'driver genes' vs
#non driver genes
samples_per_mut$driver = ""
z = which((samples_per_mut$symbol %in% all_drivers$Gene[all_drivers$type=="MCL"]) &(samples_per_mut$patient == "LY_RAP_0001"))
samples_per_mut$driver[z] = "driver"
z = which((samples_per_mut$symbol %in% all_drivers$Gene[all_drivers$type=="PMBCL"]) &(samples_per_mut$patient == "LY_RAP_0002"))
samples_per_mut$driver[z] = "driver"
z = which((samples_per_mut$symbol %in% all_drivers$Gene[all_drivers$type=="DLBCL"]) &(samples_per_mut$patient == "LY_RAP_0003"))
samples_per_mut$driver[z] = "driver"
samples_per_mut$driver[samples_per_mut$driver == ""] = "non_driver"

paste(length(unique(samples_per_mut$mut_id)), "unique mutations")
paste(length(unique(samples_per_mut$mut_id[samples_per_mut$driver == "driver"])), "unique mutations in driver genes including in non-exon regions")

samples_per_mut$morin = ""
z = which(samples_per_mut$symbol %in% genes_sum$Gene)
samples_per_mut$morin[z] = "morin"
samples_per_mut$morin[-z] = "not_in_morin"

summ_d_a = as.data.table(table(samples_per_mut$patient, samples_per_mut$phylogeny, samples_per_mut$driver))
colnames(summ_d_a) = c("Patient", "Phylogeny", "Driver", "N")
summ_d_a$Phylogeny[summ_d_a$Phylogeny %in% c("private", "shared")] = "Non-Truncal"
summ_d_a$Phylogeny[summ_d_a$Phylogeny %in% c("ancestor")] = "Truncal"
summ_d_a$N_log1p = log1p(summ_d_a$N)
summ_d_a$Phylogeny = factor(summ_d_a$Phylogeny, levels = c("Truncal", "Non-Truncal"))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_muts_mutations_drivers_vs_nondrivers_vs_ancestor_status.pdf",
height=6, width=3)
p<-ggboxplot(data=summ_d_a, x="Driver", y="N_log1p", fill="Phylogeny",
palette = get_palette("Dark2", 3)) +
  theme_classic() +
  #scale_fill_manual(values=c(""))+
	theme(strip.background = element_blank(),
  strip.text.x = element_blank(), legend.position="bottom")#+
	#scale_y_continuous(breaks=seq(0, 400000, by = 25000))
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 stat_compare_means(aes(group = Phylogeny)) +
 ylab("log1p(# of mutations)")
dev.off()

#keep only functional protein coding gene mutations
summ_d_a = as.data.table(filter(samples_per_mut, biotype == "protein_coding",
Func.ensGene %in% c("exonic", "splicing"), !(ExonicFunc.ensGene == "synonymous_SNV")))
summ_d_a = as.data.table(table(summ_d_a$patient, summ_d_a$phylogeny, summ_d_a$driver))
colnames(summ_d_a) = c("Patient", "Phylogeny", "Driver", "N")
summ_d_a$Phylogeny[summ_d_a$Phylogeny %in% c("private", "shared")] = "Non-Truncal"
summ_d_a$Phylogeny[summ_d_a$Phylogeny %in% c("ancestor")] = "Truncal"
summ_d_a$N_log1p = log1p(summ_d_a$N)
summ_d_a$Phylogeny = factor(summ_d_a$Phylogeny, levels = c("Truncal", "Non-Truncal"))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_muts_mutations_drivers_vs_nondrivers_vs_ancestor_status_exonic_only.pdf",
height=6, width=3)
p<-ggboxplot(data=summ_d_a, x="Driver", y="N_log1p", fill="Phylogeny",
palette = get_palette("Dark2", 3)) +
  theme_classic() +
  #scale_fill_manual(values=c(""))+
	theme(strip.background = element_blank(),
  strip.text.x = element_blank(), legend.position="bottom")#+
	#scale_y_continuous(breaks=seq(0, 400000, by = 25000))
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
 stat_compare_means(aes(group = Phylogeny)) +
 ylab("log1p(# of mutations)")
dev.off()

write.table(samples_per_mut, file=paste(date, "mutation_summary_occurence_phylogeny_driver_gene_status.txt", sep="_"),
quote=F, row.names=F, sep=";")

#summarize driver genes and what kind of mutations they are affected by
drivers_all = filter(samples_per_mut, driver=="driver")
drivers_all = as.data.table(table(drivers_all$patient, drivers_all$symbol, drivers_all$phylogeny, drivers_all$Func.ensGene)) %>% filter(N > 0, V4 %in% c("exonic", "intronic", "splicing", "UTR3", "UTR5"))
colnames(drivers_all) = c("Patient", "Gene", "Status", "Mutation", "Mutation_count")
drivers_all$Status = factor(drivers_all$Status, levels=c("ancestor", "shared", "private"))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/driver_genes_vs_shared_status_type.pdf", height=6,width=4)
p = ggplot(filter(drivers_all, Status == "ancestor"), aes(Status, Gene)) +
theme_bw()+
  geom_tile(aes(fill = Status), colour = "black") +
scale_fill_rickandmorty()+
  geom_text(aes(label=Mutation_count),size=1) +
	facet_grid(Patient ~ Mutation, scales="free", space='free')+
theme(axis.text =element_text(size=4),
axis.title.x=element_blank(),
legend.position = "none",
strip.text = element_text(size = 4))+xlab("")
print(p)
dev.off()

#keep only functional protein coding gene mutations
samples_per_mut = as.data.table(filter(samples_per_mut, biotype == "protein_coding",
Func.ensGene %in% c("exonic", "splicing"), !(ExonicFunc.ensGene == "synonymous_SNV")))

drivers = as.data.table(table(samples_per_mut$patient, samples_per_mut$driver, samples_per_mut$phylogeny))
colnames(drivers) = c("patient", "Driver", "Phylogeny", "Num_muts")
morin = as.data.table(table(samples_per_mut$patient, samples_per_mut$morin, samples_per_mut$phylogeny))
colnames(morin) = c("patient", "Morin", "Phylogeny", "Num_muts")
drivers$patient[drivers$patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
drivers$patient[drivers$patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
drivers$patient[drivers$patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
drivers$patient = factor(drivers$patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
drivers$Phylogeny = factor(drivers$Phylogeny, levels=c("ancestor", "shared", "private"))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/005_drivers_vs_phylogeny_exonic_splicing_nonsynon_only.pdf", width=7, height=6)
# Basic barplot
p<-ggplot(data=drivers, aes(x=Driver, y=Num_muts, fill=Phylogeny)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+theme_minimal()
p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylab("# of mutations")+
 facet_grid(. ~ patient) +
 geom_text(aes(label = Num_muts), size = 2, position=position_dodge(width=0.9), vjust=-0.25) +
scale_fill_brewer(palette="Accent")#+ggtitle("protein-coding genes mutations drivers vs non-drivers")
dev.off()

#which driver genes are found in both more than one phylogeny
convergent = as.data.table(table(samples_per_mut$patient, samples_per_mut$symbol, samples_per_mut$phylogeny))
convergent = as.data.table(filter(convergent, N >0))
convergent = as.data.table(table(convergent$V1, convergent$V2))
convergent = as.data.table(filter(convergent, N >0))
convergent = convergent[order(-N)]
convergent_cands = filter(convergent, N >=2)$V2

#check status of convergent genes
convergent_muts = as.data.table(filter(samples_per_mut, symbol %in% convergent_cands))
convergent_muts$phylogeny = factor(convergent_muts$phylogeny, levels=c("ancestor",
"shared", "private"))
convergent_muts$symbol = factor(convergent_muts$symbol, levels=unique(convergent_cands))

#summarize "convergent genes" (found in more than 2 phylogeny types)
pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/006_convergent_genes_cands_summary_exonic_splicing_nonsynon_only.pdf", width=9, height=15)
p = ggplot(convergent_muts, aes(phylogeny, symbol)) +
  geom_tile(aes(fill = driver), colour = "grey50") +
  xlab("Phylogeny") + ylab("Gene") +
	facet_grid(. ~ patient)+
	theme(axis.text=element_text(size=6))+
	ggtitle("Summary of protein-coding genes mutated in more than 1 phylogeny")
print(p)
dev.off()

#look at shared genes only
#are there examples of genes mutated in multiple ways across shared samples?
drivers = as.data.table(filter(samples_per_mut, driver=="driver"))
drivers$patient[drivers$patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
drivers$patient[drivers$patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
drivers$patient[drivers$patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
drivers$patient = factor(drivers$patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))
drivers$phylogeny = factor(drivers$phylogeny, levels=c("ancestor", "shared", "private"))

drivers$gene_mut = paste(drivers$symbol, drivers$mut_id)
pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/007_driver_genes_muts_summary_exonic_splicing_nonsynon_only.pdf", width=10, height=6)
p = ggplot(drivers, aes(phylogeny, gene_mut)) +
  geom_tile(aes(fill = mut_type), colour = "grey50") +
  xlab("Phylogeny") + ylab("Gene") +
	theme_minimal()+
	facet_grid(. ~ patient, scales="free", space='free')+coord_flip()+
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
#	ggtitle("Summary of mutations in driver genes")
print(p)
dev.off()

drivers$mut_type = factor(drivers$mut_type, levels=c("SNV", "INDEL"))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/007_driver_genes_muts_summary_exonic_splicing_nonsynon_only_001.pdf", width=6, height=6)
p = ggplot(filter(drivers, patient ==  "MCL blastoid stage IV"), aes(phylogeny, gene_mut)) +
  geom_tile(aes(fill = mut_type), colour = "grey50") +
  xlab("Phylogeny") + ylab("Gene") +
	theme_minimal()+
	#coord_flip()+
	theme(axis.text.x = element_text(angle = 0, vjust = 1, size=10))
#	ggtitle("Summary of mutations in driver genes")
print(p)
dev.off()

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/007_driver_genes_muts_summary_exonic_splicing_nonsynon_only_002.pdf", width=6, height=9)
p = ggplot(filter(drivers, patient ==  "PMBCL stage IV bulky B symptoms"), aes(phylogeny, gene_mut)) +
  geom_tile(aes(fill = mut_type), colour = "grey50") +
  xlab("Phylogeny") + ylab("Gene") +
	#theme(axis.text.x = element_text(size=3), axis.text.y=element_text(size=1))+#, angle = 0, hjust = 1, vjust = 1))+
	theme_minimal() + theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 5))

#	ggtitle("Summary of mutations in driver genes")
print(p)

dev.off()

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/007_driver_genes_muts_summary_exonic_splicing_nonsynon_only_003.pdf", width=6, height=9)
p = ggplot(filter(drivers, patient ==  "DLCBL double hit stage IV"), aes(phylogeny, gene_mut)) +
  geom_tile(aes(fill = mut_type), colour = "grey50") +
  xlab("Phylogeny") + ylab("Gene") +
	theme_minimal()+
	#coord_flip()+
	theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1))
#	ggtitle("Summary of mutations in driver genes")
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
