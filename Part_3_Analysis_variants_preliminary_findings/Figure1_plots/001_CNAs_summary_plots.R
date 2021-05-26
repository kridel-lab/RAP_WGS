#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(ggpubr)
library("ggsci")
library(annotables)
hg19_genes = as.data.table(grch37)
#hg19_genes = filter(hg19_genes, biotype == "protein_coding")
hg19_genes = unique(hg19_genes[,c("chr", "start", "end", "symbol")]) #22810 genes
z = which(hg19_genes$chr %in% c(1:22, "X"))
hg19_genes = hg19_genes[z,]

#----------------------------------------------------------------------
#overlap CNA calls with protein coding genes coordinates
#----------------------------------------------------------------------

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")
cnas$Patient = sapply(cnas$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
cnas = filter(cnas, !(ntot==2))
cnas$CNA[cnas$ntot >2] = "Amplification"
cnas$CNA[cnas$ntot <2] = "Deletion"
cnas$cna_id = paste(cnas$CHROM, cnas$Start, cnas$End, sep="_")

cnas_gr = GRanges(
  seqnames = cnas$CHROM,
  ranges = IRanges(cnas$Start, end = cnas$End),
  strand = rep("*", length(cnas$Start)),
  depth_ratio = cnas$depth.ratio,
sample = cnas$Sample,
Nmin = cnas$Nmin,
Nmaj = cnas$Nmaj,
ntot = cnas$ntot,
cna_id = cnas$cna_id)

genes_test = unique(hg19_genes$symbol)

overlap_genes_cnas = function(gene){
  print(gene)
  gene_cords = as.data.table(filter(hg19_genes, symbol==gene))

  #make granges objects
  gene_cords$CHROM = paste("chr", gene_cords$chr, sep="")

  gene_cords_gr = GRanges(
    seqnames = gene_cords$CHROM,
    ranges = IRanges(gene_cords$start, end = gene_cords$end),
    strand = rep("*", length(gene_cords$start)),
    score = 1:length(gene_cords$start))

  #intersect them
  #Then subset the original objects with the negative indices of the overlaps:
  hits <- findOverlaps(gene_cords_gr, cnas_gr, ignore.strand=TRUE)
  hits_overlap = cbind(gene_cords[queryHits(hits),], cnas[subjectHits(hits),])
  #print(head(hits_overlap))
  return(hits_overlap)
}

merge_cnas_in_patient = function(patient){
  print(patient)
  pat_coords = as.data.table(filter(cnas, Patient == patient))

  grl_my = makeGRangesListFromDataFrame(pat_coords,
    split.field ="Sample", seqnames.field = "CHROM",
    start.field = "Start", end.field = "End", keep.extra.columns=TRUE)

  #intersect them
  #Then subset the original objects with the negative indices of the overlaps:
  hits <- findOverlaps(grl_my, ignore.strand=TRUE)
  hits_overlap = cbind(grl_my[queryHits(hits),], cnas[subjectHits(hits),])
  #print(head(hits_overlap))
  return(hits_overlap)
}

#ONLY RUN ONCE

#all_genes_cnas_samples = as.data.table(ldply(llply(genes_test, overlap_genes_cnas, .progress="text")))

#keep only non-neutral CNAs
#all_genes_cnas_samples$Patient = sapply(all_genes_cnas_samples$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
#saveRDS(all_genes_cnas_samples, file="cnas_across_hg19_genomes_all_samples.rds")

all_genes_cnas_samples = readRDS("cnas_across_hg19_genomes_all_samples.rds")

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

dir.create(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))
setwd(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))

#average number of CNAs per sample/patient
t = unique(cnas[,c("Patient", "Sample", "cna_id", "CNA", "Ploidy")])
tt=as.data.table(table(t$Patient, t$Sample))
tt=filter(tt, N >0)
tt %>% group_by(V1) %>% dplyr::summarize(mean = mean(N))
tt %>% group_by(V1) %>% dplyr::summarize(sd = sd(N))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

samples_per_mut = as.data.table(table(t$cna_id, t$Patient))

mut_gene = unique(t[,c("cna_id", "Patient",  "CNA")])

colnames(samples_per_mut) = c("cna_id", "Patient", "num_of_samples_with_mut")
samples_per_mut = merge(samples_per_mut, mut_gene, by=c("cna_id", "Patient"))
samples_per_mut = samples_per_mut[order(-num_of_samples_with_mut)]
samples_per_mut = as.data.table(filter(samples_per_mut, num_of_samples_with_mut >0))

samples_per_mut$phylogeny = ""
z = which(samples_per_mut$Patient == "LY_RAP_0001")
LY_RAP_0001 = samples_per_mut[z,]
LY_RAP_0001$phylogeny[LY_RAP_0001$num_of_samples_with_mut == 3] = "ancestor"
LY_RAP_0001$phylogeny[LY_RAP_0001$num_of_samples_with_mut == 1] = "private"
LY_RAP_0001$phylogeny[LY_RAP_0001$phylogeny == ""] = "shared"

z = which(samples_per_mut$Patient == "LY_RAP_0002")
LY_RAP_0002 = samples_per_mut[z,]
LY_RAP_0002$phylogeny[LY_RAP_0002$num_of_samples_with_mut == 4] = "ancestor"
LY_RAP_0002$phylogeny[LY_RAP_0002$num_of_samples_with_mut == 1] = "private"
LY_RAP_0002$phylogeny[LY_RAP_0002$phylogeny == ""] = "shared"

z = which(samples_per_mut$Patient == "LY_RAP_0003")
LY_RAP_0003 = samples_per_mut[z,]
LY_RAP_0003$phylogeny[LY_RAP_0003$num_of_samples_with_mut == 20] = "ancestor"
LY_RAP_0003$phylogeny[LY_RAP_0003$num_of_samples_with_mut == 1] = "private"
LY_RAP_0003$phylogeny[LY_RAP_0003$phylogeny == ""] = "shared"

samples_per_mut = rbind(LY_RAP_0001, LY_RAP_0002, LY_RAP_0003)

#data table for barplot
barplot = as.data.table(table(samples_per_mut$num_of_samples_with_mut, samples_per_mut$Patient))
barplot = as.data.table(filter(barplot, N >0))
colnames(barplot) = c("num_of_samples_with_mut", "patient", "num_of_muts")
barplot$num_of_samples_with_mut = factor(barplot$num_of_samples_with_mut, levels=unique(barplot$num_of_samples_with_mut))
barplot$patient[barplot$patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
barplot$patient[barplot$patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
barplot$patient[barplot$patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
barplot$patient = factor(barplot$patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

colnames(barplot)[2]="Patient"
pdf("001_samples_per_mutation_lineplot_CNAs.pdf", width=5, height=5)
# Basic barplot
p<-ggline(barplot, x="num_of_samples_with_mut", y="num_of_muts",
palette = c("#00AFBB", "#E7B800", "#FC4E07"), color="Patient")+
xlab("# of samples sharing mutation") + ylab("CNAs count")+
scale_y_continuous(breaks=seq(0,15000,1000))
print(p)
dev.off()

#summarize number of mutations per sample
z = which((cnas$Tissue_Site == "Adrenal gland, NOS") & (cnas$Patient == "LY_RAP_0003"))
cnas$Tissue_Site[z] = "Adrenal gland"
z = which((cnas$Tissue_Site == "Aorta, ascending, not specified \n\n") & (cnas$Patient == "LY_RAP_0001"))
cnas$Tissue_Site[z] = "Aorta, ascending"

muts_per_sample = as.data.table(table(cnas$Patient, cnas$Tissue_Site, cnas$CNA))
muts_per_sample = as.data.table(filter(muts_per_sample, N >0))
muts_per_sample = muts_per_sample[order(-N)]

colnames(muts_per_sample) = c("Patient", "Sample", "type_CNA", "num_of_muts")
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0001"] = "MCL blastoid stage IV"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0002"] = "PMBCL stage IV bulky B symptoms"
muts_per_sample$Patient[muts_per_sample$Patient == "LY_RAP_0003"] = "DLCBL double hit stage IV"
muts_per_sample$Patient = factor(muts_per_sample$Patient, levels=c("MCL blastoid stage IV",
"PMBCL stage IV bulky B symptoms", "DLCBL double hit stage IV"))

#get order of samples for barplot
t = as.data.table(table(cnas$Patient, cnas$Tissue_Site))
t = as.data.table(filter(t, N >0))
t = t[order(-N)]

muts_per_sample$Sample = factor(muts_per_sample$Sample, levels=unique(t$V2))

pdf("002_mutations_per_sample_CNAs.pdf",width=5, height=5)
# Basic barplot
p<-ggplot(data=muts_per_sample, aes(x=Sample, y=num_of_muts, fill=type_CNA)) +
  geom_bar(stat="identity")+theme_classic()+#+ggtitle("Number of mutations per sample")+
	theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5, size=4),
  legend.position = "bottom") + xlab("Sample")+
	ylab("CNAs count")+
	facet_grid(. ~ Patient, scales="free", space='free')+
	theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
	legend.text = element_text(size=6))+
	scale_y_continuous(breaks=seq(0, 4000, by = 200))+
  scale_fill_manual(values=c("#FC4E07", "#00AFBB"))
print(p)
dev.off()
