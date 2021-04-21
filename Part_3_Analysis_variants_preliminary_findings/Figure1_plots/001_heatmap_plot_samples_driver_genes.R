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

#summarize driver genes and what kind of mutations they are affected by
drivers_all = filter(samples_per_mut, driver=="driver")
drivers_all = as.data.table(table(drivers_all$patient, drivers_all$symbol, drivers_all$phylogeny, drivers_all$Func.ensGene)) %>% filter(N > 0, V4 %in% c("exonic", "intronic", "splicing", "UTR3", "UTR5"))
colnames(drivers_all) = c("Patient", "Gene", "Status", "Mutation", "Mutation_count")
drivers_all$Status = factor(drivers_all$Status, levels=c("ancestor", "shared", "private"))

#keep only functional protein coding gene mutations
samples_per_mut = as.data.table(filter(samples_per_mut, biotype == "protein_coding",
driver=="driver", Func.ensGene %in% c("exonic", "splicing"),
!(ExonicFunc.ensGene == "synonymous_SNV"), phylogeny == "ancestor"))

#get mutation info for each sample and VAF/copy number info
samples_per_mut$combo = paste(samples_per_mut$mut_id, samples_per_mut$patient, sep="_")
read_only$combo = paste(read_only$mut_id, read_only$STUDY_PATIENT_ID, sep="_")

driver_ancestor_muts = filter(read_only, combo %in% samples_per_mut$combo)
driver_ancestor_muts$gene_mut = paste(driver_ancestor_muts$symbol, driver_ancestor_muts$POS, sep="_")

#MAKE HEATMAP OF HIGH CONFIDENCE ANCESTRAL MUTATIONS
dir.create(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))
setwd(file.path("/cluster/projects/kridelgroup/RAP_ANALYSIS/plots", date))

#order genes by most mutationsa across all samples
#order samples with patients by decreasing number of mutations 

pdf("FIGURE1_heatmap_high_conf_ancestral_muts.pdf", width=7, height=4)
p = ggplot(driver_ancestor_muts, aes(gene_mut, Tissue_Site)) +
  geom_tile(aes(fill = gt_AF), colour = "grey50") +
  theme_bw()+
  xlab("Gene mutation") + ylab("Sample") +
	facet_grid(STUDY_PATIENT_ID ~ ., scales="free", space="free")+
	theme(axis.text=element_text(size=6), axis.text.x = element_text(angle = 90, size=5, vjust = 0.5, hjust=1))+
   scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
print(p)
dev.off()
