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
