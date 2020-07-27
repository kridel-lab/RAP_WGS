#----------------------------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

#pyclone summary

date = Sys.Date()

options(stringsAsFactors=F)
setwd("/Users/kisaev/OneDrive - UHN/RAP_ANALYSIS")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats", "ggthemes")
lapply(packages, require, character.only = TRUE)

library(RColorBrewer)
library(openxlsx)
library(plotly)

display.brewer.all()
display.brewer.pal(9, "Set1")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize pyclone run

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data
muts = fread("2019-11-27_READ_ONLY_ALL_MERGED_MUTS.txt")
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")
genes_class = read.xlsx("Figure2_data_ABC_GCB_gene_Dave.xlsx")
colnames(genes_class)[1] = "hg19.ensemblToGeneName.value"
#sample info
samp_info = fread("RAP_samples_information.txt")

#2. pyclone data
setwd("/Users/kisaev/OneDrive - UHN/RAP_ANALYSIS/Pyclone")
list.files()

#3. concatenate all pyclone inputs into one file
inputs=list.files(pattern="pyclone_input.tsv")
all_inputs= as.data.table(ldply(llply(inputs, function(x){fread(x)})))

#4. pyclone res
pyc_cluster=fread("RAP_WGS_pyclone_table_file_cluster.txt")
pyc_table = fread("RAP_WGS_pyclone_table_file.txt")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#heatmap of clusters and mean prevalence per sample
pyc_cluster$Indiv = sapply(pyc_cluster$sample_id, function(x){paste(unlist(strsplit(x, "_"))[2:7], collapse="_")})
pyc_cluster = as.data.table(filter(pyc_cluster, size >=2))
pyc_cluster$cluster_id = paste("cluster", pyc_cluster$cluster_id, sep="_")
pyc_cluster = merge(pyc_cluster, samp_info, by="Indiv")

library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)

set.seed(01115)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
col=sample(col_vector, n)

# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
#col = list(site = c("FFPE" = "green", "FT" = "gray"),
#           organ =c("Axilla, NOS" = col_vector[1], "Cervical lymph node"= col_vector[2], "Shoulder, NOS"= col_vector[3], "Omentum"= col_vector[4],
#                    "Abdomen, NOS"= col_vector[5], "Inguinal region, NOS"= col_vector[6],
#           "Retroperitoneum, NOS"= col_vector[7], "Cecum"= col_vector[8],"Mediastinum, NOS"= col_vector[9] ,"Kidney, NOS"= col_vector[10], "Parotid gland"= col_vector[11],
#           "Pancreas, NOS"= col_vector[12],
#           "Spleen"= col_vector[13] , "Bladder, NOS"= col_vector[14], "Stomach, NOS"= col_vector[15],
#           "Adrenal gland, NOS"= col_vector[16], "right_neck_LN"= col_vector[17] , "left_axilla_LN" = col_vector[18],
#           "left_breast"= col_vector[19]))

matrix_heatmtap = as.data.frame(dcast(pyc_cluster, cluster_id ~ id, value.var = "mean"))
rownames(matrix_heatmtap) = matrix_heatmtap$cluster_id
matrix_heatmtap$cluster_id = NULL
matrix_heatmtap = as.matrix(matrix_heatmtap)
theme_update(text = element_text(size=6))
p=Heatmap(matrix_heatmtap, name = "RAP_WGS", column_names_gp = gpar(fontsize = 8),
        row_names_gp = gpar(fontsize = 8), rect_gp = gpar(col = "black", lwd = 1))
p
pdf(paste(date, "Pyclone_results_all_samples_clusters.pdf", sep="_"), width=10)
print(p)
dev.off()

#add mutation information
colnames(pyc_table)[1] = "mut_id"
muts = unique(muts[,c("mut_id", "Gene.ensGene", "ChromKey", "POS", "CHROM", "Func.ensGene" ,"GeneDetail.ensGene" ,"ExonicFunc.ensGene", "AAChange.ensGene", "hg19.ensemblToGeneName.value")])
pyc_table = merge(pyc_table, muts, by = "mut_id")
pyc_table$cluster_id =  paste("cluster", pyc_table$cluster_id, sep="_")
pyc_table = as.data.table(filter(pyc_table, cluster_id %in% pyc_cluster$cluster_id))
z = which(pyc_table$hg19.ensemblToGeneName.value  %in% morin$Gene)
pyc_table$morin = ""
pyc_table$morin[z] = "morin_WGS"
write.table(pyc_table, paste(date, "pyclone_clusters_results_wGene_names_mutations.csv", sep="_"), quote=F, row.names=F, sep=";")
