#----------------------------------------------------------------------
#exploratory_plotting_002.R
#karin isaev
#last updated: June 24th 2019
#----------------------------------------------------------------------

Sys.setenv("plotly_username"="karini925")
Sys.setenv("plotly_api_key"="pmHjbMwifJL3JzMCNM78")

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("~/Documents/RAP_analysis")

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr", "rjson", 
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom", "EnvStats", "ggthemes")
lapply(packages, require, character.only = TRUE)
library(igraph)
library(randomcoloR)
library(RColorBrewer)
library(ggplotify)
library(plotly)

display.brewer.all()
display.brewer.pal(9, "Set1")
pal = brewer.pal(9, "Set1")
#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarized snvs and cnvs from 21 sequencing folders 
#here, summarize number of mutations/sample/location
#which genes are mutated across all sites which are unique?

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data 
muts = fread("2019-11-27_READ_ONLY_ALL_MERGED_MUTS.txt") 

#Mutation data from WGS Morin et al 2013
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")
tables2 = morin

#PhyloWGS results - one folder for each sample within this folder
folders = list.files("phylowgs_post_witness/")

#-----------------------------------------------------------------------------------
#ANALYSIS 
#-----------------------------------------------------------------------------------

open_phylo_res = function(folder){
  #get sample name 
  sample_name = paste(unlist(strsplit(folder, "_"))[3:8], collapse="_")
  phylo_res_summ = fromJSON(file=paste("phylowgs_post_witness/", folder, "/", sample_name, ".summ.json", sep=""))
  phylo_tree_summary = as.data.frame(phylo_res_summ$trees$`0`$populations)
  
  phylo_res_muts = fromJSON(file=paste("phylowgs_post_witness/", folder, "/", sample_name, ".muts.json", sep=""))
  
  #need to extract zipped tree folder
  phylo_res_tree = fromJSON(file=paste("phylowgs_post_witness/", folder, "/", sample_name, ".mutass", "/0.json", sep="")) 
  phylo_res_tree = phylo_res_tree$mut_assignments
}





#[CELLULAR PREVALENCES OF CLUSTERS ACROSS SAMPLES]

#integrate patient info into tree summary file
phylo_tree_summary = cbind(phylo_tree_summary, ssm_order)
colnames(phylo_tree_summary)[ncol(phylo_tree_summary)] = "Indiv"
phylo_tree_summary = merge(phylo_tree_summary, dna, by = "Indiv")

# id, measure as character/integer/numeric vectors
tree_bar = as.data.table(melt(phylo_tree_summary, id=c("Indiv", "Tissue_Site", "Specimen_Type")))
z = which(str_detect(tree_bar$variable, "cellular"))
tree_bar = tree_bar[z,]
tree_bar$value = log1p(tree_bar$value)
p = ggplot(tree_bar, aes(variable, Indiv , fill = value))+
  geom_tile(color="black") +
  scale_fill_gradient2(low = "yellow", mid = "purple", midpoint = 0.25, high = "steelblue", na.value="transparent") +
  theme_classic()
p + rotate_x_text(90) + rremove("y.ticks")
ggsave("phylowgs_tree_nodes_cellular_prevalences_summary.png")

#[PLOT VAFS OF MUTATIONS IN EACH CLUSTER FOR EACH SAMPLE]

#PARSE MUTATION SUMMARY FILE
get_muts = function(cluster){
  cluster_name = cluster
  z = which(names(phylo_res_tree) == cluster)
  cluster = phylo_res_tree[[z]]
  cluster_muts = unlist(cluster[[1]])
  cluster_muts = as.data.table(filter(ssm_input, id %in% cluster_muts))
  colnames(cluster_muts)[2] = "mut_id"
  cluster_muts = merge(cluster_muts, muts, by = "mut_id")
  cluster_muts$cluster_name = cluster_name
  return(cluster_muts)
}

clusters = names(phylo_res_tree)
all_cluster_mutations = as.data.table(ldply(llply(clusters, get_muts)))

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#make VAF plot
ggplot(data = all_cluster_mutations, aes(x = cluster_name, y = gt_AF, colour=cluster_name)) + 
  stat_summary(fun.data=data_summary, color="blue")+ 
   stat_summary(fun.y=median, geom="point", shape=18, size=3, color="red") +
coord_flip() +facet_wrap(~id.y, ncol=4) #+
  #geom_jitter(position = position_jitter(width = 0.1, height = 0.1))
ggsave("phylowgs_tree_nodes_vaf_summary.png")

#make VAF plot comparing clusters
ggplot(data = all_cluster_mutations, aes(x = id.y, y = gt_AF, colour=id.y)) + 
  stat_summary(fun.data=data_summary, color="blue")+ 
  stat_summary(fun.y=median, geom="point", shape=18, size=3, color="red") +
  coord_flip() +facet_wrap(~cluster_name, ncol=4) #+
#geom_jitter(position = position_jitter(width = 0.1, height = 0.1))
ggsave("phylowgs_tree_nodes_vaf_summary_by_cluster.png")

#WHAT GENES ARE IN EACH CLUSTER? 
t = as.data.table(table(all_cluster_mutations$cluster_name, all_cluster_mutations$hg19.ensemblToGeneName.value))
t=as.data.table(filter(t, N >0))
colnames(t) = c("cluster", "gene", "N")
t=t[order(cluster)]

#summarize which patient has most mutations in each cluster 
pats = as.data.table(table(all_cluster_mutations$id.y, all_cluster_mutations$cluster_name))
colnames(pats) = c("Indiv", "Cluster_Name", "Num_mutations")
tots = as.data.table(table(all_cluster_mutations$cluster_name))
colnames(tots) = c("Cluster_Name", "total_mutations")
pats = merge(pats, tots, by = "Cluster_Name")
pats$prop = pats$Num_mutations/pats$total_mutations

#make barpot 
p <- ggplot(pats, aes(x=Indiv, y=prop, fill=Indiv)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() + facet_wrap(~Cluster_Name, ncol=3) + rremove("x.text") 
p 

p = ggline(pats, x = "Cluster_Name", y = "prop", 
       color = "Indiv") + scale_color_manual(values=distinctColorPalette(length(unique(pats$Indiv)))) 
print(p)
ggsave("phylowgs_tree_nodes_proportion_mutations_in_each_sample.png")
ggplotly(p)
chart_link = api_create(p, filename = "public-graph")





