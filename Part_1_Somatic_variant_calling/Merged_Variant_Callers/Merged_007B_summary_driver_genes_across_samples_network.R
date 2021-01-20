#-------------------------------------------------------------------------------
#Merged_007_summary_driver_genes_across_samples.R
#Karin Isaev
#Monday January 18th, 2020
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")
library(tidyverse)
library(tidygraph)
library(ggraph)
library(tidyverse)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#----------------------------------------------------------------------
#Analysis
#----------------------------------------------------------------------

z=which(read_only$Sample == filter(samps,Tissue_Site =="Kidney, NOS 2")$Indiv)
read_only$Tissue_Site[z] = "Kidney, NOS 2"

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

#p003
p003 = as.data.table(filter(read_only, STUDY_PATIENT_ID == "LY_RAP_0003"))

#prepare edges and nodes
#which mutations in whom
t1= as.data.table(table(p003$mut_id, p003$Tissue_Site)) %>% filter(N >0)
t1$N=NULL
colnames(t1) = c("mutation", "sample")

founders = as.data.table(table(t1$mutation)) %>% filter(N <20)

pairs = as.data.table(t1 %>%
    mutate(n = 1) %>%
    spread(sample, n, fill=0) %>%
    select(-mutation) %>%
    {crossprod(as.matrix(.))} %>%
    replace(lower.tri(., diag=T), NA) %>%
    reshape2::melt(na.rm=T) %>%
    unite('sample', c('Var1', 'Var2'), sep=":")) %>%
     separate(col = "sample",into = c("from", "to"), sep = ":")

#convert to edges object, "from" "to"
edges=pairs

nodes = as.data.table(unique(c(edges$from, edges$to)))
colnames(nodes)="id"
nodes$type="FF"
nodes$type[] = "FFPE"
#nodes$label = nodes$id
#nodes$id = 1:20

p003.net <- tbl_graph(
  nodes = nodes,
  edges = edges,
  directed = TRUE)

p003.net = p003.net %>%
  activate(edges) %>%
  rename(weight = value)

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/test.pdf")
set.seed(123)
ggraph(p003.net) +
  geom_edge_link(aes(width = weight), alpha = 0.2) +
  scale_edge_width(range = c(0.2, 1)) +
  geom_node_point(aes(color = type), size = 2) +
  geom_node_text(aes(label = id), size = 3, repel = TRUE)

p003.net %>%
    activate(nodes) %>%
    mutate(centrality = centrality_authority()) %>%
    ggraph(layout = "graphopt") +
    geom_edge_link(aes(width = weight), alpha = 0.2) +
    #geom_edge_link(width = 1, colour = "lightgray") +
    geom_node_point(aes(size = centrality, colour = centrality)) +
    geom_node_text(aes(label = id), repel = TRUE)+
    scale_color_gradient(low = "yellow", high = "red")

dev.off()
