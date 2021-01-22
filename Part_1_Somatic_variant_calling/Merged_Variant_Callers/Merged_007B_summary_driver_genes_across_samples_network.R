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

#+------------------------------------------------------------------------------
#p003
p003 = as.data.table(filter(read_only, STUDY_PATIENT_ID == "LY_RAP_0003"))
#+------------------------------------------------------------------------------

#prepare edges and nodes
#which mutations in whom
t1= as.data.table(table(p003$mut_id, p003$Tissue_Site)) %>% filter(N >0)
t1$N=NULL
colnames(t1) = c("mutation", "sample")

founders = as.data.table(table(t1$mutation)) %>% filter(N ==20)
unique = as.data.table(table(t1$mutation)) %>% filter(N ==1)
t1=as.data.table(filter(t1, !(mutation %in% founders$V1), !(mutation %in% unique$V1)))

pairs = as.data.table(t1 %>%
    mutate(n = 1) %>%
    spread(sample, n, fill=0) %>%
    select(-mutation) %>%
    {crossprod(as.matrix(.))} %>%
    replace(lower.tri(., diag=T), NA) %>%
    reshape2::melt(na.rm=T) %>%
    unite('sample', c('Var1', 'Var2'), sep=":")) %>%
     separate(col = "sample",into = c("from", "to"), sep = ":")
pairs=pairs[order(-value)]
pairs$from=factor(pairs$from, levels=unique(pairs$from))
pairs$to=factor(pairs$to, levels=unique(pairs$to))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/p003_summary_shared_mutations.pdf")
#geom_tile plot
p = ggplot(pairs, aes(from, to)) +
  geom_tile(aes(fill = value), colour = "grey50") +
  xlab("From") + ylab("To") +
	ggtitle("How many mutations each pair of samples share \nnot incluidng founders")+theme_bw()+
  #scale_fill_gradient(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_gradientn(colours = c("blue", "lightblue", "yellow", "pink", "red"),
                       breaks=c(3000,4000,5000,6000, 7000, Inf),
                       na.value = "red")
print(p)
dev.off()

#+------------------------------------------------------------------------------
#p002
p002 = as.data.table(filter(read_only, STUDY_PATIENT_ID == "LY_RAP_0002"))
#+------------------------------------------------------------------------------

#prepare edges and nodes
#which mutations in whom
t1= as.data.table(table(p002$mut_id, p002$Tissue_Site)) %>% filter(N >0)
t1$N=NULL
colnames(t1) = c("mutation", "sample")

founders = as.data.table(table(t1$mutation)) %>% filter(N ==4)
unique = as.data.table(table(t1$mutation)) %>% filter(N ==1)
t1=as.data.table(filter(t1, !(mutation %in% founders$V1), !(mutation %in% unique$V1)))

pairs = as.data.table(t1 %>%
    mutate(n = 1) %>%
    spread(sample, n, fill=0) %>%
    select(-mutation) %>%
    {crossprod(as.matrix(.))} %>%
    replace(lower.tri(., diag=T), NA) %>%
    reshape2::melt(na.rm=T) %>%
    unite('sample', c('Var1', 'Var2'), sep=":")) %>%
     separate(col = "sample",into = c("from", "to"), sep = ":")
pairs=pairs[order(-value)]
pairs$from=factor(pairs$from, levels=unique(pairs$from))
pairs$to=factor(pairs$to, levels=unique(pairs$to))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/p002_summary_shared_mutations.pdf")
#geom_tile plot
p = ggplot(pairs, aes(from, to)) +
  geom_tile(aes(fill = value), colour = "grey50") +
  xlab("From") + ylab("To") +
	ggtitle("How many mutations each pair of samples share \nnot incluidng founders")+theme_bw()+
  #scale_fill_gradient(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_fill_gradientn(colours = c("blue", "lightblue", "yellow", "pink", "red"),
                       breaks=c(19000,20000,21000,22000, 23000, Inf),
                       na.value = "red")
print(p)
dev.off()

#+------------------------------------------------------------------------------
#p001
p001 = as.data.table(filter(read_only, STUDY_PATIENT_ID == "LY_RAP_0001"))
#+------------------------------------------------------------------------------

#prepare edges and nodes
#which mutations in whom
t1= as.data.table(table(p001$mut_id, p001$Tissue_Site)) %>% filter(N >0)
t1$N=NULL
colnames(t1) = c("mutation", "sample")

founders = as.data.table(table(t1$mutation)) %>% filter(N ==3)
unique = as.data.table(table(t1$mutation)) %>% filter(N ==1)
t1=as.data.table(filter(t1, !(mutation %in% founders$V1), !(mutation %in% unique$V1)))

pairs = as.data.table(t1 %>%
    mutate(n = 1) %>%
    spread(sample, n, fill=0) %>%
    select(-mutation) %>%
    {crossprod(as.matrix(.))} %>%
    replace(lower.tri(., diag=T), NA) %>%
    reshape2::melt(na.rm=T) %>%
    unite('sample', c('Var1', 'Var2'), sep=":")) %>%
     separate(col = "sample",into = c("from", "to"), sep = ":")
pairs=pairs[order(-value)]
pairs$from=factor(pairs$from, levels=unique(pairs$from))
pairs$to=factor(pairs$to, levels=unique(pairs$to))

pdf("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/p001_summary_shared_mutations.pdf")
#geom_tile plot
p = ggplot(pairs, aes(from, to)) +
  geom_tile(aes(fill = value), colour = "grey50") +
  xlab("From") + ylab("To") +
	ggtitle("How many mutations each pair of samples share \nnot incluidng founders")+theme_bw()+
  scale_fill_gradient(low="blue", high="red")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))#+
  #scale_fill_gradientn(colours = c("blue", "lightblue", "yellow", "pink", "red"),
  #                     breaks=c(300000,300100,300200,300300, 300400, Inf),
#                       na.value = "red")
print(p)
dev.off()
