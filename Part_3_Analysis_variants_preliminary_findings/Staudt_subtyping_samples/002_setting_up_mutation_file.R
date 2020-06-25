#----------------------------------------------------------------------
#001_setting_up_sample_annotation_file.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/LymphGen")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust",
              "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr",
              "ggExtra", "broom", "ggthemes")

lapply(packages, require, character.only = TRUE)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare sample annotation file to be able to run the tool by Staudt
#https://www.sciencedirect.com/science/article/pii/S1535610820301550

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#sample annotation
samps = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/RAP_samples_information.txt")

#structural variant annotation
svs = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_RESULTS/PROCESSED_VCFs/2019-12-16_all_SVs_samples.txt")

#muts
muts = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text/2019-11-27_READ_ONLY_ALL_MERGED_MUTS.txt")

#entrez id conversion
#gene annotations
genes = unique(fread("/cluster/home/kisaev/data/annotables_grch37.txt"))
genes = as.data.table(filter(genes, biotype == "protein_coding"))
genes = as.data.table(filter(genes, !(is.na(entrez))))

#format file needs to be in
#column1 - Sample
#column2 - ENTREZ.ID
#column3 - Type (TRUNC (nonsense coding), MUTATION (missense/frameshift),
        #Synon(5'UTR or synonymous within 4kb of TSS), L265P(MYD88 gene mutation))
#column4 - Location (mutation start position)

dat = unique(muts[,c("Indiv", "Gene.ensGene", "POS", "Func.ensGene",
"ExonicFunc.ensGene",
"AAChange.ensGene", "hg19.ensemblToGeneName.value")])

#add entrez id
colnames(dat)[2] = "ensgene"
dat = merge(dat, genes, by="ensgene")
dat$Type = ""
dat$dist = dat$POS - dat$start

##figure out mutation type for each mutation
#1. L265P(MYD88 gene mutation))
z = which(str_detect(dat$AAChange.ensGene, "L265P"))#none
#2. TRUNC = nonsense coding
dat$Type[dat$ExonicFunc.ensGene == "stopgain"] = "TRUNC"
#3. MUTATION = missense/frameshift
dat$Type[dat$ExonicFunc.ensGene == "stoploss"] = "MUTATION"
dat$Type[dat$ExonicFunc.ensGene == "nonsynonymous_SNV"] = "MUTATION"
#4. Synon(5'UTR or synonymous within 4kb of TSS)
dat$Type[dat$Func.ensGene == "UTR5"] = "Synon"
dat$Type[(dat$Func.ensGene == "downstream") & (abs(dat$dist < 4000))] = "Synon"
dat$Type[(dat$Func.ensGene == "upstream") & (abs(dat$dist < 4000))] = "Synon"

#final columns select
dat = unique(dat[,c("Indiv", "entrez", "Type", "POS")])
dat = as.data.table(filter(dat, !(Type == "")))
colnames(dat)=c("Sample", "ENTREZ.ID", "Type", "Location")

write.table(dat, paste(date, "LymphGen_Sample_Flat_Mutations_File.txt", sep="_"), quote=F, row.names=F, sep="\t")
