#----------------------------------------------------------------------
#sciclone_001_make_input_files.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()
print(date)

options(stringsAsFactors=F)
setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone")

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom")
lapply(packages, require, character.only = TRUE)
library(RColorBrewer)
library(openxlsx)
library(plotly)
library(readxl)
library(GenomicRanges)
library(sciClone)

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#here try running using full mutation matrix produced via bamreadcount
#which should include read counts for mutations not called in cases
#try running sciclone on all samples and all mutations to see if it works
#otherwise might need to reduce number of mutations

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. SNV list
snv_list = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_snvs_df.rds")

#2. Copy Number List
cnas_list = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_cnas_df.rds")

#3. sample names (use ids for now)
samples = names(snv_list)

#4. annotation file
annotations = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/annotation_file.rds")

#5. loh regions
loh = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone/list_loh_df.rds")
loh = unique(as.data.table(ldply(loh)))
loh = unique(loh[,-1])
#loh$CHROM=sapply(loh$CHROM, function(x){unlist(strsplit(x, "chr"))[2]})

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#------------------------------------
#clustering using all samples:
pairs_done = c("0")

for(i in 1:20){
  for(j in 1:20){
      if(!(i==j)){
      print(i)
      print(j)
      if(!(paste(i,j,sep="_")) %in% pairs_done){

      sc = sciClone(vafs=snv_list[i:j],
       #copyNumberCalls=cnas_list[i:j],
       annotation=annotations,
       #regionsToExclude=loh,
        sampleNames=samples[i:j], maximumClusters=10,
        #copyNumberMargins=0.95,
        minimumDepth=30)
        #create output
        p = paste(samples[i:j], collapse="_")

        #res file
        writeClusterTable(sc, paste(p, "sciclone_table.txt", sep=""))

        #plot
        sc.plot2d(sc, paste(p, ".pdf", sep=""))

        pairs_done = c(pairs_done, paste(i,j,sep="_"))
        pairs_done = c(pairs_done, paste(j,i,sep="_"))
  }
}
}
}

print("done")

#------------------------------------
