#----------------------------------------------------------------------
#exploratory_plotting_002.R
#karin isaev
#last updated: June 24th 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

date = Sys.Date()

options(stringsAsFactors=F)
setwd("~/Documents/RAP_analysis")

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

#generate intermutation plot 

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data 
muts = fread("2019-11-27_READ_ONLY_ALL_MERGED_MUTS.txt") 
morin = read.xlsx("supp_blood-2013-02-483727_TableS3.xlsx")

#----------------------------------------------------------------------
#analysis 
#----------------------------------------------------------------------

#for each chromosome, order mutations and calculate intermutation distance 

get_dists_full = function(sample){
  
  samp_dat = as.data.table(filter(muts, Indiv == sample))
  samp_dat$chr_n = sapply(samp_dat$CHROM, function(x){unlist(strsplit(x, "chr"))[2]})
  samp_dat$chr_n = as.numeric(samp_dat$chr_n)
  samp_dat = samp_dat[order(chr_n)]
  chrs = unique(samp_dat$CHROM)
  samp_dat$mut_type = paste(samp_dat$REF, samp_dat$ALT, sep=">")
  
  get_dists = function(chr){
    samp_dat_chr = as.data.table(filter(samp_dat, CHROM == chr))
    samp_dat_chr = samp_dat_chr[order(POS)]
    test = unique(samp_dat_chr[,c("CHROM", "POS", "Indiv", "mut_type", "Func.ensGene")])
    #get intermutation distances 
    get_distances = as.data.table(test %>% mutate(Diff = POS - lag(POS)))
    return(get_distances)
  }
  
  all_chrs_dists = as.data.table(ldply(llply(chrs, get_dists)))
  all_chrs_dists$mut_num = 1:nrow(all_chrs_dists)
  all_chrs_dists$mb = all_chrs_dists$Diff/1000
  
  #make plot
  p = ggplot(all_chrs_dists, aes(mut_num, mb))
  p = p + geom_point(aes(colour = factor(Func.ensGene))) + theme_base() + ylab("log10 KB") + geom_hline(yintercept=37.529, linetype="dashed", color = "red")
  print(p + scale_y_continuous(trans='log10') + ggtitle(samp_dat$id[samp_dat$Indiv == sample][1]))
  return(all_chrs_dists)
  }

pdf(paste(date, "intermutation_distances_plots_all_samples.pdf"), width=12, height=9)
all_muts_dats = as.data.table(ldply(llply(unique(muts$Indiv), get_dists_full, .progress="text")))
dev.off()


