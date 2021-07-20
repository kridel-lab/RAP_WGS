#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(dplyr)

#-------------------------------------------------------------------------------
#Purpose
#-------------------------------------------------------------------------------

#for each patient see if can identify tumour clonal CNA events
#these will have same major minor allele calls across all clones

setwd("/Users/kisaev/Documents/Hatchet_analysis/p001/results")

#-------------------------------------------------------------------------------
#Analysis
#-------------------------------------------------------------------------------

#test file
hatchet_res = fread("chosen.tetraploid.bbc.ucn")

#patient specific
pats = unique(hatchet_res$SAMPLE)

#function define
get_clonal_cnas = function(pat){

  pat_cnas = filter(hatchet_res, SAMPLE==pat)

  check_cnas= transform(pat_cnas,
  same = apply(pat_cnas[,c("cn_clone1", "cn_clone2", "cn_clone3", "cn_clone4")], 1, function(x) length(unique(x)) == 1))

}
