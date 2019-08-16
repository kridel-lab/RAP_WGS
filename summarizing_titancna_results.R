#----------------------------------------------------------------------
#summarize_titancna_results.R
#karin isaev
#last updated: august 16 2019
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries 
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries 
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table", 
              "plyr", "ggpubr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "cowplot")
lapply(packages, require, character.only = TRUE)
library(ggthemes)

date = Sys.Date()

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------


#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

folders  = list.files(pattern="tumor_sample")

get_dat = function(folder){
  all = list.files(folder)
}