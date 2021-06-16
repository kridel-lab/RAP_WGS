#----------------------------------------------------------------------
#pyclone_007_mapscape.R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "cowplot")
lapply(packages, require, character.only = TRUE)
library("readxl")
library(openxlsx)
library(mapscape)

date = Sys.Date()

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Analysis-Files/Pairtree/Mapscape")
dir_main = "/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Analysis-Files/Pairtree/Mapscape"

sample_data = read.xlsx("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files/27_RAP_samples_used_in_study.xlsx")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#set up mapscape and test out with limited set of clusters/mutations for now
#later replace with more representative dataset

patients = c("LY_RAP_0001", "LY_RAP_0002", "LY_RAP_0003")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

get_mapscape = function(patient){

  setwd(patient)

  #adjaceny matrix (made manually for now) harder to extract from clonevol
  citup_ad = read.xlsx("RAP_001_adjacency_matrix.xlsx")

  #clone and sample data
  clone = fread("mapscape_input_from_pairtree_pyclonevi.txt")

  #sample information
  z = which(sample_data$Indiv %in% clone$sample_id)
  samps = sample_data[z,c("Indiv", "Confirmed_by_Robert")]
  colnames(samps) = c("sample_id", " location_id")
  samps$x = c(745, 685, 540)
  samps$y = c(495, 510, 310)

  #image of body and tumours
  img_ref = "mapscape_001.png"

  #run tool
  mapscape(clone, citup_ad, samps, img_ref)

  setwd(dir_main)

}
