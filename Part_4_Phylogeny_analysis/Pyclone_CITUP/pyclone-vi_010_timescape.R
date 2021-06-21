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
library(timescape)

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

get_mapscape_patient_1 = function(patient){

  setwd(patient)

  colors = read.xlsx("clones_colours.xlsx")

  #adjaceny matrix (made manually for now) harder to extract from clonevol
  citup_ad = read.xlsx("RAP_001_adjacency_matrix.xlsx")

  #clone and sample data
  clone = fread("mapscape_input_from_pairtree_pyclonevi.txt")

  #sample information
  z = which(sample_data$Indiv %in% clone$sample_id)
  samps = sample_data[z,c("Indiv", "Confirmed_by_Robert")]
  colnames(samps) = c("sample_id", "location_id")
  samps = as.data.frame(samps)

  clone = merge(clone, samps, by = "sample_id")
  clone = as.data.frame(clone)
  clone$sample_id = clone$location_id
  clone$location_id = NULL
  colnames(clone)[1] = "timepoint"

  #run tool
  m = timescape(clone, citup_ad, clone_colours = colors)
  return(m)

  print("done yay")
  setwd(dir_main)
}

get_mapscape_patient_2 = function(patient){

  setwd(patient)

  colors = read.xlsx("clones_colours.xlsx")

  #adjaceny matrix (made manually for now) harder to extract from clonevol
  citup_ad = read.xlsx("RAP_002_adjacency_matrix.xlsx")

  #clone and sample data
  clone = fread("mapscape_input_from_pairtree_pyclonevi.txt")

  #sample information
  z = which(sample_data$Indiv %in% clone$sample_id)
  samps = sample_data[z,c("Indiv", "Confirmed_by_Robert")]
  colnames(samps) = c("sample_id", "location_id")
  samps = as.data.frame(samps)
  samps$location_id[samps$location_id == "Adrenal, right"] = "Adrenal"

  clone = merge(clone, samps, by = "sample_id")
  clone = as.data.frame(clone)
  clone$sample_id = clone$location_id
  clone$location_id = NULL
  colnames(clone)[1] = "timepoint"

  #run tool
  m = timescape(clone, citup_ad, clone_colours = colors)
  return(m)
  print("done yay")
  setwd(dir_main)

}

get_mapscape_patient_3 = function(patient){

  setwd(patient)

  colors = read.xlsx("clones_colours.xlsx")

  #adjaceny matrix (made manually for now) harder to extract from clonevol
  citup_ad = read.xlsx("RAP_003_adjacency_matrix.xlsx")

  #clone and sample data
  clone = fread("mapscape_input_from_pairtree_pyclonevi.txt")

  #sample information
  z = which(sample_data$Indiv %in% clone$sample_id)
  samps = sample_data[z,c("Indiv", "New_Tissue_Site")]
  colnames(samps) = c("sample_id", "location_id")
  samps = as.data.frame(samps)

  clone = merge(clone, samps, by = "sample_id")
  clone = as.data.frame(clone)
  clone$sample_id = clone$location_id
  clone$location_id = NULL
  colnames(clone)[1] = "timepoint"

  #run tool
  m = timescape(clone, citup_ad, clone_colours = colors)
  return(m)
  print("done yay")
  setwd(dir_main)

}

get_mapscape_patient_1(patients[1])
setwd(dir_main)

get_mapscape_patient_2(patients[2])
setwd(dir_main)

get_mapscape_patient_3(patients[3])
setwd(dir_main)
