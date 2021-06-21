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

  samps$x = c(448, 483, 443)
  samps$y = c(832, 437, 490)
  samps = as.data.frame(samps)
  samps$sample_id = samps$location_id
  colnames(samps) = c("sample_id", "location_id", "x", "y")

  #image of body and tumours
  img_ref = "mapscape_001.png"

  #run tool
  m = mapscape(clone, citup_ad, samps, img_ref, n_cells =300, clone_colours = colors)
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

  samps$x = c(449, 522, 373, 452)
  samps$y = c(835, 445, 806, 733)
  samps = as.data.frame(samps)
  samps$sample_id = samps$location_id
  colnames(samps) = c("sample_id", "location_id", "x", "y")

  #image of body and tumours
  img_ref = "mapscape_002.png"

  #run tool
  m = mapscape(clone, citup_ad, samps, img_ref, n_cells =300, clone_colours = colors)
  return(m)
  print("done yay")
  setwd(dir_main)

}

get_mapscape_patient_3_autopsies = function(patient){

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

  samps = samps[1:17,]
  samps$x = c(645, 371, 558, 470, 410, 549, 454, 361, 456, 557, 556, 380, 431, 547, 459, 451, 533)
  samps$y = c(455, 318, 386, 877, 879, 956, 839, 900, 487, 827, 860, 231, 753, 720, 949, 691, 793)
  samps = as.data.frame(samps)
  samps$sample_id = samps$location_id
  colnames(samps) = c("sample_id", "location_id", "x", "y")

  clone = filter(clone, sample_id %in% samps$sample_id)
  #image of body and tumours
  img_ref = "mapscape_003.png"

  #run tool
  m = mapscape(clone, citup_ad, samps, img_ref, n_cells = 100, clone_colours = colors)
  return(m)
  print("done yay")
  setwd(dir_main)

}

get_mapscape_patient_3_diagnosis = function(patient){

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

  samps = samps[18:20,]
  samps$x = c(550, 389, 620)
  samps$y = c(540, 296, 488)
  samps = as.data.frame(samps)
  samps$sample_id = samps$location_id
  colnames(samps) = c("sample_id", "location_id", "x", "y")

  #image of body and tumours
  clone = filter(clone, sample_id %in% samps$sample_id)
  img_ref = "mapscape_003.png"

  #run tool
  m = mapscape(clone, citup_ad, samps, img_ref, n_cells = 100, clone_colours = colors)
  return(m)
  print("done yay")
  setwd(dir_main)

}

get_mapscape_patient_1(patients[1])
setwd(dir_main)

get_mapscape_patient_2(patients[2])
setwd(dir_main)

get_mapscape_patient_3_autopsies(patients[3])
setwd(dir_main)

get_mapscape_patient_3_diagnosis(patients[3])
setwd(dir_main)
