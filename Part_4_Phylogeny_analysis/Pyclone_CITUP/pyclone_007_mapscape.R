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
#cd /Users/kisaev/UHN/kridel-lab - Documents/
setwd("RAP_WGS/Analysis-Files/Pyclone-vi/All-mutations")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#set up mapscape and test out with limited set of clusters/mutations for now
#later replace with more representative dataset

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#adjaceny matrix (made manually for now) harder to extract from clonevol
citup_ad = read.xlsx("adjacency_matrix_mapscape.xlsx")

#citup_ad[,1] = NULL
#colnames(citup_ad) = c("source", "target")

clone = fread("2020-09-28_clonevol_output_for_mapscape.txt")
clone = as.data.table(filter(clone, !(clone_id==3)))
#reshape so that there are three columns
#sample_id
#clone_id
#clonal_prev

#sample information
samps = fread("/Users/kisaev/UHN/kridel-lab - Documents/RAP_WGS/Data-Files/RAP_samples_information.txt")
samps$sample_id = unique(clone$sample_id)
samps$location_id = unique(clone$sample_id)
samps$x = c(745, 685, 540, 580, 685, 755, 600, 520, 540,
  685, 685, 520, 600, 620, 585, 550, 600, 685, 680, 600)
samps$y = c(495, 510, 310, 800, 750, 450, 850, 820, 350,
  885, 785, 780, 480, 800, 700, 265, 760, 400, 680, 650)
samps = samps[,c("sample_id", "location_id", "x", "y")]

img_ref = "biorender.png"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

mapscape(clone, citup_ad, samps, img_ref)

#mapscape(clonal_prev, tree_edges, sample_locations, img_ref,
#       clone_colours = "NA", mutations = "NA", sample_ids = c("NA"),
#       n_cells = 100, show_low_prev_gtypes = FALSE,
#       phylogeny_title = "Clonal Phylogeny", anatomy_title = "Anatomy",
#       classification_title = "Phylogenetic Classification",
#       show_warnings = TRUE, width = 960, height = 960)

#clonal_prev: ‘data.frame’ Clonal prevalence.  Required columns are:
#sample_id: ‘character()’ id for the tumour sample.
#clone_id: ‘character()’ clone id.
#clonal_prev: ‘numeric()’ clonal prevalence.

#tree_edges: ‘data.frame’ Tree edges of a rooted tree.  Required columns are:
#source: ‘character()’ source node id.
#target: ‘character()’ target node id.

#sample_locations: ‘data.frame’ Anatomic locations for each tumour sample. Required columns are:
#sample_id: ‘character()’ id for the tumour sample.
#location_id: ‘character()’ name of anatomic location for this tumour sample.
#x: ‘numeric()’ x-coordinate (in pixels) for anatomic location on anatomic image.
#y: ‘numeric()’ y-coordinate (in pixels) for anatomic location on anatomic image.

#img_ref: ‘character()’ A reference for the custom anatomical image to
#use, *** in PNG format ***, either a URL to an image hosted
#online or a path to the image in local file system.
