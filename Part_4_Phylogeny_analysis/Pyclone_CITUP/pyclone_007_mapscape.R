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
library(mapscape)

date = Sys.Date()
setwd("/Users/kisaev/UHN/kridel-lab - Documents/RAP_WGS/Analysis-Files/Pyclone")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#set up mapscape and test out with limited set of clusters/mutations for now
#later replace with more representative dataset

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#adjaceny matrix
citup_ad = fread("pat003_adjacency.csv")
citup_ad[,1] = NULL
colnames(citup_ad) = c("source", "target")

clone = fread("pat003_clone.csv")
#reshape so that there are three columns
#sample_id
#clone_id
#clonal_prev
colnames(clone) = c("sample_id", "0", "1", "2", "3", "4")
clone = melt(clone, id = "sample_id")
colnames(clone) = c("sample_id", "clone_id", "clonal_prev")

#sample information
samps = fread("/Users/kisaev/UHN/kridel-lab - Documents/RAP_WGS/Data-Files/RAP_samples_information.txt")
samps$sample_id = unique(clone$sample_id)
samps$location_id = samps$Tissue_Site
samps$x = c(1,2,3,4,5,6,7,8,9,10,11,12,11,8,9,3,4,5,6,7)
samps$y = c(1,2,3,4,5,6,7,8,9,10,11,12,11,8,9,3,4,5,6,7)
samps = samps[,c("sample_id", "location_id", "x", "y")]

img_ref = "pic1-68.png"

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

mapscape(clone, citup_ad, samps, img_ref, ncells=20)

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
