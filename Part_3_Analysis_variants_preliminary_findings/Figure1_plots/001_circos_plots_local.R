#-------------------------------------------------------------------------------
#Make circos plots for each patient
#Karin Isaev
#-------------------------------------------------------------------------------

#load packages and data
library(dplyr)
library(ggpubr)
library(data.table)
library(BioCircos)
library(plyr)

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

snvs = readRDS("SNVs/Mutect2_Strelka_merged_mutations_wCNA_status/2021-04-06_Mutect2_Strelka_merged_mutations_wCNA_status.rds")
svs = readRDS("SVs/2021-04-20_RAP_WGS_all_SVs_heatmap_plot.rds")
nontransloc = fread("SVs/2021-05-17_RAP_WGS_nontranslocations_ONLY.txt", sep="}")
cnas = readRDS("CNAs/all_CNAs_by_Sequenza.rds")

cnas$Patient = sapply(cnas$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
cnas$cna_id = paste(cnas$CHROM, cnas$Start, cnas$End, sep="_")
cnas$width = cnas$End - cnas$Start

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

dir.create(file.path("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/circos_plots"))
setwd(file.path("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/circos_plots"))

patient = "LY_RAP_0001"

#1. SNVs

generate_circos = function(patient){

  #1. SNPs keep

  snps = filter(snvs, STUDY_PATIENT_ID == patient, biotype == "protein_coding",
  Func.ensGene %in% c("exonic", "splicing", "exonic\\x3bsplicing", "UTR5", "UTR3"))

  #keep only ancestral ones

  if(patient == "LY_RAP_0001"){
    num_check = 3
  }

  if(patient == "LY_RAP_0002"){
    num_check = 4
  }

  if(patient == "LY_RAP_0003"){
    num_check = 20
  }

  t=as.data.table(table(snps$mut_id)) %>% filter(N == num_check)
  snps = filter(snps, mut_id %in% t$V1)

  # Chromosomes on which the points should be displayed
  points_chromosomes = unique(snps$CHROM)
  points_chromosomes = (sapply(points_chromosomes, function(x){unlist(strsplit(x, "chr"))[2]}))

  # Coordinates on which the points should be displayed
  points_coordinates = unique(snps$POS)

  #2. SVs keep

  svs_keep = filter(svs, STUDY_PATIENT_ID == patient, SVTYPE == "BND")
  t=as.data.table(table(svs_keep$id)) %>% filter(N ==num_check)
  svs_keep = filter(svs_keep, id %in% t$V1)
  svs_keep = unique(svs_keep[,c("STUDY_PATIENT_ID", "gene", "id")])
  links_chromosomes_1 = sapply(svs_keep$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[2], "_"))[1]})
  links_chromosomes_2 = sapply(svs_keep$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[3], "_"))[1]})

  links_pos_1 = as.numeric(sapply(svs_keep$id, function(x){unlist(strsplit(x, "_"))[2]}))
  links_pos_2 = as.numeric(sapply(svs_keep$id, function(x){unlist(strsplit(x, "_"))[5]}))

  links_labels = as.character(svs_keep$gene)

  #3. CNAs keep

  cnas_pat = filter(cnas, Patient == patient, ntot > 2 | ntot < 2)
  t=as.data.table(table(cnas_pat$cna_id)) %>% filter(N ==num_check)
  cnas_pat = filter(cnas_pat, cna_id %in% t$V1)
  cnas_pat = unique(cnas_pat[,c("cna_id", "Patient", "ntot", "width")])
  z = which(duplicated(cnas_pat$cna_id))
  if(!(length(z) == 0)){
    cnas_pat = cnas_pat[-z,]
  }

  # Arcs coordinates
  arcs_chromosomes = sapply(cnas_pat$cna_id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[2], "_"))[1]})
  arcs_begin = as.numeric(sapply(cnas_pat$cna_id, function(x){unlist(strsplit(x, "_"))[2]}))
  arcs_end = as.numeric(sapply(cnas_pat$cna_id, function(x){unlist(strsplit(x, "_"))[3]}))

  #Make tracks for circos plot

  #first SNPs
  points_values = 0:1

  tracklist = BioCircosSNPTrack('mySNPTrack', points_chromosomes, points_coordinates,
    points_values, colors = c("tomato2", "darkblue"), minRadius = 0.8, maxRadius = 0.95)

  # Background are always placed below other tracks
  tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack",
    minRadius = 0.8, maxRadius = 0.95,
    borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF")

  #CNAs
  tracklist = tracklist + BioCircosArcTrack('myArcTrack', arcs_chromosomes, arcs_begin, arcs_end,
      minRadius = 0.68, maxRadius = 0.78, colors=c("green"))

  tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', links_chromosomes_1, links_pos_1,
      links_pos_1 + 50000000, links_chromosomes_2, links_pos_2, links_pos_2 + 750000,
      maxRadius = 0.65, color = "orange", labels = links_labels, labelSize=0.2, labelColor="blue", displayLabel=F)

  tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.65,
          borderSize = 0, fillColors = "#EEFFEE")

  BioCircos(tracklist, genomeFillColor = "PuOr",
    displayGenomeBorder = FALSE, yChr =  FALSE,
    chrPad = 0,
       genomeTicksLen = 3, genomeTicksTextSize = 0,
       genomeTicksScale = 40000000,
       genomeLabelTextSize = 18, genomeLabelDy = 0)

}
