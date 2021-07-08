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
library(Signac)
library(GenomicRanges)

setwd("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/Data-Files")

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#patient <- commandArgs(trailingOnly = TRUE)
#patient_id = patient[1]
#print(patient_id)

snvs = readRDS("SNVs/Mutect2_Strelka_merged_mutations_wCNA_status/2021-04-06_Mutect2_Strelka_merged_mutations_wCNA_status.rds")
svs = readRDS("SVs/2021-04-20_RAP_WGS_all_SVs_heatmap_plot.rds")
nontransloc = fread("SVs/2021-05-17_RAP_WGS_nontranslocations_ONLY.txt", sep="}")
cnas = readRDS("CNAs/all_CNAs_by_Sequenza.rds")

cnas$Patient = sapply(cnas$Sample, function(x){paste(unlist(strsplit(x, "_"))[1:3], collapse="_")})
cnas$cna_id = paste(cnas$CHROM, cnas$Start, cnas$End, sep="_")
cnas$width = cnas$End - cnas$Start

#function to merge CNAs
UnifyPeaks = function (object.list, mode = "reduce"){
    peak.ranges <- list()
    for (i in seq_along(along.with = object.list)) {
        peak.ranges[[i]] <- object.list[[i]]
    }
    peak.ranges <- Reduce(f = c, x = peak.ranges)
    if (mode == "reduce") {
        return(reduce(x = peak.ranges))
    }
    else if (mode == "disjoin") {
        return(disjoin(x = peak.ranges))
    }
    else {
        stop("Unknown mode requested")
    }
}

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#summarize mutation patterns across samples and driver genes
#check which mutations occur in all samples versus only 1 or several

#dir.create(file.path("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/circos_plots"))
setwd(file.path("/Users/kisaev/UHN/kridel-lab - Documents (1)/RAP_WGS/circos_plots"))

#patient = "LY_RAP_0001"

#1. SNVs

generate_circos = function(patient){

  #1. SNPs keep
  print(paste(patient, "started"))

  snps = filter(snvs, STUDY_PATIENT_ID == patient, biotype == "protein_coding",
  (Func.ensGene %in% c("exonic", "splicing", "exonic\\x3bsplicing") &
  (ExonicFunc.ensGene %in% c("nonsynonymous_SNV", "stopgain",
 "frameshift_deletion", "frameshift_insertion", "stoploss", "."))))

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
  snps = unique(snps[,c("CHROM", "POS")])

  # Chromosomes on which the points should be displayed
  points_chromosomes = snps$CHROM
  points_chromosomes = sapply(points_chromosomes, function(x){unlist(strsplit(x, "chr"))[2]})

  # Coordinates on which the points should be displayed
  points_coordinates = snps$POS

  #2. SVs keep

  svs_keep = filter(svs, STUDY_PATIENT_ID == patient, SVTYPE == "BND")
  svs_keep$chr1 = sapply(svs_keep$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[2], "_"))[1]})
  svs_keep$chr2 = sapply(svs_keep$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[3], "_"))[1]})

  z = which(svs_keep$chr1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
"13", "14", "15", "16", "17", "18", "19", "20"))
  svs_keep = svs_keep[z,]
  z = which(svs_keep$chr2 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
"13", "14", "15", "16", "17", "18", "19", "20"))
  svs_keep = svs_keep[z,]

  svs_keep$gene1 = sapply(svs_keep$gene, function(x){unlist(strsplit(as.character(x), "_"))[1]})
  svs_keep$gene2 = sapply(svs_keep$gene, function(x){unlist(strsplit(as.character(x), "_"))[2]})

  svs_keep$new_id = paste(svs_keep$chr1, svs_keep$chr2, svs_keep$gene2, sep="_")
  z = which(svs_keep$gene2 == "NA")
  if(!(length(z)==0)){
  svs_keep$new_id[z] = paste(svs_keep$chr1[z], svs_keep$chr2[z], svs_keep$gene1[z], sep="_")}

  t=as.data.table(table(svs_keep$new_id)) %>% filter(N >= num_check)
  svs_keep = filter(svs_keep, new_id %in% t$V1)
  svs_keep = unique(svs_keep[,c("STUDY_PATIENT_ID", "gene", "id", "new_id")])
  links_chromosomes_1 = sapply(svs_keep$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[2], "_"))[1]})
  links_chromosomes_2 = sapply(svs_keep$id, function(x){unlist(strsplit(unlist(strsplit(x, "chr"))[3], "_"))[1]})

  links_pos_1 = as.numeric(sapply(svs_keep$id, function(x){unlist(strsplit(x, "_"))[2]}))
  links_pos_2 = as.numeric(sapply(svs_keep$id, function(x){unlist(strsplit(x, "_"))[5]}))

  links_labels = as.character(svs_keep$gene)

  #3. CNAs keep
  cnas_pat = filter(cnas, Patient == patient, ntot > 2 | ntot < 2, !(CHROM=="chrX"))
  cnas_pat = unique(cnas_pat[,c("CHROM", "Start", "End", "ntot", "Sample")])

  amps = filter(cnas_pat, ntot > 2)
  dels = filter(cnas_pat, ntot < 2)

  #AMPS
  grl_cna_list_amps = makeGRangesListFromDataFrame(amps,
    split.field ="Sample", seqnames.field = "CHROM",
    start.field = "Start", end.field = "End", keep.extra.columns=TRUE)
  #merged data
  merged_data = UnifyPeaks(grl_cna_list_amps)
  merged_results = as.data.table(findOverlaps(grl_cna_list_amps, merged_data))
  indices_keep = as.data.table(table(merged_results$subjectHits)) %>% filter(N == num_check)
  merged_results_amps = as.data.table(merged_data)[as.numeric(indices_keep$V1),]
  merged_results_amps$color = "firebrick2"

  #DELS
  grl_cna_list_dels = makeGRangesListFromDataFrame(dels,
    split.field ="Sample", seqnames.field = "CHROM",
    start.field = "Start", end.field = "End", keep.extra.columns=TRUE)
  #merged data
  merged_data = UnifyPeaks(grl_cna_list_dels)
  merged_results = as.data.table(findOverlaps(grl_cna_list_dels, merged_data))
  indices_keep = as.data.table(table(merged_results$subjectHits)) %>% filter(N == num_check)
  merged_results_dels = as.data.table(merged_data)[as.numeric(indices_keep$V1),]
  merged_results_dels$color = "dodgerblue2"

  cnas_pat = rbind(merged_results_amps, merged_results_dels)

  # Arcs coordinates
  cnas_pat$seqnames = as.character(cnas_pat$seqnames)
  arcs_chromosomes = sapply(cnas_pat$seqnames, function(x){unlist(strsplit(x, "chr"))[2]})
  arcs_begin = as.numeric(cnas_pat$start)
  arcs_end = as.numeric(cnas_pat$end)

  arcs_cols = cnas_pat$color

  #Make tracks for circos plot

  #first SNPs
  points_values = 0:4

  tracklist = BioCircosSNPTrack('mySNPTrack', as.numeric(points_chromosomes)-1, points_coordinates,
    points_values, colors = c("black"), minRadius = 0.8, maxRadius = 0.95, size=)

  # Background are always placed below other tracks
  tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack",
    minRadius = 0.8, maxRadius = 0.95,
    borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF")

  #CNAs
  tracklist = tracklist + BioCircosArcTrack('myArcTrack', as.numeric(arcs_chromosomes)-1, arcs_begin, arcs_end,
      minRadius = 0.68, maxRadius = 0.78, colors=arcs_cols)

  #SVs
  tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', as.numeric(links_chromosomes_1)-1, links_pos_1,
      links_pos_1 + 40000000, as.numeric(links_chromosomes_2)-1, links_pos_2, links_pos_2 + 550000,
      maxRadius = 0.65, color = "orange", labels = links_labels, labelSize=0.2,
      labelColor="blue", displayLabel=F)

  tracklist = tracklist + BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.65,
          borderSize = 0, fillColors = "#EEFFEE")

  circosplot = BioCircos(tracklist, genomeFillColor = "YlGnBu",
    displayGenomeBorder = TRUE,
    chrPad = 0.03, genome=c(248956422, 242193529, 198295559,
    190214555, 181538259, 170805979, 159345973, 145138636,
  138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
50818468),
       genomeTicksLen = 3, genomeTicksTextSize = 0,
       genomeTicksScale = 20000000,
       genomeLabelTextSize = 18, genomeLabelDy = 0, genomeBorderSize=0.7)

  print(circosplot)
  return(circosplot)
  print("done")
}

generate_circos("LY_RAP_0001")
generate_circos("LY_RAP_0002")
generate_circos("LY_RAP_0003")
