#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

#purpose: complete final soft filters
#mainly intersect with copy number file
#generate filters required for Treeomics, PhyloWGS and Pyclone
#so that VCFs can be filtered accordingly

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

date = Sys.Date()
options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "readr", "ggplot2", "vcfR", "tidyr", "mclust", "data.table",
              "plyr", "openxlsx", "readxl", "GenomicRanges",
              "ggrepel", "stringr", "maftools", "magrittr", "ggExtra", "broom")
lapply(packages, require, character.only = TRUE)

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#----------------------------------------------------------------------
#purpose
#----------------------------------------------------------------------

#prepare input mutation files for pyclone, treeomics and phylowgs

#----------------------------------------------------------------------
#data
#----------------------------------------------------------------------

#1. Summary SNV data
muts = readRDS(list.files(pattern="mut_calls.rds")[length(list.files(pattern="mut_calls.rds"))]) #get most recent mutation file

#blacklist mutations were already removed
#TLOD filter and minimum # of alternative alleles was also filtered
table(muts$Indiv)
unique(muts$Indiv)

#2. Sample summary
#dna = fread("RAP_DNA.txt") ; dna=dna[,1:3] ; colnames(dna)[2] = "barcode"; dna$barcode = as.numeric(dna$barcode)
#biops = fread("RAP_FFbiopsies_extracted.txt" ); biops = biops[,1:6] ; colnames(biops)[4] = "barcode"
#dna = merge(dna, biops, by="barcode")
#colnames(dna)[2] = "Indiv"
#colnames(dna)[7] = "Tissue_Site"
#colnames(dna)[8] = "Specimen_Type"
#dna$Specimen_Type = "FT"

#dna = as.data.table(filter(dna, Indiv %in% muts$Indiv))
#ffpe = as.data.table(matrix(ncol=ncol(dna), nrow=3))
#colnames(ffpe) = colnames(dna)
#ffpe = as.data.frame(ffpe)

#ffpe$Indiv = c("LY_RAP_0003_Dia_FoT_05", "LY_RAP_0003_Dia_FoT_01" ,"LY_RAP_0003_Dia_FoT_03")
#ffpe$barcode =c("15:S12966E", "15:S12966A", "15:S12966C")
#ffpe$Tissue_Site = c("left_breast", "right_neck_LN", "left_axilla_LN")
#ffpe$Specimen_Type = "FFPE"
#ffpe$DNA = "DNA"
#ffpe$STUDY_PATIENT_ID = "LY_RAP_0003"
#dna = rbind(dna, ffpe)
#dna$id = paste(dna$Specimen_Type, dna$Tissue_Site, dna$barcode, sep="_")

samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")
colnames(samps)[4] ="Indiv"
z = which(samps$Indiv %in% muts$Indiv)
samps = samps[z,]

muts = merge(muts, samps, by="Indiv", all=TRUE)
unique(muts$Indiv)

#Mutation data from WGS Morin et al 2013
#Table S3. All somatic SNVs identified from 40 genome pairs and 13 cell lines (XLSX, 918 KB)
#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))
treeomics_drivers = reddy
colnames(treeomics_drivers)[1] = "Gene_Symbol"
write.csv(treeomics_drivers, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/data/reddy_drivers.csv",
quote=F, row.names=F)

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_TITAN.rds")

#----------------------------------------------------------------------
#analysis
#----------------------------------------------------------------------

#overlap SNVs with CNAs

overlap_snvs_cnas = function(sample){
  print(sample)
  snvs = as.data.table(filter(muts, Indiv==sample))
  cnas_pat = as.data.table(filter(cnas, Sample == sample))

  #make granges objects
  snvs$CHROM = paste("chr", snvs$CHROM, sep="")
  snvs_gr = GRanges(
    seqnames = snvs$CHROM,
    ranges = IRanges(snvs$POS, end = snvs$POS),
    strand = rep("*", length(snvs$POS)),
    score = 1:length(snvs$POS))

  cnas_gr = GRanges(
    seqnames = cnas_pat$CHROM,
    ranges = IRanges(cnas_pat$Start, end = cnas_pat$End),
    strand = rep("*", length(cnas_pat$Start)),
    logR = cnas_pat$LogR)

  #intersect them
  #Then subset the original objects with the negative indices of the overlaps:
  hits <- findOverlaps(snvs_gr, cnas_gr, ignore.strand=TRUE)
  hits_overlap = cbind(snvs[queryHits(hits),], cnas_pat[subjectHits(hits),])
  print(head(hits_overlap))
  return(hits_overlap)
}

muts_wCNAs = as.data.table(ldply(llply(unique(muts$Indiv), overlap_snvs_cnas, .progress = "text")))
muts_wCNAs = muts_wCNAs[order(CHROM, POS)]
table(muts_wCNAs$Copy_Number)

#----------------------------------------------------------------------
#SEPERATE FILES INTO FINAL SETS BASED ON TOOL INPUT
#----------------------------------------------------------------------

#1. READ-ONLY FILE = FINAL MUTATIONS = KEEP THIS WAY UNLESS MAJOR CHANGE NEEDED
read_only = as.data.table(muts_wCNAs)
read_only$tot_counts = read_only$Ref_counts+read_only$alt_counts

#2. PYCLONE = REMOVE UNIQUE MUTATIONS
#REMOVE MUTATION WITH NMAJ OF 0
#KEEP WES GENE MUTATIONS TO SIMPLIFY

read_only$isdriver=""
read_only$isdriver[which((read_only$symbol %in% reddy$Gene) | (read_only$symbol %in% genes_sum$Gene))] = "yes"

saveRDS(read_only, file=paste(date, "READ_ONLY_ALL_MERGED_MUTS.rds", sep="_"))
patients = unique(read_only$STUDY_PATIENT_ID)

#for clonal evolution analysis only keep founder mutations that are in DLBCL genes

get_pyclone_input = function(patient){

  print(patient)

    if(patient == "LY_RAP_0001"){
    read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
    pat_founds_keep = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 3, V2=="yes")
    pat_founds_remove = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 3, V2=="")
  }

  if(patient == "LY_RAP_0002"){
    read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
    pat_founds_keep = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 4, V2=="yes")
    pat_founds_remove = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 4, V2=="")
  }

  if(patient == "LY_RAP_0003"){
    read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
    pat_founds_keep = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 20, V2=="yes")
    pat_founds_remove = filter(as.data.table(table(read_only_pat$mut_id, read_only_pat$isdriver)), N == 20, V2=="")
  }

  print(dim(read_only_pat))
  #filter out unique mutations those only in one sample for clonal evolution analysis
  unique = filter(as.data.table(table(read_only_pat$mut_id)), N == 1)

  #run one version of pyclone with all mutations except for unique ones
  #and those with major copy number greater than 0
  pyclone_full = as.data.table(filter(read_only_pat, !(mut_id %in% unique$V1),
  !(mut_id %in% pat_founds_remove$V1),
  MajorCN > 0, Copy_Number >=2))

  print(patient)
  print(length(unique(pyclone_full$mut_id)))

  #create small subset for pyclone input for testing purposes
  pyclone_input = as.data.table(filter(pyclone_full,
  !(Func.ensGene %in% c("ncRNA_intronic", "intergenic", "intronic"))))

  print(patient)
  print(length(unique(pyclone_input$mut_id)))

  #for mutations that are not present in all samples need to generate an entry for them
  #ideally need to get count of reads mapping there but for now just gonna put in zeros

  all_pat_datas = list(pyclone_full, pyclone_input)
  names(all_pat_datas) = c("pyclone_full", "pyclone_small_subset")
  print("done")
  return(all_pat_datas)
}

all_pyclone_input = llply(patients, get_pyclone_input)
names(all_pyclone_input) = patients

#get input for bamreadcount
get_bam_read_full_dat = function(dat){

  name_analysis = names(dat)[1]
  patient = unique(dat[[1]]$STUDY_PATIENT_ID)
  print(name_analysis)

  if(patient == "LY_RAP_0001"){
    t = filter(as.data.table(table(dat[[1]]$mut_id)), (N==3))
    muts_all = as.data.table(filter(dat[[1]], mut_id %in% t$V1))
    muts_some = as.data.table(filter(dat[[1]], !(mut_id %in% t$V1)))
  }

  if(patient == "LY_RAP_0002"){
    t = filter(as.data.table(table(dat[[1]]$mut_id)), (N==4))
    muts_all = as.data.table(filter(dat[[1]], mut_id %in% t$V1))
    muts_some = as.data.table(filter(dat[[1]], !(mut_id %in% t$V1)))
  }

  if(patient == "LY_RAP_0003"){
    t = filter(as.data.table(table(dat[[1]]$mut_id)), (N==20))
    muts_all = as.data.table(filter(dat[[1]], mut_id %in% t$V1))
    muts_some = as.data.table(filter(dat[[1]], !(mut_id %in% t$V1)))
  }

  print(patient)
  print(dim(muts_all))
  print(dim(muts_some))

  #save muts_some so that can run bam readcount and extract counts in those mutations
  #across all samples
  #need chr, start, end, ref and alt save as bed file, no colnames, tab sep
  muts_some_bam_readcount = unique(muts_some[,c("CHROM", "POS", "mut_id", "REF", "ALT")])
  muts_some_bam_readcount$mut_id = muts_some_bam_readcount$POS
  #remove chr from CHROM
  muts_some_bam_readcount$CHROM = sapply(muts_some_bam_readcount$CHR, function(x){
    unlist(strsplit(x, "chr"))[2]
  })
  print(dim(muts_some_bam_readcount))
  file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/", name_analysis, "_", patient, "_pyclone_bam_readcount_input.bed", sep="")
  write.table(muts_some_bam_readcount,
    file,
    col.names=F, quote=F, row.names=F, sep="\t")
  print("done")
}

get_bam_read_partial_dat = function(dat){

  name_analysis = names(dat)[2]
  patient = unique(dat[[1]]$STUDY_PATIENT_ID)
  print(name_analysis)

  if(patient == "LY_RAP_0001"){
    t = filter(as.data.table(table(dat[[2]]$mut_id)), (N==3))
    muts_all = as.data.table(filter(dat[[2]], mut_id %in% t$V1))
    muts_some = as.data.table(filter(dat[[2]], !(mut_id %in% t$V1)))
  }

  if(patient == "LY_RAP_0002"){
    t = filter(as.data.table(table(dat[[2]]$mut_id)), (N==4))
    muts_all = as.data.table(filter(dat[[2]], mut_id %in% t$V1))
    muts_some = as.data.table(filter(dat[[2]], !(mut_id %in% t$V1)))
  }

  if(patient == "LY_RAP_0003"){
    t = filter(as.data.table(table(dat[[2]]$mut_id)), (N==20))
    muts_all = as.data.table(filter(dat[[2]], mut_id %in% t$V1))
    muts_some = as.data.table(filter(dat[[2]], !(mut_id %in% t$V1)))
  }

  print(patient)
  print(dim(muts_all))
  print(dim(muts_some))

  #save muts_some so that can run bam readcount and extract counts in those mutations
  #across all samples
  #need chr, start, end, ref and alt save as bed file, no colnames, tab sep
  muts_some_bam_readcount = unique(muts_some[,c("CHROM", "POS", "mut_id", "REF", "ALT")])
  muts_some_bam_readcount$mut_id = muts_some_bam_readcount$POS
  #remove chr from CHROM
  muts_some_bam_readcount$CHROM = sapply(muts_some_bam_readcount$CHR, function(x){
    unlist(strsplit(x, "chr"))[2]
  })
  print(dim(muts_some_bam_readcount))
  file=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/", name_analysis, "_", patient, "_pyclone_bam_readcount_input.bed", sep="")
  write.table(muts_some_bam_readcount,
    file,
    col.names=F, quote=F, row.names=F, sep="\t")
  print("done")
}

llply(all_pyclone_input, get_bam_read_full_dat, .progress="text")
llply(all_pyclone_input, get_bam_read_partial_dat, .progress="text")
