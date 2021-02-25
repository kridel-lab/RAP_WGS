#-------------------------------------------------------------------------------
#pyclone-vi_001_prepare_subset_and_all_mutation_files.R
#Karin Isaev
#Thursday February 25, 20221
#load R/4.0.0
#-------------------------------------------------------------------------------

#load packages and data
source("/cluster/home/kisaev/RAP_WGS/config-file.R")

#prepare input data for pyclone and bamreadcount

#1. define patients
patients = unique(read_only$STUDY_PATIENT_ID)

#2. function to choose which mutations to include in the analysis

get_pyclone_input = function(patient){

  print(patient)

  #founding mutations don't provide much information for clonal reconstruction
  #but should include any founding mutations that fall ind driver genes
  #potentially functional ones

  if(patient == "LY_RAP_0001"){
    read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
    pat_founds = filter(as.data.table(table(read_only_pat$mut_id)), N == 3)

    pat_founds_keep = filter(read_only_pat, mut_id %in% pat_founds$V1,
      (ExonicFunc.ensGene %in% c("frameshift_deletion", "frameshift_insertion",
      "nonframeshift_deletion", "nonsynonymous_SNV", "stopgain",
      "stoploss") | (Func.ensGene %in% c("splicing", "exonic", "UTR3", "UTR5"))))$mut_id

    pat_founds_rm = unique(pat_founds[which(!(pat_founds$V1 %in% pat_founds_keep))]$V1)
  }

  if(patient == "LY_RAP_0002"){
    read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
    pat_founds = filter(as.data.table(table(read_only_pat$mut_id)), N == 4)

    pat_founds_keep = filter(read_only_pat, mut_id %in% pat_founds$V1,
      (ExonicFunc.ensGene %in% c("frameshift_deletion", "frameshift_insertion",
      "nonframeshift_deletion", "nonsynonymous_SNV", "stopgain",
      "stoploss") | (Func.ensGene %in% c("splicing", "exonic", "UTR3", "UTR5"))))$mut_id

    pat_founds_rm = unique(pat_founds[which(!(pat_founds$V1 %in% pat_founds_keep))]$V1)
  }

  if(patient == "LY_RAP_0003"){
    read_only_pat = as.data.table(filter(read_only, STUDY_PATIENT_ID == patient))
    pat_founds = filter(as.data.table(table(read_only_pat$mut_id)), N == 20)

    pat_founds_keep = filter(read_only_pat, mut_id %in% pat_founds$V1,
      (ExonicFunc.ensGene %in% c("frameshift_deletion", "frameshift_insertion",
      "nonframeshift_deletion", "nonsynonymous_SNV", "stopgain",
      "stoploss") | (Func.ensGene %in% c("splicing", "exonic", "UTR3", "UTR5"))))$mut_id

    pat_founds_rm = unique(pat_founds[which(!(pat_founds$V1 %in% pat_founds_keep))]$V1)
  }

  print(length(pat_founds_rm))
  print(dim(read_only_pat))

  print(length(unique(read_only_pat$mut_id)))

  #filter out unique mutations those only in one sample for clonal evolution analysis
  unique = filter(as.data.table(table(read_only_pat$mut_id)), N == 1)
  unique_muts_keep = filter(read_only_pat, mut_id %in% unique$V1,
    ExonicFunc.ensGene %in% c("frameshift_deletion", "frameshift_insertion",
    "nonframeshift_deletion", "nonsynonymous_SNV", "stopgain", "stoploss"))$mut_id
  unique_rm = unique[which(!(unique$V1 %in% unique_muts_keep))]$V1

  #run one version of pyclone with all mutations except for unique ones
  #and those with major copy number greater than 0
  pyclone_full = as.data.table(filter(read_only_pat, !(mut_id %in% unique_rm),
  Nmaj > 0, ntot >=2, !(mut_id %in% pat_founds_rm)))

  print(table(as.data.table(table(pyclone_full$mut_id))$N))

  print(patient)
  print(length(unique(pyclone_full$mut_id)))

  #create small subset for pyclone input for testing purposes
  pyclone_input = as.data.table(filter(pyclone_full,
  !(Func.ensGene %in% c("ncRNA_intronic", "intergenic", "intronic"))))

  print(table(as.data.table(table(pyclone_input$mut_id))$N))

  print(patient)
  print(length(unique(pyclone_input$mut_id)))

  #for mutations that are not present in all samples need to generate an entry for them
  #ideally need to get count of reads mapping there but for now just gonna put in zeros

  all_pat_datas = list(pyclone_full, pyclone_input)
  names(all_pat_datas) = c("pyclone_full", "pyclone_small_subset")
  print("done")
  return(all_pat_datas)
}

#get list of mutation files for each patient
all_pyclone_input = llply(patients, get_pyclone_input)
names(all_pyclone_input) = patients

#3. function to prepare bam readcount input files (full dataset)
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

#4. function to prepare bam readcount input files (subset of dataset)
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
