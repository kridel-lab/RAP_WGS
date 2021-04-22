#config-file.R

#----------------------------------------------------------------------
#karin isaev
#----------------------------------------------------------------------

date = Sys.Date()
print(date)
options(scipen=999) #no scientific notation

#----------------------------------------------------------------------
#load functions and libraries
#----------------------------------------------------------------------

options(stringsAsFactors=F)

#load libraries
packages <- c("dplyr", "ggplot2", "tidyr", "data.table", "plyr",
	"stringr")
lapply(packages, require, character.only = TRUE)
library(readxl)
library(GenomicRanges)
library(openxlsx)

#----------------------------------------------------------------------
#load mutation data
#----------------------------------------------------------------------

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text")

#DLBCL driver genes from Reddy et al 2017
reddy = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/Reddyetal_2017_driver_mutations.xlsx"))
reddy$type = "DLBCL"
treeomics_driver_3 = reddy[,1]
colnames(treeomics_driver_3)[1] = "Gene_Symbol"
write.csv(treeomics_driver_3, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/p003_drivers.csv", quote=F, row.names=F)

#PMBCL genes
pmbcl = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/PMBCL_genes.xlsx"))
pmbcl$type = "PMBCL"
pmbcl = unique(pmbcl)
treeomics_driver_2 = pmbcl[,1]
colnames(treeomics_driver_2)[1] = "Gene_Symbol"
write.csv(treeomics_driver_2, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/p002_drivers.csv", quote=F, row.names=F)

#MCL genes
mcl = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/MCL_genes.xlsx"))
mcl$type = "MCL"
mcl = unique(mcl)
treeomics_driver_1 = mcl[,1]
colnames(treeomics_driver_1)[1] = "Gene_Symbol"
write.csv(treeomics_driver_1, "/cluster/projects/kridelgroup/RAP_ANALYSIS/data/p001_drivers.csv", quote=F, row.names=F)

all_drivers = rbind(reddy, pmbcl, mcl)

#Copy number data
cnas = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_CNAs_by_Sequenza.rds")

#save driver genes for Treeomics
#column_name="Gene_Symbol"

#DLBCL mutations from Morin Blood 2013
morin = as.data.table(read_excel("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/supp_blood-2013-02-483727_TableS3.xlsx"))
genes_sum=as.data.table(table(morin$Gene))
genes_sum = as.data.table(filter(genes_sum, N > 5))
colnames(genes_sum)=c("Gene", "num_samples_w_mut")

#Our mutations
read_only_snvs = readRDS(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.rds")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS.rds"))])
read_only_indels = readRDS(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS_INDELS.rds")[length(list.files(pattern="READ_ONLY_ALL_MERGED_MUTS_INDELS.rds"))])
read_only_snvs$mut_type = "SNV"
read_only_indels$mut_type = "INDEL"

read_only = rbind(read_only_snvs, read_only_indels)

#sample info
samps = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/copy_RAP_masterlist_samples.rds")
colnames(samps)[4] ="Indiv"
z = which(samps$Indiv %in% read_only$Indiv)
samps = samps[z,]
samps[18,2] = "Kidney, NOS 2"

read_only$Tissue_Site = NULL
read_only$STUDY_PATIENT_ID = NULL
read_only$Specimen_Type = NULL

read_only = merge(read_only, samps, by = "Indiv")
#save file and download locally
#file_name=paste("/cluster/projects/kridelgroup/RAP_ANALYSIS/data/", date, "_", "Mutect2_Strelka_merged_mutations_wCNA_status.rds", sep="")
#saveRDS(read_only, file_name)

#mut-gene summary table for downstream use
mut_gene = unique(read_only[,c("mut_id", "symbol", "biotype", "Func.ensGene", "ExonicFunc.ensGene",
"AAChange.ensGene", "cosmic68", "mut_type")])

#1. how many patient samples is each mutation found in?
samples_per_mut = as.data.table(table(read_only$mut_id, read_only$STUDY_PATIENT_ID))

colnames(samples_per_mut) = c("mut_id", "patient", "num_of_samples_with_mut")
samples_per_mut = merge(samples_per_mut, mut_gene, by="mut_id")
samples_per_mut = samples_per_mut[order(-num_of_samples_with_mut)]
samples_per_mut = as.data.table(filter(samples_per_mut, num_of_samples_with_mut >0))

samples_per_mut$phylogeny = ""
z = which(samples_per_mut$patient == "LY_RAP_0001")
LY_RAP_0001 = samples_per_mut[z,]
LY_RAP_0001$phylogeny[LY_RAP_0001$num_of_samples_with_mut == 3] = "ancestor"
LY_RAP_0001$phylogeny[LY_RAP_0001$num_of_samples_with_mut == 1] = "private"
LY_RAP_0001$phylogeny[LY_RAP_0001$phylogeny == ""] = "shared"

z = which(samples_per_mut$patient == "LY_RAP_0002")
LY_RAP_0002 = samples_per_mut[z,]
LY_RAP_0002$phylogeny[LY_RAP_0002$num_of_samples_with_mut == 4] = "ancestor"
LY_RAP_0002$phylogeny[LY_RAP_0002$num_of_samples_with_mut == 1] = "private"
LY_RAP_0002$phylogeny[LY_RAP_0002$phylogeny == ""] = "shared"

z = which(samples_per_mut$patient == "LY_RAP_0003")
LY_RAP_0003 = samples_per_mut[z,]
LY_RAP_0003$phylogeny[LY_RAP_0003$num_of_samples_with_mut == 20] = "ancestor"
LY_RAP_0003$phylogeny[LY_RAP_0003$num_of_samples_with_mut == 1] = "private"
LY_RAP_0003$phylogeny[LY_RAP_0003$phylogeny == ""] = "shared"

samples_per_mut = rbind(LY_RAP_0001, LY_RAP_0002, LY_RAP_0003)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

get_known_signatures <- function(muttype = c("snv", "dbs", "indel", "tsb_snv"),
                                 source = c("COSMIC", "SIGNAL", "SPARSE"),
                                 sig_type = c("reference", "exposure", "tissue"),
                                 incl_poss_artifacts = FALSE,
                                 tissue_type = c(
                                   NA, "Biliary", "Bladder", "Bone",
                                   "Breast", "Cervix", "CNS",
                                   "Colorectal", "Esophagus", "Head",
                                   "Kidney", "Liver", "Lung",
                                   "Lymphoid", "Myeloid", "Ovary",
                                   "Pancreas", "Prostate", "Skin",
                                   "Stomach", "Thyroid", "Uterus"
                                 )) {

  # Validate arguments
  muttype <- match.arg(muttype)
  source <- match.arg(source)
  sig_type <- match.arg(sig_type)
  tissue_type <- match.arg(tissue_type)

  if (!.is_na(tissue_type) & sig_type != "tissue") {
    stop("tissue_type can only be used with `sig_type == 'tissue'`",
      call. = FALSE
    )
  }

  # Determine signature file name
  basename_sig <- paste0(muttype, "_", source, "_", sig_type, ".txt")
  fname_sig <- file.path("extdata", "signatures", basename_sig)
  fname_sig <- system.file(fname_sig, package = "MutationalPatterns")

  # Give error if file doesn't exist.
  if (!file.exists(fname_sig)) {
    stop(paste0(
      "The signature file: ", fname_sig, " does not exist.\n",
      "Look at the documentation of 'get_known_signatures()' for",
      " all the possible combinations of arguments."
    ),
    call. = FALSE
    )
  }

  # Read in signature file
  signatures <- read.table(fname_sig, sep = "\t", header = TRUE)


  # Remove meta columns
  if (muttype == "snv") {
    meta_cols <- c(1, 2)
  } else if (muttype == "tsb_snv") {
    meta_cols <- c(1, 2, 3)
  } else {
    meta_cols <- 1
  }
  signatures <- as.matrix(signatures[, -meta_cols, drop = FALSE])

  # Remove possible artifacts
  if (!incl_poss_artifacts) {
    if (source == "SIGNAL" & sig_type == "reference") {
      good_cols <- grep("Ref.Sig.N[0-9]{0-2}",
        colnames(signatures),
        invert = TRUE
      )
      signatures <- signatures[, good_cols, drop = FALSE]
    }

    if (source == "COSMIC" & muttype == "snv") {
      bad_sigs <- paste0("SBS", c(27, 43, seq(45, 60)))
      good_cols <- !colnames(signatures) %in% bad_sigs
      signatures <- signatures[, good_cols, drop = FALSE]
    }
  }

  # Select signatures of the specified tissue type
  if (!.is_na(tissue_type)) {
    tissue_cols <- grep(paste0("^", tissue_type, "_"), colnames(signatures))
    signatures <- signatures[, tissue_cols, drop = FALSE]
  }

  return(signatures)
}
