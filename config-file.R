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

signatures = readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/known_signatures_MutationalPatterns.rds")

fit_to_signatures_strict <- function(mut_matrix, signatures, max_delta = 0.004) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  rowname <- . <- NULL

  #Set colnames if absent, to prevent duplicate names later.
  if (is.null(colnames(mut_matrix))){
    colnames(mut_matrix) <- seq_len(ncol(mut_matrix))
  }

  # Remove signatures with zero contribution across samples
  fit_res <- fit_to_signatures(mut_matrix, signatures)
  sig_pres <- rowSums(fit_res$contribution) != 0
  my_signatures_total <- signatures[, sig_pres, drop = FALSE]
  nsigs <- ncol(my_signatures_total)

  # perform signature selection per sample
  all_results <- vector("list", ncol(mut_matrix))
  for (i in seq(1, ncol(mut_matrix))) {
    my_signatures <- my_signatures_total
    mut_mat_sample <- mut_matrix[, i, drop = FALSE]

    # Fit again
    fit_res <- fit_to_signatures(mut_mat_sample, my_signatures)
    sim <- .get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)

    # Keep track of the cosine similarity and which signatures are removed.
    sims <- vector("list", nsigs)
    sims[[1]] <- sim
    removed_sigs <- vector("list", nsigs)
    removed_sigs[[1]] <- "None"

    # Sequentially remove the signature with the lowest contribution
    for (j in seq(2, nsigs)) {

      # Remove signature with the weakest relative contribution
      contri_order <- fit_res$contribution %>%
        prop.table(2) %>%
        rowSums() %>%
        order()
      weakest_sig_index <- contri_order[1]
      weakest_sig <- colnames(my_signatures)[weakest_sig_index]
      removed_sigs[[j]] <- weakest_sig
      signatures_sel <- my_signatures[, -weakest_sig_index, drop = FALSE]


      # Fit with new signature selection
      fit_res <- fit_to_signatures(mut_mat_sample, signatures_sel)
      sim_new <- .get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)

      if (is.nan(sim_new) == TRUE) {
        sim_new <- 0
        warning("New similarity between the original and the reconstructed
                        spectra after the removal of a signature was NaN.
                        It has been converted into a 0.
                        This happened with the following fit_res:")
        print(fit_res)
      }
      sims[[j]] <- sim_new

      # Check if the loss in cosine similarity between the original vs reconstructed after removing the signature is below the cutoff.
      delta <- sim - sim_new
      if (delta <= max_delta) {
        my_signatures <- signatures_sel
        sim <- sim_new
      }
      else {
        break
      }
    }

    # Plot how the cosine similarities decayed
    sim_decay_fig <- .plot_sim_decay(sims, removed_sigs, max_delta)

    # Perform final fit on selected signatures
    fit_res <- fit_to_signatures(mut_mat_sample, my_signatures)

    # Add data of sample to list.
    results <- list("sim_decay_fig" = sim_decay_fig, "fit_res" = fit_res)
    all_results[[i]] <- results
  }

  # Get decay figs and fit_res in separate lists
  decay_figs <- purrr::map(all_results, "sim_decay_fig")
  fit_res <- purrr::map(all_results, "fit_res")

  # Combine the contribution of all samples
  contribution <- purrr::map(fit_res, "contribution") %>%
    purrr::map(function(x) tibble::rownames_to_column(as.data.frame(x))) %>%
    purrr::reduce(dplyr::full_join, by = "rowname")

  # Fix signature order of contribution and add absent sigs to
  # keep the legend colors consistent for plotting.
  sig_ref <- tibble::tibble("rowname" = colnames(signatures))
  contribution <- dplyr::left_join(sig_ref, contribution, by ="rowname") %>%
    as.data.frame()

  # Turn contribution into matrix and remove NAs
  rownames(contribution) <- contribution$rowname
  contribution <- contribution %>%
    dplyr::select(-rowname) %>%
    as.matrix()
  contribution[is.na(contribution)] <- 0

  # Combine the reconstructed of all samples
  reconstructed <- purrr::map(fit_res, "reconstructed") %>%
    do.call(cbind, .)

  # Combine all and return
  fit_res <- list("contribution" = contribution, "reconstructed" = reconstructed)
  results <- list("sim_decay_fig" = decay_figs, "fit_res" = fit_res)
  return(results)
}



#' Get the cosine similarity between a reconstructed mutation matrix and the original
#'
#' @param mut_matrix mutation count matrix (dimensions: x mutation types
#' X n samples)
#' @param fit_res Named list with signature contributions and reconstructed
#' mutation matrix
#'
#' @return Cosine similarity
#' @noRd
#'
.get_cos_sim_ori_vs_rec <- function(mut_matrix, fit_res) {
  cos_sim_all <- cos_sim_matrix(mut_matrix, fit_res$reconstructed)
  cos_sim <- diag(cos_sim_all)
  mean_cos_sim <- mean(cos_sim)
  return(mean_cos_sim)
}


#' Plot decay in cosine similarity as signatures are removed.
#'
#' This function is called by fit_to_signatures_strict
#'
#' @param sims List of cosine similarities
#' @param removed_sigs List of iteratively removed signatures
#' @param max_delta The maximum difference in original vs reconstructed cosine similarity.
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @noRd
#' @return ggplot object
#'
.plot_sim_decay <- function(sims, removed_sigs, max_delta) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  Removed_signatures <- Cosine_similarity <- NULL

  # Prepare data
  sims <- sims[!S4Vectors::isEmpty(sims)] %>%
    unlist()
  removed_sigs <- removed_sigs[!S4Vectors::isEmpty(removed_sigs)] %>%
    unlist()
  tb <- tibble::tibble(
    "Cosine_similarity" = sims,
    "Removed_signatures" = factor(removed_sigs, levels = removed_sigs)
  )

  # Determine if the final removed signature exceeded the cutoff.
  sims_l <- length(sims)
  col <- rep("low_delta", sims_l)
  final_delta <- sims[sims_l - 1] - sims[sims_l]
  if (final_delta > max_delta) {
    col[sims_l] <- "high_delta"
  }

  fig <- ggplot(data = tb, aes(x = Removed_signatures, y = Cosine_similarity, fill = col)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      limits = c("low_delta", "high_delta"),
      values = c("grey", "red"),
      guide = FALSE
    ) +
    labs(
      x = "Removed signatures",
      y = paste0("Cosine similarity (max delta: ", max_delta, ")")
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
      text = element_text(size = 12)
    )
  return(fig)
}

merge_signatures <- function(signatures, cos_sim_cutoff = 0.8, merge_char = ";", verbose = TRUE) {

  # These variables use non standard evaluation.
  # To avoid R CMD check complaints we initialize them to NULL.
  . <- NULL

  # Validate signature names
  nr_lowercase <- colnames(signatures) %>%
    grepl(merge_char, .) %>%
    sum()
  if (nr_lowercase > 0) {
    stop(paste0("Please remove all ", merge_char, " characters from your signature names."), call. = FALSE)
  }

  # Determine max similarity between signatures
  sim_m <- cos_sim_matrix(signatures, signatures)
  diag(sim_m) <- 0
  max <- max(sim_m)

  # Merge signatures while max similarity is higher than cutoff.
  while (max > cos_sim_cutoff) {

    # Find signatures that need to be merged
    max_index <- order(sim_m, decreasing = TRUE)[1]
    max_loc <- arrayInd(max_index, dim(sim_m), useNames = TRUE)
    sigs_left <- signatures[, -max_loc, drop = FALSE]
    sigs_to_combi <- signatures[, max_loc, drop = FALSE]

    # Signatures that have already been merged and thus exist
    # of multiple signatures are weighted accordingly.
    weights <- sigs_to_combi %>%
      colnames() %>%
      stringr::str_count(merge_char) %>%
      magrittr::add(1)
    combi_sig <- sigs_to_combi %*% diag(weights)

    # Merge signatures
    combi_sig <- combi_sig %>%
      rowSums() %>%
      matrix()
    combi_sig <- combi_sig / sum(weights)
    colnames(combi_sig) <- paste(colnames(sigs_to_combi), collapse = merge_char)

    # Add merged signature to the rest.
    signatures <- cbind(sigs_left, combi_sig)

    # Print which signatures have been merged
    if (verbose) {
      merged_sig_names <- paste0(colnames(sigs_to_combi), collapse = ", ")
      message(paste0("Combined the following two signatures: ", merged_sig_names))
    }

    # Determine max similarity between signatures for next loop
    sim_m <- cos_sim_matrix(signatures, signatures)
    diag(sim_m) <- 0
    max <- max(sim_m)
  }
  return(signatures)
}

#library(randomcoloR)
#n <- 50
#mypal <- distinctColorPalette(n)
#saveRDS(mypal, file="/cluster/projects/kridelgroup/RAP_ANALYSIS/mypal.rds")
mypal=readRDS("/cluster/projects/kridelgroup/RAP_ANALYSIS/mypal.rds")
