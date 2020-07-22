#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_all_scripts

#1. merge variants from strelka and mutect2


#----Pyclone

#prepare mutations for bamreadcount
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_prepare_bed_files_to_fillter_VCFs_by_soft_filters.R

#run bamreadcount
sbatch /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Pyclone_get_missing_mutations.sh

#assemble results from bamreadcount into pyclone mutation list
Rscript /cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_008_get_pyclone_missing_counts.R

#prepare individual sample mutations
Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/Pyclone_CITUP/pyclone_001_make_input_files.R
