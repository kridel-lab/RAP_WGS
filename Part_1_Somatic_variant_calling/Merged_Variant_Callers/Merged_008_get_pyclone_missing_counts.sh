#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=31440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J get_missing_counts

script=/cluster/home/kisaev/RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_008_get_pyclone_missing_counts.R
Rscript $script LY_RAP_0001
Rscript $script LY_RAP_0002
Rscript $script LY_RAP_0003
