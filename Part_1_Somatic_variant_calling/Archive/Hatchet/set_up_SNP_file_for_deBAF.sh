#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=60000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J get_snps

module load R/4.0.0

Rscript /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/Hatchet/set_up_SNP_file_for_deBAF.R
