#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J sciclone

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone

Rscript /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/SciClone/sciclone_002_run.R
