#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J falcon_002

#falcon_pipeline_002.sh
#2019-10-10 

module load samtools
module load vcftools
module load R/3.5.0

#use list of VCF files with haployupe calls 
#to run FALCON R scripts 

Rscript falcon_R_script_001.R

echo "done"


