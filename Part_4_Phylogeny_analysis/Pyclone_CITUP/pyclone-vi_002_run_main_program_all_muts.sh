#!/bin/bash
#SBATCH -p long
#SBATCH --mem=30720M
#SBATCH -J pyclone-vi-all
#SBATCH -c 10
#SBATCH -t 21-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

source activate pyclone-vi
module load python3

pyclone-vi fit -i all_samples_pyclonevi_all_muts_LY_RAP_0001_pyclone_input.tsv \
-o LY_RAP_0001_rap_wgs_all_muts.h5 -c 80 -d binomial -r 100

pyclone-vi write-results-file -i LY_RAP_0001_rap_wgs_all_muts.h5 -o LY_RAP_0001_rap_wgs_all_muts.tsv

pyclone-vi fit -i all_samples_pyclonevi_all_muts_LY_RAP_0002_pyclone_input.tsv \
-o LY_RAP_0002_rap_wgs_all_muts.h5 -c 80 -d binomial -r 100

pyclone-vi write-results-file -i LY_RAP_0002_rap_wgs_all_muts.h5 -o LY_RAP_0002_rap_wgs_all_muts.tsv

#pyclone-vi fit -i all_samples_pyclonevi_all_muts_LY_RAP_0003_pyclone_input.tsv \
#-o LY_RAP_0003_rap_wgs_all_muts.h5 -c 80 -d binomial -r 100

#pyclone-vi write-results-file -i LY_RAP_0003_rap_wgs_all_muts.h5 -o LY_RAP_0003_rap_wgs_all_muts.tsv
