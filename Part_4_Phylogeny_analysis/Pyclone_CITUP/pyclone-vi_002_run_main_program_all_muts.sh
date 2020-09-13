#!/bin/bash
#SBATCH -p long
#SBATCH --mem=30720M
#SBATCH -J pyclone-vi-all
#SBATCH -c 10
#SBATCH -t 21-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

source activate pyclone-vi
module load python3

pyclone-vi fit -i all_samples_pyclonevi_all_muts_pyclone_input.tsv -o rap_wgs_all_muts.h5 -c 40 -d binomial -r 50
pyclone-vi write-results-file -i rap_wgs_all_muts.h5 -o rap_wgs_all_muts.tsv
