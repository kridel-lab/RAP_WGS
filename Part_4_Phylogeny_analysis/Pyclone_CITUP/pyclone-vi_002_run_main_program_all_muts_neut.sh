#!/bin/bash
#SBATCH -p long
#SBATCH --mem=30720M
#SBATCH -J pyclone-vi-neut
#SBATCH -c 10
#SBATCH -t 21-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

source activate pyclone-vi
module load python3

pyclone-vi fit -i all_samples_pyclonevi_all_muts_neutral_pyclone_input.tsv -o rap_wgs_all_muts_neut.h5 -c 80 -d binomial -r 100
pyclone-vi write-results-file -i rap_wgs_all_muts_neut.h5 -o rap_wgs_all_muts_neut.tsv
