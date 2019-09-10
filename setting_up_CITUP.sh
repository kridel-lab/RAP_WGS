#!/bin/bash
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -J citup
#SBATCH -c 8
#SBATCH -t 5-00:00 # Runtime in D-HH:MM

pwd
#/cluster/projects/kridelgroup/RAP_ANALYSIS/pyclone/rap_pyclone/rap_pyclone

#running citup
module load citup
module load python3

run_citup_qip.py 2019-09-10_rap_citup_input_freqs.txt 2019-09-10_rap_citup_input_clusters.txt p001_results.h5 --submit local 
