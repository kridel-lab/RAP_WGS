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

#less 2019-09-10_rap_citup_input_freqs.txt | head -500 > test_citup_input_freqs.txt
#less 2019-09-10_rap_citup_input_clusters.txt | head -500 > test_citup_input_clusters.txt

#run_citup_qip.py test_citup_input_freqs.txt test_citup_input_clusters.txt test.h5 --submit local 