#!/bin/bash
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -J citup
#SBATCH -c 8
#SBATCH -t 5-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

#running citup
module load citup
module load python3

run_citup_qip.py 2020-07-29_rap_citup_input_freqs.txt \
2020-07-29_rap_citup_input_clusters.txt p003_results.h5 --submit local

#less 2020-07-27_rap_citup_input_freqs.txt | head -100 > test_citup_input_freqs.txt
#less 2020-07-27_rap_citup_input_clusters.txt | head -100 > test_citup_input_clusters.txt

#run_citup_qip.py test_citup_input_freqs.txt \
t#est_citup_input_clusters.txt test.h5 --submit local
