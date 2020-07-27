#!/bin/bash
#SBATCH -p veryhimem
#SBATCH --mem=184320M
#SBATCH -J citup
#SBATCH -c 8
#SBATCH -t 5-00:00 # Runtime in D-HH:MM

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

#running citup
module load citup
module load python3

run_citup_qip.py 2020-07-27_rap_citup_input_freqs.txt \
2020-07-27_rap_citup_input_clusters.txt p003_results.h5 --submit local

#less 2019-12-04_rap_citup_input_freqs.txt | head -500 > test_citup_input_freqs.txt
#less 2019-12-04_rap_citup_input_clusters.txt | head -500 > test_citup_input_clusters.txt

#run_citup_qip.py test_citup_input_freqs.txt test_citup_input_clusters.txt test.h5 --submit local
