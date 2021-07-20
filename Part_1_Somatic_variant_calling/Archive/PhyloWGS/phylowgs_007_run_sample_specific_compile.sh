#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_phylowgs
#SBATCH -c 8

module load python2
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS

#define list of samples
accessory_script=/cluster/projects/kridelgroup/RAP_ANALYSIS/scripts/phylowgs_008_run_sample_specific_phyloWGS.sh

for sample in $(cat sample_list.txt) 
do
	sbatch $accessory_script $sample
done