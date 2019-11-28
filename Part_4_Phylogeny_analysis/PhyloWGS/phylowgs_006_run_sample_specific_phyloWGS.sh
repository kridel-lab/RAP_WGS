#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ssm_parse
#SBATCH -c 8
	
module load python2
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS

#define list of samples

for sample in $(cat sample_list.txt) 
do
	echo $sample
	sbatch /cluster/projects/kridelgroup/RAP_ANALYSIS/scripts/phylowgs_006_run_sample_specific_phyloWGS_accessory_script.sh $sample 

done






