#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J process_manta
#SBATCH -c 8
	
module load R/3.5.0
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/

#define list of samples

for sample in $(cat patient_ids.txt) 
do
	echo $sample
	Rscript /cluster/projects/kridelgroup/RAP_ANALYSIS/scripts/Strelka_006_processing_manta_results.R $sample 

done
