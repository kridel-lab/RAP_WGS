#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J sequenza
#SBATCH -c 8
#SBATCH --array=0 # job array index

#get list of input .seqz files
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Sequenza
#ls *.small.seqz.gz > sequenza_small_seqz_files.txt

names=($(cat sequenza_small_seqz_files.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
sample=${names[${SLURM_ARRAY_TASK_ID}]}

#Run Rscript to get output files for each sample
module load R/3.5.0
Rscript /cluster/home/kisaev/RAP_WGS/Part_2_Somatic_copy_number_calling/sequenza/003_R_script_processing_patient1_sample5_only.R $sample
