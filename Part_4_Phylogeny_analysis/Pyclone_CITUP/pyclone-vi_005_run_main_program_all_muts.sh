#!/bin/bash
#SBATCH -p long
#SBATCH --mem=30720M
#SBATCH -J pyclone-vi-all
#SBATCH -c 10
#SBATCH -t 21-00:00 # Runtime in D-HH:MM
#SBATCH --array=0-2 # job array index

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone

source activate pyclone-vi
module load python3

#------------------------------------------------------------------
#Set up input data
#------------------------------------------------------------------

#ls *all_samples_pyclonevi_* > pyclone_input_files.tsv

names=($(cat pyclone_input_files.tsv))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}

patient_name=${MYVAR##*all_samples_pyclonevi_}
patient_name=${MYVAR%_pyclone_input.tsv*}

echo $patient_name

density=beta-binomial

#------------------------------------------------------------------
#Run pyclone-vi
#------------------------------------------------------------------

folder="$(date +"%d-%m-%Y")"
mkdir $folder

pyclone-vi fit -i $MYVAR \
-o $folder/${patient_name}_${density}_rap_wgs_all_muts.h5 -c 20 -d $density -r 50 --seed 123 --print-freq 2

pyclone-vi write-results-file -i $folder/${patient_name}_${density}_rap_wgs_all_muts.h5 -o $folder/${patient_name}_${density}_rap_wgs_all_muts.tsv
