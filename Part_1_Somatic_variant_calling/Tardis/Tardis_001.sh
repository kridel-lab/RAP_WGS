#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J manta
#SBATCH -c 8
#SBATCH --array=0-19 # job array index
	
module load strelka/2.9.10  
module load python
module load manta/1.6.0 
module load tardis 

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

names=($(cat tum_samples_input_STRELKA_MANTA.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}

tum_loc=${MYVAR%/*}
tum_name=${MYVAR##*/}

tardis -i ${names[${SLURM_ARRAY_TASK_ID}]} \
--ref /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
--sonic /cluster/home/kisaev/sonic-prebuilt/human_g1k_v37.sonic \
--out $tum_name






