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

#MANTA 
MANTA_INSTALL_PATH=/cluster/tools/software/centos7/manta/1.6.0

names=($(cat tum_samples_input_STRELKA_MANTA.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

mkdir MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}
MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}

${MANTA_INSTALL_PATH}/bin/configManta.py \
--normalBam /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam \
--tumorBam ${names[${SLURM_ARRAY_TASK_ID}]} \
--referenceFasta /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
--runDir ${MANTA_ANALYSIS_PATH}

#After succesfful configuration run the following:
/cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_WORKDIR_${names[${SLURM_ARRAY_TASK_ID}]}/runWorkflow.py -j 20


