#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J strelka
#SBATCH -c 8
#SBATCH --array=0-19 # job array index

module load strelka/2.9.10
module load python
module load manta/1.6.0

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#MANTA
STRELKA_INSTALL_PATH=/cluster/tools/software/centos7/strelka/2.9.10

names=($(cat tum_samples_input_STRELKA_MANTA.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}

tum_loc=${MYVAR%/*}
tum_name=${MYVAR##*/}

mkdir STRELKA_WORKDIR_$tum_name
STRELKA_ANALYSIS_PATH=/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR_$tum_name

MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_WORKDIR_$tum_loc/$tum_name

${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
--normalBam /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam \
--tumorBam ${names[${SLURM_ARRAY_TASK_ID}]} \
--referenceFasta /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
--indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz \
--runDir ${STRELKA_ANALYSIS_PATH}

#After succesfful configuration run the following:
/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR_$tum_name/runWorkflow.py -j 20 -m local
