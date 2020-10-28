#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J strelka
#SBATCH -c 8
#SBATCH --array=0-26 # job array index because 27 total tumour samples

module load strelka/2.9.10
module load python
module load manta/1.6.0

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#MANTA
STRELKA_INSTALL_PATH=/cluster/tools/software/centos7/strelka/2.9.10

names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}
control_file=$(ls -d ${patient_name}_Ctl*)
str="LY_RAP_0003"

if [ "$patient_name" == "$str" ]; then
  control_file=$(ls $control_file/gatk/*.bam)
else
  control_file=$(ls $control_file/gatk/*.cram)
fi

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

mkdir /cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/STRELKA_WORKDIR_${tum_name}
STRELKA_ANALYSIS_PATH=/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR/STRELKA_WORKDIR_${tum_name}
MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_WORKDIR/MANTA_WORKDIR_${tum_name}

${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
--normalBam $control_file \
--tumorBam ${names[${SLURM_ARRAY_TASK_ID}]} \
--referenceFasta $fasta \
--indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz \
--runDir ${STRELKA_ANALYSIS_PATH}

#After succesfful configuration run the following:
${STRELKA_ANALYSIS_PATH}/runWorkflow.py -j 20 -m local
