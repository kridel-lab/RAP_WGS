#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J manta
#SBATCH -c 8
#SBATCH --array=0-27 # job array index because 28 total tumour samples

module load strelka/2.9.10
module load python
module load manta/1.6.0

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#MANTA
MANTA_INSTALL_PATH=/cluster/tools/software/centos7/manta/1.6.0

names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}
control_file=$(ls -d ${patient_name}_Ctl*)
control_file=$(ls $control_file/gatk/*.cram)
fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

mkdir /cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_WORKDIR/MANTA_WORKDIR_${tum_name}
MANTA_ANALYSIS_PATH=/cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_WORKDIR/MANTA_WORKDIR_${tum_name}

${MANTA_INSTALL_PATH}/bin/configManta.py \
--normalBam $control_file \
--tumorBam ${names[${SLURM_ARRAY_TASK_ID}]} \
--referenceFasta $fasta \
--runDir ${MANTA_ANALYSIS_PATH}

#After succesfful configuration run the following:
${MANTA_ANALYSIS_PATH}/runWorkflow.py -j 20

#move final files into /cluster/projects/kridelgroup/RAP_ANALYSIS/MANTA_RESULTS
