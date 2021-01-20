#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J sequenza
#SBATCH -c 8
#SBATCH --array=0-26 # job array index because 27 total tumour samples

module load strelka/2.9.10
module load python
module load manta/1.6.0

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

names=($(cat input_sequenza_bam_files.txt))
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
