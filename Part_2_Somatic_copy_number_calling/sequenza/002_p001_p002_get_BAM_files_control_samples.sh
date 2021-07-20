#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_vars
#SBATCH --array=0-1 # job array index

#mutations were filtered to include only those to be analyzed in pyclone
#some of those mutations are not present across all patients
#read counts need to be obtained at those regions

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load tabix
module load bam-readcount

#pwd
cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#pwd
names=($(cat control_samples_001_002_only.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

#sample CRAM file
bam_file=${names[${SLURM_ARRAY_TASK_ID}]}

samtools view -T $fasta -bh ${names[${SLURM_ARRAY_TASK_ID}]} > ${tum_name}.bam
samtools index ${tum_name}.bam

#then BAM files were moved into /cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files
