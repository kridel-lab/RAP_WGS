#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-19 # job array index

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar

#pwd
#/cluster/projects/kridelgroup/RAP_ANALYSIS/chr
names=($(cat all_patient_ids.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

gatk GatherVcfsCloud \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_1.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_2.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_3.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_4.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_5.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_6.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_7.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_8.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_9.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_10.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_11.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_12.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_13.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_14.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_15.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_16.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_17.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_18.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_19.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_20.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_21.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_22.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_X.cram.vcf.gz \
-I ${names[${SLURM_ARRAY_TASK_ID}]}_Y.cram.vcf.gz \
-O ${names[${SLURM_ARRAY_TASK_ID}]}_mutect_patient_results_all.vcf





