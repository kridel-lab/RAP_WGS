#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-47 # job array index

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

#using output from samtools (split original bam files into chrosome based files)

# get file name
#find -L . -name "*.cram.bai" > mutect_jobs
#remove the control bam file from the list using nano

#get file name just for patients 02 and 18
#find -L . -name "*LY_RAP_0003_Aut_FzT_02*" > mutect021802
#find -L . -name "*LY_RAP_0003_Aut_FzT_18*" > mutect021818
#cat mutect021802 mutect021818 > mutect_jobs_new
#sed '/bai/d' ./mutect_jobs_new > mutect0218 

#pwd
names=($(cat mutect0218))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

chr=($(echo ${names[${SLURM_ARRAY_TASK_ID}]} | awk -F'[_.]' '{print $8}'))
echo "${chr}"

gatk Mutect2 \
-R /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
-I "${names[${SLURM_ARRAY_TASK_ID}]}" \
-I /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam \
-L "${chr}" \
-tumor "${tum}" \
-O ${names[${SLURM_ARRAY_TASK_ID}]}.vcf.gz \
--germline-resource /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz

