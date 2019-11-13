#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J mutect2
#SBATCH --array=0-19 # job array index

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar

#pwd
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/chr

#collect all vcf files generated by strelka for each sample
#find -L . -name "*_filtered.vcf.gz" > mutect2_default_output

#pwd
names=($(cat mutect2_default_output))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

#gatk IndexFeatureFile \
#     -F ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum=${MYVAR:2:22}

gatk SelectVariants \
   -R /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
   -V ${names[${SLURM_ARRAY_TASK_ID}]} \
   -O mutect2_filtered/${tum}_mutect2_selectvariants.vcf.gz \
   --exclude-filtered true \
   --select-type-to-exclude INDEL\
   --exclude-intervals Y \

#normalize file 
input_vcf=mutect2_filtered/${tum}_mutect2_selectvariants.vcf.gz

#normalize variants, send to standard out and remove duplicates.
ref=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta
vt normalize ${input_vcf} -r $ref -o ${input_vcf}.normalized.vcf.gz

#DONE 
