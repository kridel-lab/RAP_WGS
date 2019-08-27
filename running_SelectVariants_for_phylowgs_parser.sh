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
module load python3
module load gatk
module load annovar

#pwd
#/cluster/projects/kridelgroup/RAP_ANALYSIS/chr

#final final of contamination results for each patient 
#find -L . -name "*vcf.gz_norm.vcf" > normalized_vcf_files_pre_annovar
#ls *vcf.gz_norm.vcf > normalized_vcf_files_pre_annovar

#for annovar annotated files try
ls *annovar_new.vcf.gz.hg19_multianno.vcf > normalized_vcf_files_after_annovar

#pwd
names=($(cat normalized_vcf_files_pre_annovar))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

tum_sample=$(echo ${names[${SLURM_ARRAY_TASK_ID}]} | cut -d'_' -f 1,2,3,4,5,6)
echo ${tum_sample}

gatk IndexFeatureFile \
     -F ${names[${SLURM_ARRAY_TASK_ID}]}

gatk SelectVariants \
   -R /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
   -V ${names[${SLURM_ARRAY_TASK_ID}]} \
   -O vcfs_final/${names[${SLURM_ARRAY_TASK_ID}]}_final_selected_variants.vcf.gz \
   --exclude-filtered true \
   --select-type-to-exclude INDEL\
   -select "DP > 100 && POP_AF < 0.001" \
   -sn ${tum_sample} \
   --exclude-intervals Y \
   --exclude-intervals X \



