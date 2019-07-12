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

#this is the last part of the best protocols GATK steps for identifying SNVs and INDELS
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146

#find -L . -name "*vcf_filtered.vcf.gz" > annovar_jobs #20

#pwd
names=($(cat annovar_jobs))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

bcftools norm -m -any ${names[${SLURM_ARRAY_TASK_ID}]} -f /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta > ${names[${SLURM_ARRAY_TASK_ID}]}_norm.vcf
table_annovar.pl --buildver hg19 ${names[${SLURM_ARRAY_TASK_ID}]}_norm.vcf /cluster/tools/software/annovar/humandb --protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f --nastring . --outfile ${names[${SLURM_ARRAY_TASK_ID}]}_annovar_new.vcf.gz --vcfinput
