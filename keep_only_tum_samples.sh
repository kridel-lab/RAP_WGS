#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J STRELKA
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
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs

for filename in *.vcf 
	do bcftools view -S /cluster/projects/kridelgroup/RAP_ANALYSIS/patient_ids.txt -o filtered_$filename.vcf $filename --force-samples
done


