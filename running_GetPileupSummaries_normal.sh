#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2

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

#this is part 2 of the best protocols GATK steps for identifying SNVs and INDELS
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146

#Summarizes counts of reads that support reference, alternate and 
#other alleles for given sites. Results can be used with CalculateContamination.

#get pileup summary table for matched normal sample to use with 

gatk GetPileupSummaries \
-R /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
-I /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam \
-O match_normal.pileups.table \
-V /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz \
-L /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz \


