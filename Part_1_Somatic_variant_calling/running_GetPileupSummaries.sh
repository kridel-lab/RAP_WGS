#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-479 # job array index

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

#using chromosome specific bam files

# get file name
#find -L . -name "*.cram.bai" > mutect_jobs #480
#sed 's/.bai//' mutect_jobs > mutect_jobs_clean

names=($(cat mutect_jobs_clean))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

chr=($(echo ${names[${SLURM_ARRAY_TASK_ID}]} | awk -F'[_.]' '{print $8}'))
echo "${chr}"

gatk GetPileupSummaries \
-R /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
-I "${names[${SLURM_ARRAY_TASK_ID}]}" \
-L "${chr}" \
-O ${names[${SLURM_ARRAY_TASK_ID}]}.pileups.table \
-V /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz


