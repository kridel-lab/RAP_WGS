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

#this is part 2 of the best protocols GATK steps for identifying SNVs and INDELS
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146

#Calculates the fraction of reads coming from cross-sample contamination, given results 
#from GetPileupSummaries. The resulting contamination table is used with FilterMutectCalls.

#using chromosome specific bam files

# get file name
#find -L . -name "*final.table*" > final_cram_tables 
#remove control sample from list
#nano final_cram_tables 

names=($(cat final_cram_tables))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

gatk CalculateContamination \
  -I "${names[${SLURM_ARRAY_TASK_ID}]}" \
  -matched /cluster/projects/kridelgroup/RAP_ANALYSIS/chr/control/chr/control_pileups_final.table \
  -O ${names[${SLURM_ARRAY_TASK_ID}]}_contamination.table

