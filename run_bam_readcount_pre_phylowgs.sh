#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J bamreadcount
#SBATCH --array=0-19 # job array index

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount

#pwd
#/cluster/projects/kridelgroup/RAP_ANALYSIS/

#this is the last part of the best protocols GATK steps for identifying SNVs and INDELS
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146

#find -L . -name "*cram" > jobs #20

#pwd
names=($(cat jobs))

echo ${names[${SLURM_ARRAY_TASK_ID}]}

bam-readcount -f /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta ${names[${SLURM_ARRAY_TASK_ID}]} -l /cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/mutations_to_retrieve_from_bam_files.bed

