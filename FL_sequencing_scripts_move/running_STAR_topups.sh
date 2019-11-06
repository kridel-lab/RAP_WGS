#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_fastq
#SBATCH --array=0-54 # job array index - number of jobs = numb of unique samples with top up runs 

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar
module load bam-readcount
module load gcc
module load star 

#Aligning reads using STAR is a two step process:

#Create a genome index
#Map reads to the genome

#pwd
cd /cluster/projects/kridelgroup/TGL13_transfer/data/fastq/topups_new_fastq_files

#---------------------------------------------------------------------

#pwd
names=($(cat topups_unique_sample_IDs.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

STAR --genomeDir /n/groups/hbctraining/intro_rnaseq_hpc/reference_data_ensembl38/ensembl38_STAR_index/ \
--runThreadN 6 \
--readFilesIn Mov10_oe_1.subset.fq \
--outFileNamePrefix ../results/STAR/Mov10_oe_1_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard 




