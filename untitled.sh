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

#this is part 1 of the best protocols GATK steps for identifying SNVs and INDELS
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146

#using output from samtools (split original bam files into chrosome based files)

# get file name
#find -L . -name "*.cram.bai" > mutect_jobs
#remove the control bam file from the list using nano

#total VCF files at the end = 20 * 24 (chromosomes) = 480
