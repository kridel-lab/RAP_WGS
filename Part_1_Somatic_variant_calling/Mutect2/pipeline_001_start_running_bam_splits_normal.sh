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

for chrom in `seq 1 22` X Y

do 
	samtools view -T /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta -bh /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam ${chrom} > chr/control_sample_${chrom}.bam 
	samtools index chr/control_sample_${chrom}.bam 

done

