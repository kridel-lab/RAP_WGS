#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=60000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J sequenza

module load python2
module load samtools
module load tabix

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#The sequenza-utils command provides various tools; here we highlight only the basic usage:

#Process a FASTA file to produce a GC Wiggle track file:
fasta="/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta"

#From BAM files
#Normal and tumor BAM files

normal="/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam"
tumor="/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Dia_FoT_01_files/gatk/LY_RAP_0003_Dia_FoT_01.sorted.dup.recal.cram.bam"
sample="LY_RAP_0003_Dia_FoT_01"
out_folder="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Sequenza"

#1.
sequenza-utils gc_wiggle --fasta $fasta -o ${out_folder}/hg19.gc50Base.wig.gz -w 50
