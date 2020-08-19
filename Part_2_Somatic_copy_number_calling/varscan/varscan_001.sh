#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J varscan

module load samtools
module load python3
module load R/3.6.1
module load varscan

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/varscan

index=$1

#VarScan expects its input in SAMtools pileup
#format, which is obtained from a BAM file via the samtools pileup command.
#For example:

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta

samtools pileup -f $fasta $index > ${index}.pileup
java -jar VarScan.jar pileup2snp ${index}.pileup
