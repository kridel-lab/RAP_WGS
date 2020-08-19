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

#VarScan expects its input in SAMtools pileup
#format, which is obtained from a BAM file via the samtools pileup command.
#For example:

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta
normal_sample=/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam
outdir=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/varscan

samtools pileup -f $fasta $index > ${outdir}/${sample}.pileup
java -jar VarScan.jar pileup2snp ${outdir}/${sample}.pileup
