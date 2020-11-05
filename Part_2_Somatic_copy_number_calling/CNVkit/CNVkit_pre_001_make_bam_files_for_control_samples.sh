#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J norm_samps
#SBATCH -c 8

#needs to run on bam files

module load strelka/2.9.10
module load python
module load manta/1.6.0
module load CNVkit
module load samtools

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle
output=/cluster/projects/burst2/CNVkit_WORKDIR/Normal_Samples

tum_name=LY_RAP_0001
file_name=/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0001_Ctl_FzG_01_files/gatk/LY_RAP_0001_Ctl_FzG_01.sorted.dup.recal.cram

samtools view -b  -T $fasta -o ${output}/${tum_name}.bam ${file_name}

tum_name=LY_RAP_0002
file_name=/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0002_Ctl_FzG_07_files/gatk/LY_RAP_0002_Ctl_FzG_07.sorted.dup.recal.cram

samtools view -b  -T $fasta -o ${output}/${tum_name}.bam ${file_name}

tum_name=LY_RAP_003
scp /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam ${output}/${tum_name}.bam
