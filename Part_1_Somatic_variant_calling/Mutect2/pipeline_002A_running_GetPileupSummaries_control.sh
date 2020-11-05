#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=41440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J controls
#SBATCH --array=0-71 # job array index 72 total files 24 chrs * 3 samples

#need to run getpileupsummaries on all bam files from normal samples
#author: Karin Isaev
#date started: June 25, 2019
#date updated: November 5th, 2020

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar

#pwd
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/MUTECT2_WORKDIR
#ls */*.bam > all_normal_chrs_bams.txt

#pwd
names=($(cat all_normal_chrs_bams.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
patient_name=${tum%_*_*_*}
chr=${MYVAR%.*}
chr="${test##*_}"

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle
output=/cluster/projects/burst2/MUTECT2_pileup

gatk GetPileupSummaries \
-R $fasta \
-I "${names[${SLURM_ARRAY_TASK_ID}]}" \
-L "${chr}" \
-O ${output}/${patient_name}_${chr}.pileups.table \
-V /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz
