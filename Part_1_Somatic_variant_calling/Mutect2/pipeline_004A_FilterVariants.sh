#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=41440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-647 # job array index

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019
#date updated: November 4th, 2020

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar

#pwd
cd /cluster/projects/burst2
#ls */*.bam > all_chrs_bams.txt

#pwd
names=($(cat all_chrs_bams.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

chr=($(echo ${names[${SLURM_ARRAY_TASK_ID}]} | awk -F'[_.]' '{print $8}'))
echo "${chr}"

patient_name=${tum%_*_*_*}
normal_file=/cluster/projects/kridelgroup/RAP_ANALYSIS/MUTECT2_WORKDIR/CHR_split/${patient_name}_${chr}.bam

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle
output=/cluster/projects/burst2/MUTECT2_contamination

gatk FilterMutectCalls \
    -V /cluster/projects/burst2/MUTECT2_raw_VCFs/${tum}_${chr}.vcf.gz \
    --contamination-table ${output}/${tum}_${chr}_contamination.table \
    -O /cluster/projects/burst2/MUTECT2_filtered_VCFs/${tum}_${chr}.vcf.gz
