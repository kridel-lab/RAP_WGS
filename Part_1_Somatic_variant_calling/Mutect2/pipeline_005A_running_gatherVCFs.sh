#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-26 # job array index because 27 total tumour samples

#needs to run on bam files

module load strelka/2.9.10
module load python
module load manta/1.6.0
module load CNVkit
module load samtools

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}
control_file=/cluster/projects/burst2/CNVkit_WORKDIR/Normal_Samples/${patient_name}.bam

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle
mkdir MUTECT2_final_VCFs
output=/cluster/projects/burst2/MUTECT2_final_VCFs
filtered=/cluster/projects/burst2/MUTECT2_filtered_VCFs

gatk GatherVcfsCloud \
-I ${filtered}/${tum_name}_1.vcf.gz \
-I ${filtered}/${tum_name}_2.vcf.gz \
-I ${filtered}/${tum_name}_3.vcf.gz \
-I ${filtered}/${tum_name}_4.vcf.gz \
-I ${filtered}/${tum_name}_5.vcf.gz \
-I ${filtered}/${tum_name}_6.vcf.gz \
-I ${filtered}/${tum_name}_7.vcf.gz \
-I ${filtered}/${tum_name}_8.vcf.gz \
-I ${filtered}/${tum_name}_9.vcf.gz \
-I ${filtered}/${tum_name}_10.vcf.gz \
-I ${filtered}/${tum_name}_11.vcf.gz \
-I ${filtered}/${tum_name}_12.vcf.gz \
-I ${filtered}/${tum_name}_13.vcf.gz \
-I ${filtered}/${tum_name}_14.vcf.gz \
-I ${filtered}/${tum_name}_15.vcf.gz \
-I ${filtered}/${tum_name}_16.vcf.gz \
-I ${filtered}/${tum_name}_17.vcf.gz \
-I ${filtered}/${tum_name}_18.vcf.gz \
-I ${filtered}/${tum_name}_19.vcf.gz \
-I ${filtered}/${tum_name}_20.vcf.gz \
-I ${filtered}/${tum_name}_21.vcf.gz \
-I ${filtered}/${tum_name}_22.vcf.gz \
-I ${filtered}/${tum_name}_X.vcf.gz \
-I ${filtered}/${tum_name}_Y.vcf.gz \
-O ${output}/${tum_name}.filtered.vcf.gz
