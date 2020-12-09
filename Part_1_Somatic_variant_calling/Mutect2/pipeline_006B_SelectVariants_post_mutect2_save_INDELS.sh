#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=21440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-26 # job array index because 27 total tumour samples

#needs to run on bam files
module load strelka/2.9.10
module load python
module load manta/1.6.0
module load CNVkit
module load samtools
module load vt

cd /cluster/projects/burst2
mkdir MUTECT2_selected_VCFs

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
filtered=/cluster/projects/burst2/MUTECT2_final_VCFs
selected=/cluster/projects/burst2/MUTECT2_selected_VCFs

#make index for VCF file
gatk IndexFeatureFile \
   --input ${filtered}/${tum_name}.filtered.vcf.gz

gatk SelectVariants \
   -R $fasta \
   -V ${filtered}/${tum_name}.filtered.vcf.gz \
   -O ${selected}/${tum_name}.indels.selected.vcf.gz \
   --exclude-filtered true \
   --select-type-to-exclude SNP \
   --exclude-intervals Y \

#normalize file
input_vcf=${selected}/${tum_name}.indels.selected.vcf.gz

#normalize variants, send to standard out and remove duplicates.
vt normalize ${input_vcf} -r $fasta -o ${selected}/${tum_name}.indels.selected.normalized.vcf.gz

#DONE
