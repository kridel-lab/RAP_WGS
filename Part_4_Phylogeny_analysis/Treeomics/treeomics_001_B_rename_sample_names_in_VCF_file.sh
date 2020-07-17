#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_vars
#SBATCH --array=0-19 # job array index

#this follow the script
#RAP_WGS/Part_1_Somatic_variant_calling/Merged_Variant_Callers/Merged_005_prepare_bed_files_to_fillter_VCFs_by_soft_filters.R

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load tabix

#pwd
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/treeomics/src/input/mutect2_strelka_vcf_pcg_only

#pwd
names=($(cat /cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_vcf_files))
vcf_file=${names[${SLURM_ARRAY_TASK_ID}]}
echo $vcf_file

new_names=($(cat /cluster/projects/kridelgroup/RAP_ANALYSIS/data/samples_names))
pat=${new_names[${SLURM_ARRAY_TASK_ID}]}
echo $pat

#bgzip
bgzip $vcf_file

#index file
bcftools index $vcf_file.gz
gatk IndexFeatureFile -F $vcf_file.gz

#rename sample name
gatk RenameSampleInVcf \
    -I $vcf_file.gz \
    -O ${pat}.vcf \
    --NEW_SAMPLE_NAME $pat

rm $vcf_file.gz
rm $vcf_file.gz.tbi
rm $vcf_file.gz.csi
