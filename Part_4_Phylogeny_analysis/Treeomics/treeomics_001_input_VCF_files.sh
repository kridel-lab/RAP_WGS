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
cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#pwd
names=($(cat /cluster/projects/kridelgroup/RAP_ANALYSIS/patient_ids.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
pat=${names[${SLURM_ARRAY_TASK_ID}]}

vcf_file=MUTECT2_RESULTS/mutect2_filtered/${pat}_mutect2_selectvariants.vcf.gz.normalized.vcf.gz
variant_file=merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text/2019-11-14_TREEOMICS_INPUT_MUTS_int_with_VCFs.bed

#first make tbi for vcf
#tabix -p vcf $vcf_file 

#extract filtered variants from mutect2 VCF files so that Treeomics has less variants to analyze and hopefully 
tabix -fhB $vcf_file $variant_file > treeomics/src/input/STRELKA_VCF_INPUTS/${pat}_mutect2_treeomics_input.vcf
bcftools view -s $pat treeomics/src/input/STRELKA_VCF_INPUTS/${pat}_mutect2_treeomics_input.vcf > treeomics/src/input/STRELKA_VCF_INPUTS/${pat}_mutect2_treeomics_input_just_tumour.vcf
rm treeomics/src/input/STRELKA_VCF_INPUTS/${pat}_mutect2_treeomics_input.vcf