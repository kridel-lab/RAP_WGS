#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_vars
#SBATCH --array=0-26 # job array index

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
echo ${pat}

#use MUTECT2 VCF file since its format is more appropriate for Treeomics
vcf_file=/cluster/projects/kridelgroup/RAP_ANALYSIS/MUTECT2_selected_VCFs/${pat}.selected.normalized.vcf.gz

#variants specific to each patient found here:
variant_file_all_muts=ANALYSIS/Treeomics/input_dat/${pat}_treeomics_input_all_muts_.bed
variant_file_pcgs_only=ANALYSIS/Treeomics/input_dat/${pat}_treeomics_input_pcgs_only_.bed

#RUN ONLY ONCE
#first make tbi for vcf
#tabix -p vcf $vcf_file

#extract filtered variants from mutect2 VCF files so that Treeomics has
#less variants to analyze and hopefully

#first ALL MUTATIONS
#tabix -fhB $vcf_file $variant_file_all_muts > treeomics/src/input/mutect2_strelka_all_muts/${pat}_mutect2_treeomics_input.vcf

intersectBed -a ${vcf_file} -b ${variant_file_all_muts} -header > treeomics/src/input/mutect2_strelka_all_muts/${pat}_mutect2_treeomics_input.vcf
bcftools view -s $pat treeomics/src/input/mutect2_strelka_all_muts/${pat}_mutect2_treeomics_input.vcf > treeomics/src/input/mutect2_strelka_all_muts/${pat}_mutect2_treeomics_input_just_tumour.vcf
rm treeomics/src/input/mutect2_strelka_all_muts/${pat}_mutect2_treeomics_input.vcf

#first ONLY PCGs
tabix -fhB $vcf_file $variant_file_pcgs_only > treeomics/src/input/mutect2_strelka_pcgs_only/${pat}_mutect2_treeomics_input.vcf

intersectBed -a ${vcf_file} -b ${variant_file_pcgs_only} -header > treeomics/src/input/mutect2_strelka_pcgs_only/${pat}_mutect2_treeomics_input.vcf
bcftools view -s $pat treeomics/src/input/mutect2_strelka_pcgs_only/${pat}_mutect2_treeomics_input.vcf > treeomics/src/input/mutect2_strelka_pcgs_only/${pat}_mutect2_treeomics_input_just_tumour.vcf
rm  treeomics/src/input/mutect2_strelka_pcgs_only/${pat}_mutect2_treeomics_input.vcf
