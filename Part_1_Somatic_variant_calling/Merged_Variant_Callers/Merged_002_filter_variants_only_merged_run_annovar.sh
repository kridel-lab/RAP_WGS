#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --mem=30000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_vars
#SBATCH --array=0-26 # job array index

#how did we get here?
#mutect2 was run with basic filters
#strelka was run with basic filters
#identified mutations that were called in both --> regions file
#go back to VCF files produeced by MUTECT2 (since they have fields required for Annvoar to work)
#first extract only those common variants from those vcf files (the ones that were called by both tools)
#via tabix
#then using this new reduced VCF. file --> run annovar to annotate variants
#then in next script
#filter these variants using soft filters like population frequencies and genes

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load tabix

#pwd
cd /cluster/projects/kridelgroup/RAP_ANALYSIS

names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}

vcf_file=/cluster/projects/burst2/MUTECT2_selected_VCFs/${tum_name}.selected.normalized.vcf.gz
variant_file=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/_${tum_name}_merged_mutations.bed

tabix -fhB $vcf_file $variant_file > merged_MUTECT2_STRELKA/merged_variants_vcfs/${tum_name}_merged_variants.vcf

#RUN ANNOVAR
anno_input=merged_MUTECT2_STRELKA/merged_variants_vcfs/${tum_name}_merged_variants.vcf

table_annovar.pl --buildver hg19 ${anno_input} /cluster/tools/software/annovar/humandb --protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f --outfile /cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcfs_annovar_annotated/${tum_name}_annovar.vcf.gz --vcfinput
