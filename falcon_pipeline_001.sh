#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J falcon_001

#falcon_pipeline_001.sh
#2019-10-10 

module load samtools
module load vcftools

#obtain list of VCFs that have been filtered 
#mutect2 pipeline -> additional filters through selectvariants (no annovar)

find -L . -name "*phylowgs_parser_input_no_indels.vcf.gz" > names_somatic_muts_vcf_files.txt
bcftools merge -l names_somatic_muts_vcf_files.txt -Oz -o multisample.vcf.gz -m none


