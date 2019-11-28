#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ssm_parse
#SBATCH -c 8
	
module load python2
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS

#define list of samples

for sample in $(cat sample_list.txt) 
do
	echo $sample
	mkdir phyloWGS_$sample 

	#define VCF file for sample 
	vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/phyloWGS_VCF_files/${sample}_mutect2_phyloWGS_input_just_tumour.vcf
	echo $vcf

	#define CNA file for sample 
	cnv=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/CNAs_titan/_${sample}_cna_datat_input_phyloWGS_CNV_parser.txt.parsed.txt
	echo $cnv

	sbatch /cluster/projects/kridelgroup/RAP_ANALYSIS/scripts/phylowgs_005_accessory_script.sh $vcf $cnv $sample 
done






