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

#define list of VCFs in same order to match CNV files
tumor_sample_1_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Dia_FoT_05_merged_variants.vcf.vcf
tumor_sample_2_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Dia_FoT_03_merged_variants.vcf.vcf
tumor_sample_3_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Dia_FoT_01_merged_variants.vcf.vcf
tumor_sample_4_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_18_merged_variants.vcf.vcf
tumor_sample_5_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_17_merged_variants.vcf.vcf
tumor_sample_6_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_16_merged_variants.vcf.vcf
tumor_sample_7_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_15_merged_variants.vcf.vcf
tumor_sample_8_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_14_merged_variants.vcf.vcf
tumor_sample_9_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_13_merged_variants.vcf.vcf
tumor_sample_10_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_12_merged_variants.vcf.vcf
tumor_sample_11_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_11_merged_variants.vcf.vcf
tumor_sample_12_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_10_merged_variants.vcf.vcf
tumor_sample_13_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_09_merged_variants.vcf.vcf
tumor_sample_14_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_07_merged_variants.vcf.vcf
tumor_sample_15_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_06_merged_variants.vcf.vcf
tumor_sample_16_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_05_merged_variants.vcf.vcf
tumor_sample_17_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_04_merged_variants.vcf.vcf
tumor_sample_18_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_03_merged_variants.vcf.vcf
tumor_sample_19_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_02_merged_variants.vcf.vcf
tumor_sample_20_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/filtered_LY_RAP_0003_Aut_FzT_01_merged_variants.vcf.vcf

#python /cluster/home/kisaev/phylowgs/parser/create_phylowgs_inputs.py --cnvs tumor_sample_15=$tumor_sample_15 --cnvs tumor_sample_16=$tumor_sample_16 --vcf-type tumor_sample_15=mutect_smchet --vcf-type tumor_sample_16=mutect_smchet tumor_sample_15=$tumor_sample_15_vcf tumor_sample_16=$tumor_sample_16_vcf --verbose 
parser=/cluster/home/kisaev/phylowgs/parser/create_phylowgs_inputs.py 

python2 $parser --vcf-type s1=mutect_smchet --vcf-type s2=mutect_smchet --vcf-type s3=mutect_smchet --vcf-type s4=mutect_smchet \
--vcf-type s5=mutect_smchet --vcf-type s6=mutect_smchet --vcf-type s7=mutect_smchet --vcf-type s8=mutect_smchet \
--vcf-type s9=mutect_smchet --vcf-type s10=mutect_smchet --vcf-type s11=mutect_smchet --vcf-type s12=mutect_smchet \
--vcf-type s13=mutect_smchet --vcf-type s14=mutect_smchet --vcf-type s15=mutect_smchet --vcf-type s16=mutect_smchet \
--vcf-type s17=mutect_smchet --vcf-type s18=mutect_smchet --vcf-type s19=mutect_smchet --vcf-type s20=mutect_smchet \
s1=$tumor_sample_1_vcf s2=$tumor_sample_2_vcf s3=$tumor_sample_3_vcf s4=$tumor_sample_4_vcf \
s5=$tumor_sample_5_vcf s6=$tumor_sample_6_vcf s7=$tumor_sample_7_vcf s8=$tumor_sample_8_vcf \
s9=$tumor_sample_9_vcf s10=$tumor_sample_10_vcf s11=$tumor_sample_11_vcf s12=$tumor_sample_12_vcf \
s13=$tumor_sample_13_vcf s14=$tumor_sample_14_vcf s15=$tumor_sample_15_vcf s16=$tumor_sample_16_vcf \
s17=$tumor_sample_17_vcf s18=$tumor_sample_18_vcf s19=$tumor_sample_19_vcf s20=$tumor_sample_20_vcf --regions=all --verbose 

#awk 'BEGIN {srand()} !/^$/ { if (rand() <= .05 || FNR==1) print > "ssm_data_sample_10000_.txt"}' ssm_data.txt


#sample random 10000 (aprox 10% of file) variants for phylowgs to run 
filename=ssm_data.txt

#DO NOT RUN 

#awk 'BEGIN {srand()} !/^$/ { if (rand() <= .1 || FNR==1) print > "ssm_data_sample_10000.txt"}' $filename

#need to rename the mutation IDs in SSM file for phylowgs to run 
#module load R/3.5.0
#R

#library(data.table)
#f = fread("ssm_data_sample_10000.txt")
#f$id = paste("s", 0:(nrow(f)-1), sep="")
#write.table(f, file="ssm_data_sample_10000_new_ids.txt", sep="\t", quote=F, row.names=F)




