#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_phylowgs
#SBATCH -c 8
	
module load python2

#the way the files are listed below is the order in which the samplse are represented in the SSM input file 

#define list of VCFs in same order to match CNV files
tumor_sample_1_vcf=LY_RAP_0003_Aut_FzT_01_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_2_vcf=LY_RAP_0003_Aut_FzT_02_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_3_vcf=LY_RAP_0003_Aut_FzT_03_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_4_vcf=LY_RAP_0003_Aut_FzT_04_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_5_vcf=LY_RAP_0003_Aut_FzT_05_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_6_vcf=LY_RAP_0003_Aut_FzT_06_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_7_vcf=LY_RAP_0003_Aut_FzT_07_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_8_vcf=LY_RAP_0003_Aut_FzT_11_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_9_vcf=LY_RAP_0003_Aut_FzT_09_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_10_vcf=LY_RAP_0003_Aut_FzT_10_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_11_vcf=LY_RAP_0003_Aut_FzT_12_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_12_vcf=LY_RAP_0003_Aut_FzT_13_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_13_vcf=LY_RAP_0003_Aut_FzT_14_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_14_vcf=LY_RAP_0003_Aut_FzT_15_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_15_vcf=LY_RAP_0003_Aut_FzT_16_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_16_vcf=LY_RAP_0003_Aut_FzT_17_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_17_vcf=LY_RAP_0003_Aut_FzT_18_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_18_vcf=LY_RAP_0003_Dia_FoT_01_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_19_vcf=LY_RAP_0003_Dia_FoT_03_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz
tumor_sample_20_vcf=LY_RAP_0003_Dia_FoT_05_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz

#make empty CNV file
#run phylowgs
phylowgs=/cluster/home/kisaev/phylowgs/multievolve.py 
#less ssm_data.txt | head -100 > test_ssm_data_input.txt
#python2 $phylowgs --num-chains 4 --ssms ssm_data_additional_soft_filters.txt --cnvs test_cnv_data.txt 
#python2 $phylowgs --num-chains 4 --ssms test_data_ssm_cleaned.txt --cnvs cnv_data.txt --burnin-samples 1 --mcmc-samples 1

python2 $phylowgs --num-chains 4 --ssms ssm_data.txt --cnvs cnv_data.txt

python2 $phylowgs --num-chains 2 --ssms ssm_data.txt --cnvs cnv_data.txt --burnin-samples 1 --mcmc-samples 1 -O testing_all_data
