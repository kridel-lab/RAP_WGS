#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ssm_parse
	
module load python2

#define list of CNV parsed files
tumor_sample_1=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_1_cluster3.segs.txt.parsed.txt_plusone.txt 
tumor_sample_2=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_2_cluster4.segs.txt.parsed.txt_plusone.txt 
tumor_sample_3=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_3_cluster3.segs.txt.parsed.txt_plusone.txt 
tumor_sample_4=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_4_cluster2.segs.txt.parsed.txt_plusone.txt 
tumor_sample_5=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_5_cluster4.segs.txt.parsed.txt_plusone.txt 
tumor_sample_6=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_6_cluster2.segs.txt.parsed.txt_plusone.txt 
tumor_sample_7=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_7_cluster4.segs.txt.parsed.txt_plusone.txt 
tumor_sample_8=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_8_cluster2.segs.txt.parsed.txt_plusone.txt 
tumor_sample_9=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_9_cluster2.segs.txt.parsed.txt_plusone.txt 
tumor_sample_10=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_10_cluster1.segs.txt.parsed.txt_plusone.txt 
tumor_sample_11=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_11_cluster2.segs.txt.parsed.txt_plusone.txt 
tumor_sample_12=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_12_cluster4.segs.txt.parsed.txt_plusone.txt 
tumor_sample_13=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_13_cluster1.segs.txt.parsed.txt_plusone.txt 
tumor_sample_14=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_14_cluster4.segs.txt.parsed.txt_plusone.txt 
tumor_sample_15=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_15_cluster1.segs.txt.parsed.txt_plusone.txt 
tumor_sample_16=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_16_cluster1.segs.txt.parsed.txt_plusone.txt 
tumor_sample_17=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_17_cluster2.segs.txt.parsed.txt_plusone.txt 
tumor_sample_18=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_18_cluster1.segs.txt.parsed.txt_plusone.txt 
tumor_sample_19=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_19_cluster1.segs.txt.parsed.txt_plusone.txt 
tumor_sample_20=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_20_cluster1.segs.txt.parsed.txt_plusone.txt 

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

awk 'BEGIN {srand()} !/^$/ { if (rand() <= .05 || FNR==1) print > "ssm_data_sample_10000_.txt"}' ssm_data.txt




