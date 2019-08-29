#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J make_ccfs
#SBATCH -c 8
	
module load python2
module load bedtools

#the way the files are listed below is the order in which the samplse are represented in the SSM input file 

#in R make bedtool input 
module load R/3.5.0

setwd("/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution")
library(data.table)
library(stringr)

files = list.files(pattern = ".segs.txt")
z = which(str_detect(files, "parsed")) ; z1 = which(str_detect(files, "bedtools"))
files=files[-c(z, z1)]

for(i in 1:length(files)){
	f = fread(files[i])
	#make bed input format 
	#chr	start	end 	name 	score(1-1000)	strand 
	f = f[,c(2:4)]
	write.table(f, file=paste(files[i], "bedtools_input.bed",sep=""), quote=F, row.names=F, col.names=F, sep="\t")
	print("done")
}

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/vcfs_annovar_annotated

#define list of CNV parsed files
tumor_sample_1=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_1_cluster3.segs.txtbedtools_input.bed
tumor_sample_2=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_2_cluster4.segs.txtbedtools_input.bed 
tumor_sample_3=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_3_cluster3.segs.txtbedtools_input.bed 
tumor_sample_4=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_4_cluster2.segs.txtbedtools_input.bed 
tumor_sample_5=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_5_cluster4.segs.txtbedtools_input.bed 
tumor_sample_6=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_6_cluster2.segs.txtbedtools_input.bed 
tumor_sample_7=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_7_cluster4.segs.txtbedtools_input.bed 
tumor_sample_8=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_8_cluster2.segs.txtbedtools_input.bed 
tumor_sample_9=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_9_cluster2.segs.txtbedtools_input.bed 
tumor_sample_10=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_10_cluster1.segs.txtbedtools_input.bed 
tumor_sample_11=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_11_cluster2.segs.txtbedtools_input.bed 
tumor_sample_12=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_12_cluster4.segs.txtbedtools_input.bed 
tumor_sample_13=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_13_cluster1.segs.txtbedtools_input.bed 
tumor_sample_14=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_14_cluster4.segs.txtbedtools_input.bed
tumor_sample_15=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_15_cluster1.segs.txtbedtools_input.bed 
tumor_sample_16=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_16_cluster1.segs.txtbedtools_input.bed 
tumor_sample_17=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_17_cluster2.segs.txtbedtools_input.bed 
tumor_sample_18=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_18_cluster1.segs.txtbedtools_input.bed 
tumor_sample_19=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_19_cluster1.segs.txtbedtools_input.bed 
tumor_sample_20=/cluster/projects/kridelgroup/RAP_ANALYSIS/results/titan/hmm/optimalClusterSolution/tumor_sample_20_cluster1.segs.txtbedtools_input.bed


cd /cluster/projects/kridelgroup/RAP_ANALYSIS/chr/vcfs_final/vcfs_annovar_annotated


library(data.table)
library(stringr)

files = list.files(pattern = "_multianno.txt")

for(i in 1:length(files)){
	f = fread(files[i])
	#make bed input format 
	#chr	start	end 	name 	score(1-1000)	strand 
	#f = f[,c(2:4)]
	write.table(f, file=paste(files[i], "bedtools_input.bed",sep=""), quote=F, row.names=F, col.names=F, sep="\t")
	print("done")
}


#define list of VCFs in same order to match CNV files
tumor_sample_1_vcf=LY_RAP_0003_Aut_FzT_01_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_2_vcf=LY_RAP_0003_Aut_FzT_02_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_3_vcf=LY_RAP_0003_Aut_FzT_03_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_4_vcf=LY_RAP_0003_Aut_FzT_04_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_5_vcf=LY_RAP_0003_Aut_FzT_05_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_6_vcf=LY_RAP_0003_Aut_FzT_06_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_7_vcf=LY_RAP_0003_Aut_FzT_07_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_8_vcf=LY_RAP_0003_Aut_FzT_11_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_9_vcf=LY_RAP_0003_Aut_FzT_09_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_10_vcf=LY_RAP_0003_Aut_FzT_10_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_11_vcf=LY_RAP_0003_Aut_FzT_12_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_12_vcf=LY_RAP_0003_Aut_FzT_13_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_13_vcf=LY_RAP_0003_Aut_FzT_14_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_14_vcf=LY_RAP_0003_Aut_FzT_15_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_15_vcf=LY_RAP_0003_Aut_FzT_16_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_16_vcf=LY_RAP_0003_Aut_FzT_17_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_17_vcf=LY_RAP_0003_Aut_FzT_18_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_18_vcf=LY_RAP_0003_Dia_FoT_01_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_19_vcf=LY_RAP_0003_Dia_FoT_03_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed
tumor_sample_20_vcf=LY_RAP_0003_Dia_FoT_05_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf_phylowgs_parser_input_no_indels.vcf.gz_annovar.vcf.gz.hg19_multianno.txtbedtools_input.bed

bedtools intersect -a $tumor_sample_1 -b $tumor_sample_1_vcf -wa -wb > {$tumor_sample_1_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_2 -b $tumor_sample_2_vcf -wa -wb > {$tumor_sample_2_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_3 -b $tumor_sample_3_vcf -wa -wb > {$tumor_sample_3_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_4 -b $tumor_sample_4_vcf -wa -wb > {$tumor_sample_4_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_5 -b $tumor_sample_5_vcf -wa -wb > {$tumor_sample_5_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_6 -b $tumor_sample_6_vcf -wa -wb > {$tumor_sample_6_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_7 -b $tumor_sample_7_vcf -wa -wb > {$tumor_sample_7_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_8 -b $tumor_sample_8_vcf -wa -wb > {$tumor_sample_8_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_9 -b $tumor_sample_9_vcf -wa -wb > {$tumor_sample_9_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_10 -b $tumor_sample_10_vcf -wa -wb > {$tumor_sample_10_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_11 -b $tumor_sample_11_vcf -wa -wb > {$tumor_sample_11_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_12 -b $tumor_sample_12_vcf -wa -wb > {$tumor_sample_12_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_13 -b $tumor_sample_13_vcf -wa -wb > {$tumor_sample_13_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_14 -b $tumor_sample_14_vcf -wa -wb > {$tumor_sample_14_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_15 -b $tumor_sample_15_vcf -wa -wb > {$tumor_sample_15_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_16 -b $tumor_sample_16_vcf -wa -wb > {$tumor_sample_16_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_17 -b $tumor_sample_17_vcf -wa -wb > {$tumor_sample_17_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_18 -b $tumor_sample_18_vcf -wa -wb > {$tumor_sample_18_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_19 -b $tumor_sample_19_vcf -wa -wb > {$tumor_sample_19_vcf}bedtools_output.bed
bedtools intersect -a $tumor_sample_20 -b $tumor_sample_20_vcf -wa -wb > {$tumor_sample_20_vcf}bedtools_output.bed










