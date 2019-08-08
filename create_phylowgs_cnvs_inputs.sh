#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J cnv_parse
	
module load python2

#define list of CNV parsed files
tumor_sample_15=tumor_sample_15_cluster1.segs.txt.parsed.txt
tumor_sample_16=tumor_sample_16_cluster1.segs.txt.parsed.txt

#define list of VCFs in same order to match CNV files
tumor_sample_15_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/LY_RAP_0003_Aut_FzT_16_mutect_patient_results_all.vcf_filtered.vcf.gz_norm.vcf
tumor_sample_16_vcf=/cluster/projects/kridelgroup/RAP_ANALYSIS/chr/LY_RAP_0003_Aut_FzT_17_mutect_patient_results_all.vcf_filtered.vcf.gz_annovar_new.vcf.gz.hg19_multianno.vcf

#python /cluster/home/kisaev/phylowgs/parser/create_phylowgs_inputs.py --cnvs tumor_sample_15=$tumor_sample_15 --cnvs tumor_sample_16=$tumor_sample_16 --vcf-type tumor_sample_15=mutect_smchet --vcf-type tumor_sample_16=mutect_smchet tumor_sample_15=$tumor_sample_15_vcf tumor_sample_16=$tumor_sample_16_vcf --verbose 

python /cluster/home/kisaev/phylowgs/parser/create_phylowgs_inputs.py --cnvs tumor_sample_15=$tumor_sample_15  --vcf-type tumor_sample_15=mutect_smchet tumor_sample_15=$tumor_sample_15_vcf --verbose