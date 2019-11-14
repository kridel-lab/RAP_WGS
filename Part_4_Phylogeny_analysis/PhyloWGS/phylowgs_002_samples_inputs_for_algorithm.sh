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

#parsed SSM file from VCF files that were filtered for PhyloWGS SSMs
#first remove variants that are in this file that we decided to remove 

#DO NOT RUN
#module load R/3.5.0
#R
#library(data.table)
#phylo_keep = fread("/cluster/projects/kridelgroup/RAP_ANALYSIS/merged_MUTECT2_STRELKA/merged_variants_vcfs/vcf_summary_text/2019-11-14_PHYLOWGS_INPUT_MUTS.txt")
#ssm = fread("ssm_data.txt")
#z = which(ssm$gene %in% phylo_keep$mut_id)
#ssm = ssm[-z,]
#write.table(ssm, file="ssm_data.txt", quote=F, row.names=F, sep="\t")

filename=ssm_data.txt

for i in {1..10}
do
 echo $i 
 #need to make sure mutation IDs are labelled 0 to n... 
 #run accessory R script that will do this 
 Rscript /cluster/projects/kridelgroup/RAP_ANALYSIS/scripts/phylowgs_002_accessory_script_make_SSM_order.R $i
 sbatch /cluster/projects/kridelgroup/RAP_ANALYSIS/scripts/run_phyloWGS_ssms_algorithm.sh $i

done








