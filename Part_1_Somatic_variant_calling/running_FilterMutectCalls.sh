#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-19 # job array index

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar

#pwd
#/cluster/projects/kridelgroup/RAP_ANALYSIS/chr

#this is part 3 of the best protocols GATK steps for identifying SNVs and INDELS
#https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146

#final vcf files results of mutect2 for each patient 
#find -L . -name "*mutect_patient_results_all.vcf" > filter_jobs #20

#final final of contamination results for each patient 
#find -L . -name "*contamination.table" > contamination_tables_final #20

#re-order those two files so patient ids match 
#library(data.table)
#filters = fread("filter_jobs", header =F)
#conts = fread("contamination_tables_final", header =F)
#filters$pats = sapply(filters$V1, function(x){paste(unlist(strsplit(x, "_"))[2:6], collapse="_")})
#conts$pats = sapply(conts$V1, function(x){paste(unlist(strsplit(x, "_"))[1:5], collapse="_")})
#conts$pats = sapply(conts$pats, function(x){unlist(strsplit(x, "/"))[2]})

#filters = filters[order(pats)]
#conts = conts[order(pats)]

#filters$pats = NULL
#conts$pats = NULL
#write.table(filters, file="filter_jobs_ordered", quote=F, row.names=F, col.names=F)
#write.table(conts, file="contamination_tables_final_ordered", quote=F, row.names=F, col.names=F)

#pwd
names=($(cat filter_jobs_ordered))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

#contamination tables 
conts=($(cat contamination_tables_final_ordered))
echo ${conts[${SLURM_ARRAY_TASK_ID}]}

gatk FilterMutectCalls \
   -R /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
   -V ${names[${SLURM_ARRAY_TASK_ID}]} \
   --contamination-table ${conts[${SLURM_ARRAY_TASK_ID}]} \
   -O ${names[${SLURM_ARRAY_TASK_ID}]}_filtered.vcf.gz




