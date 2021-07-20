#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ssm_parse
#SBATCH -c 8

module load python2

vcf=$1 #this refers to the VCF file 
echo $vcf

cnv=$2 #this refers o the CNA file entered as input in previous script 
echo $cnv

sample=$3
echo $sample

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/phyloWGS_$sample 

parser=/cluster/home/kisaev/phylowgs/parser/create_phylowgs_inputs.py 

python2 $parser --cnvs sample1=$cnv --vcf-type sample1=mutect_smchet sample1=$vcf

