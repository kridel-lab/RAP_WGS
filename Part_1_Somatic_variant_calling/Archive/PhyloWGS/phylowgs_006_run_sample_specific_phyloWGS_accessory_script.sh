#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_phylowgs
#SBATCH -c 8

module load python2

sample=$1
echo $sample

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/phyloWGS_$sample 

phylowgs=/cluster/home/kisaev/phylowgs/multievolve.py
rm -r phylowgs_output_$sample
mkdir phylowgs_output_$sample

python2 $phylowgs --num-chains 2 --ssms ssm_data.txt --cnvs cnv_data.txt -O phylowgs_output_$sample --burnin-samples 200 --mcmc-samples 500


