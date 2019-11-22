#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p long
#SBATCH --mem=20000M
#SBATCH -t 9-00:00 # Runtime in D-HH:MM
#SBATCH -J run_phylowgs_faster
#SBATCH -c 8
        

#date stamp = November 14
module load python2
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS
index=$1
echo $index

#run phylowgs
phylowgs=/cluster/home/kisaev/phylowgs/multievolve.py 

head ssm_data_ordered_$index.txt
mkdir full_run_1000_muts_$index

python2 $phylowgs --num-chains 4 --ssms ssm_data_ordered_$index.txt --cnvs cnv_data.txt -O full_run_1000_muts_$index

