#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J run_phylowgs
#SBATCH -c 8

module load python2
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS

#define list of samples
python_mod=/cluster/tools/software/centos7/python/2.7.15/bin/python2
writeresults=/cluster/home/kisaev/phylowgs/write_results.py

sample=$1
echo $sample
#rm /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/phyloWGS_$sample/phylowgs_output_$sample/chain_0/trees.zip /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/phyloWGS_$sample/phylowgs_output_$sample/chain_1/trees.zip
	
mkdir phyloWGS_output_$sample
#To write JSON results, please run: 
cd phyloWGS_output_$sample
$python_mod $writeresults $sample /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/phyloWGS_$sample/phylowgs_output_$sample/trees.zip $sample.summ.json.gz $sample.muts.json.gz $sample.mutass.zip

