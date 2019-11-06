#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J phylowgs_filt_reduced
#SBATCH -c 8
	
module load python2

#the way the files are listed below is the order in which the samplse are represented in the SSM input file 

#make empty CNV file
#run phylowgs
phylowgs=/cluster/home/kisaev/phylowgs/multievolve.py 
#less ssm_data.txt | head -100 > test_ssm_data_input.txt
#python2 $phylowgs --num-chains 4 --ssms ssm_data_additional_soft_filters.txt --cnvs test_cnv_data.txt 
#python2 $phylowgs --num-chains 4 --ssms test_data_ssm_cleaned.txt --cnvs cnv_data.txt --burnin-samples 1 --mcmc-samples 1

#python2 $phylowgs --num-chains 4 --ssms ssm_data.txt --cnvs cnv_data.txt

python2 $phylowgs --num-chains 4 --ssms q --cnvs cnv_data.txt -O testing_all_data_runtime_low --burnin-samples 1 --mcmc-samples 1 