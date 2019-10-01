#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J phylowgs
	
module load python2

#test_ssm.txt just has 1000 of the total unique mutations 
#phylowgs run time is proportional to number of mutations 
#less 2019-07-16_ssm_data_RAP_WGS_input.txt | head -2500 > test_ssm.txt

python2 /cluster/home/kisaev/phylowgs/multievolve.py --num-chains 4 --ssms test_ssm.txt --cnvs cnv_data.txt --burnin-samples 1 --mcmc-samples 1

sed -n "$({ echo 1; seq $(wc -l <ssm_data.txt) | sed 1d | shuf | head -n 500 | sort -n; } | sed 's/$/p/')" ssm_data.txt > test_data.txt

python2 /cluster/home/kisaev/phylowgs/multievolve.py --num-chains 4 --ssms test_data.txt --cnvs cnv_data.txt --burnin-samples 1 --mcmc-samples 1
