#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J treeomics
#SBATCH -c 8

module load python3/3.6.5
export PYTHONPATH=/cluster/tools/software/centos7/cplex/12.8.0/cplex/python/3.6/x86-64_linux/
python3 -c 'import cplex'

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/treeomics/src

#RUN

python treeomics -d input/VCFS_treeomics --wes_filtering 


