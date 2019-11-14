#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p veryhimem
#SBATCH --mem=175440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J treeomics
#SBATCH -c 12

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/treeomics/src

module load python3/3.6.5
module load texlive/2019

export PYTHONPATH=/cluster/tools/software/centos7/cplex/12.8.0/cplex/python/3.6/x86-64_linux/
export TEXTLIVEPATH=/cluster/tools/software/centos7/texlive/2019/

python3 -c 'import cplex'

normal=LY_RAP_0003_Ctl_FzG_01.hc.vqsr.vcf.gz #4,763,150 predicted germline variants 

#RUN
python treeomics -d input/STRELKA_VCF_INPUTS -n $normal -l 5 -f 0.1 --wes_filtering -o /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Treeomics/Treeomics_WES
