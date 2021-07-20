#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p veryhimem
#SBATCH --mem=175440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J treeomics
#SBATCH -c 2
#SBATCH --array=0-2 # job array index

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/treeomics/src

module load python3/3.6.5
module load texlive/2019

export PYTHONPATH=/cluster/tools/software/centos7/cplex/12.8.0/cplex/python/3.6/x86-64_linux/
export TEXTLIVEPATH=/cluster/tools/software/centos7/texlive/2019/

python3 -c 'import cplex'

names=($(cat /cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_treomics_pats.txt))
pat=${names[${SLURM_ARRAY_TASK_ID}]}
echo $pat

driver_genes=/cluster/projects/kridelgroup/RAP_ANALYSIS/data/${pat}_drivers.csv

str="p001"

if [ "$pat" == "$str" ]; then
  #RUN
  python treeomics -d input/mutect2_strelka_all_muts/${pat} \
  -l 5 --driver_genes=$driver_genes \
  -o /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Treeomics/Treeomics_WGS_5_mp/${pat}
else
  #RUN
  python treeomics -d input/mutect2_strelka_pcgs_only/${pat} \
  -l 5 --driver_genes=$driver_genes --wes_filtering \
  -o /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Treeomics/Treeomics_WGS_5_mp/${pat}
fi
