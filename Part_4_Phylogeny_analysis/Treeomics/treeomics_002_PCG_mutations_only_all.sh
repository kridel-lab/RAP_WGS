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

#normal=LY_RAP_0003_Ctl_FzG_01.hc.vqsr.vcf.gz #4,763,150 predicted germline variants

#purity_info=/cluster/projects/kridelgroup/RAP_ANALYSIS/TITAN_CNA/results/titan/hmm/optimalClusterSolution_files/titanCNA_ploidy2/annotation_data_palimpsest_input.txt

#remove header/column names from purity file
#tail -n +2 $purity_info > purities_no_cols.txt

#sample names
#awk -F"\t" '{print $1}' purities_no_cols.txt > treeomics_samples_include.txt

#purities
#awk -F"\t" '{print $3}' purities_no_cols.txt > treeomics_samples_purities.txt

driver_genes=/cluster/projects/kridelgroup/RAP_ANALYSIS/data/${pat}_drivers.csv

#RUN
python treeomics -d input/mutect2_strelka_all_muts/${pat} \
-l 3 --wes_filtering \
--driver_genes=$driver_genes \
-o /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Treeomics/Treeomics_WGS/${pat}

#--include `cat treeomics_samples_include.txt` \
#--purities `cat treeomics_samples_purities.txt` \
