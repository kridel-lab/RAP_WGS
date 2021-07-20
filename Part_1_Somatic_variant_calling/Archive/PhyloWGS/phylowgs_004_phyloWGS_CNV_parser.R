#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J cnv_parse
#SBATCH --array=0-19 # job array index

module load python3

#pwd
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/phyloWGS/CNAs_titan

#find -L . -name "*parser.txt" > cnv_parse_jobs #20

#pwd
names=($(cat cnv_parse_jobs))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

purs=($(cat _phylogwgs_CNV_parser_purities.txt))
echo ${purs[${SLURM_ARRAY_TASK_ID}]}

python /cluster/home/kisaev/phylowgs/parser/parse_cnvs.py ${names[${SLURM_ARRAY_TASK_ID}]} -f titan --cnv-output ${names[${SLURM_ARRAY_TASK_ID}]}.parsed.txt -c ${purs[${SLURM_ARRAY_TASK_ID}]}

