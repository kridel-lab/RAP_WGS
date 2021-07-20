#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J sciclone

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/SciClone

for i in {1..20}
do
echo $i
#need to make sure mutation IDs are labelled 0 to n...
#run accessory R script that will do this
sbatch /cluster/home/kisaev/RAP_WGS/Part_4_Phylogeny_analysis/SciClone/sciclone_001_job_submit.sh $i
done
