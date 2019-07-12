#!/bin/bash
#SBATCH -t 100:0:0
#SBATCH -J proc_vcfs

#setting_up_individuakl_pyclone_inputs.sh
#author: Karin Isaev
#date started: May 28, 2019

module load R 

#loop over 8/10 annotated vcf files and turn them into simplified bed files for intersection with cnvkit output

for i in `seq 1 20`; do
    echo $i
    sbatch -n 1 -t 2-00:00 --mem=61440M -p himem -J proc_vcfs$i --wrap="Rscript /cluster/home/kisaev/scripts/processing_annovar_results.R $i"
    sleep 1 # pause to be kind to the scheduler
done

