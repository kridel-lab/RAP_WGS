#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J titan

module load samtools
module load python3
module load R 

###-------PART 1---------------------------------------------------------------------------------------
cluster=/cluster/home/kisaev/TitanCNA/scripts/snakemake/config/cluster_slurm.yaml
#icho=/cluster/home/kisaev/TitanCNA/scripts/snakemake/ichorCNA.snakefile
#snakemake -s $icho --cluster-config $cluster --cluster "sbatch -p himem --mem=51440M -t 5-00:00" -j 50

#done! after finished running, comment out "icho" and snakemake command above and run the icho and 
#snakemake below then repeat for part 3 


###-------PART 2---------------------------------------------------------------------------------------
icho=/cluster/home/kisaev/TitanCNA/scripts/snakemake/getAlleleCounts.snakefile
snakemake -s $icho --cluster-config $cluster --cluster "sbatch -p himem --mem=51440M -t 5-00:00" -j 50


###-------PART 3---------------------------------------------------------------------------------------
#icho=/cluster/home/kisaev/TitanCNA/scripts/snakemake/TitanCNA.snakefile
#snakemake -s $icho --cluster-config $cluster --cluster "sbatch -p himem --mem=51440M -t 5-00:00" -j 50


