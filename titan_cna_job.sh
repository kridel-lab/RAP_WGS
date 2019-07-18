icho=/cluster/home/kisaev/TitanCNA/scripts/snakemake/ichorCNA.snakefile
cluster=/cluster/home/kisaev/TitanCNA/scripts/snakemake/config/cluster_slurm.yaml

#snakemake -s $icho --cluster-config $cluster --cluster "sbatch -p himem --mem=51440M -t 5-00:00" -j 50
module load samtools
module load python3

icho=/cluster/home/kisaev/TitanCNA/scripts/snakemake/getAlleleCounts.snakefile
snakemake -s $icho --cluster-config $cluster --cluster "sbatch -p himem --mem=51440M -t 5-00:00" -j 50



