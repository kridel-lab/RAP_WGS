#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=60000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ConsensusCruncher
#SBATCH -o /cluster/projects/kridelgroup/LIBERATE/gpersad/slurm_outputs/%x-%j.out #redirect job output (both stdout and stderr) to a file called “<job name>-<job id>.out”

#load required modules
module load picard/2.10.9
module load bwa/0.7.15
module load python3/3.7.2
module load samtools
module load snakemake/5.20.1

#set temp directories
export TMPDIR=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
export TEMP=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
export TMP=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/

snakemake -s /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Consensus_Cruncher_Pipeline.snakefile -j56 --latency-wait 86400 --cluster-config /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Config/cluster.yaml --cluster-sync "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -o {cluster.out}" --rerun-incomplete --nolock
