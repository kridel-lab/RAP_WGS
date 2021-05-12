#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=20000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J ConsensusCruncher
#SBATCH -o /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/slurm_outputs/%x-%j.out #redirect job output (both stdout and stderr) to a file called “<job name>-<job id>.out”

#load required modules
module purge
module load picard/2.10.9
module load bwa/0.7.15
module load python3/3.7.2
module load samtools
module load snakemake/5.20.1

#set temp directories
export TMPDIR=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
export TEMP=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
export TMP=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/

#one time make symlinks
#ln -s /cluster/projects/kridelgroup/210419_A00469_0169_BH53HYDRXY.KLCS.fastqs/* /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/test_files

#then manually renamed files using this conversion provided by Michael Hong
#KLCS_0084=LY_RAP_0001, KLCS_0085=LY_RAP_0002, KLCS_0086=LY_RAP_0003

snakemake -s /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Consensus_Cruncher_Pipeline.snakefile -j56 --latency-wait 432000 --cluster-config /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/config/cluster.yaml --cluster-sync "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -o {cluster.out}" --rerun-incomplete --nolock
