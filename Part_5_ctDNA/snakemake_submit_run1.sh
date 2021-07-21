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
#manually renamed files using this conversion provided by Michael Hong
#KLCS_0084=LY_RAP_0001, KLCS_0085=LY_RAP_0002, KLCS_0086=LY_RAP_0003

#input_folder=/cluster/projects/kridelgroup/210617_A00469_0183_AHCNJCDRXY.KLCS.fastqs
#output_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload

#ln -s ${input_folder}/KLCS_0084_Pl_n_PE_415_TS_210617_A00469_0183_AHCNJCDRXY_1_CTGATCGT-GCGCATAT_R1.fastq.gz ${output_folder}/LY_0001_R1.fastq.gz
#ln -s ${input_folder}/KLCS_0084_Pl_n_PE_415_TS_210617_A00469_0183_AHCNJCDRXY_1_CTGATCGT-GCGCATAT_R2.fastq.gz ${output_folder}/LY_0001_R2.fastq.gz

#ln -s ${input_folder}/KLCS_0085_Pl_n_PE_350_TS_210617_A00469_0183_AHCNJCDRXY_1_ACTCTCGA-CTGTACCA_R1.fastq.gz ${output_folder}/LY_0002_R1.fastq.gz
#ln -s ${input_folder}/KLCS_0085_Pl_n_PE_350_TS_210617_A00469_0183_AHCNJCDRXY_1_ACTCTCGA-CTGTACCA_R2.fastq.gz ${output_folder}/LY_0002_R2.fastq.gz

#ln -s ${input_folder}/KLCS_0086_Pl_n_PE_381_TS_210617_A00469_0183_AHCNJCDRXY_1_TGAGCTAG-GAACGGTT_R1.fastq.gz ${output_folder}/LY_0003_R1.fastq.gz
#ln -s ${input_folder}/KLCS_0086_Pl_n_PE_381_TS_210617_A00469_0183_AHCNJCDRXY_1_TGAGCTAG-GAACGGTT_R2.fastq.gz ${output_folder}/LY_0003_R2.fastq.gz

snakemake -s /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Consensus_Cruncher_Pipeline_Run1.snakefile -j56 --latency-wait 432000 --cluster-config /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/config/cluster.yaml --cluster-sync "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -o {cluster.out}" --rerun-incomplete --nolock
