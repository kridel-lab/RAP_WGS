#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -J combinefastq

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher

#old files (May 2021)
ls -lt /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/test_files

#newer files (July 2021)
ls -lt /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload

#mkdir RAP_ctDNA_combined_two_runs/fastq_files
output_dir=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/RAP_ctDNA_combined_two_runs/fastq_files

cat test_files/LY_0001_R1.fastq.gz June30OICRupload/LY_0001_R1.fastq.gz > $output_dir/LY_0001_R1.fastq.gz
cat test_files/LY_0001_R2.fastq.gz June30OICRupload/LY_0001_R2.fastq.gz > $output_dir/LY_0001_R2.fastq.gz

cat test_files/LY_0002_R1.fastq.gz June30OICRupload/LY_0002_R1.fastq.gz > $output_dir/LY_0002_R1.fastq.gz
cat test_files/LY_0002_R2.fastq.gz June30OICRupload/LY_0002_R2.fastq.gz > $output_dir/LY_0002_R2.fastq.gz

cat test_files/LY_0003_R1.fastq.gz June30OICRupload/LY_0003_R1.fastq.gz > $output_dir/LY_0003_R1.fastq.gz
cat test_files/LY_0003_R2.fastq.gz June30OICRupload/LY_0003_R2.fastq.gz > $output_dir/LY_0003_R2.fastq.gz
