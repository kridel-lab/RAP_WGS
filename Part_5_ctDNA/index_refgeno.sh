#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J Index_RefGenome

#May 5
module purge
module load bwa/0.7.15
module load samtools

#symlink to ref geno here
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/

#ln -s /cluster/projects/kridelgroup/genome_files/gatk_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta /cluster/projects/kridelgroup/LIBERATE/gpersad/reference_genomes/hg38/hg38.fa
bwa index human_g1k_v37_decoy.fasta
