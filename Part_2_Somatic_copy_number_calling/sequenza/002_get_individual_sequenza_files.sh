#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J sequenza
#SBATCH -c 8
#SBATCH --array=0-26 # job array index because 27 total tumour samples

module load python2
module load samtools
module load tabix

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files

names=($(cat sequenza_input_bam_files_tum.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
#tum_loc=${MYVAR%/*}
tum_name=${MYVAR%.bam*}
patient_name=${MYVAR%_*_*_*}
control_file=$(ls ${patient_name}_Ctl*.bam)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#run sequenza
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Process a FASTA file to produce a GC Wiggle track file:
fasta="/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta"

#From BAM files
#Normal and tumor BAM files

normal=$control_file
tumor=$MYVAR
sample=$tum_name
out_folder="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Sequenza"

#2.
sequenza-utils bam2seqz --normal $normal --tumor $tumor \
    --fasta $fasta -gc hg19.gc50Base.wig.gz --output ${out_folder}/${sample}.seqz.gz

#3.
sequenza-utils seqz_binning --seqz ${out_folder}/${sample}.seqz.gz -o ${out_folder}/${sample}.small.seqz.gz
