#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_vars
#SBATCH --array=0-19 # job array index

#mutations were filtered to include only those to be analyzed in pyclone
#some of those mutations are not present across all patients
#read counts need to be obtained at those regions

module load java/8  #8
module load samtools
module load python3
module load gatk
module load annovar
module load bedtools
module load tabix
module load bam-readcount

#pwd
cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#pwd
names=($(cat /cluster/projects/kridelgroup/RAP_ANALYSIS/patient_ids.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
pat=${names[${SLURM_ARRAY_TASK_ID}]}
echo $pat

#file with mutations that need to be identified from each sample (bed_file)
muts=

#sample BAM file
bam_file=

#output file to save 
$out_put_file

#get counts
bam-readcount -f human_g1k_v37_decoy.fasta \
$bam_file \
-l $muts > $out_put_file
