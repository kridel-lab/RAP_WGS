#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_vars
#SBATCH --array=0-26 # job array index

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
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files

#pwd
names=($(cat sequenza_input_bam_files_tum.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.bam*}
patient_name=${MYVAR%_*_*_*}

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

#file with mutations that need to be identified from each sample (bed_file)
muts_all=/cluster/projects/kridelgroup/RAP_ANALYSIS/data/all_mutations_${patient_name}_pyclone_bam_readcount_input.bed

#sample BAM file
bam_file=${names[${SLURM_ARRAY_TASK_ID}]}

#output file to save
out_put_file_all=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone/${tum_name}_missing_muts_all.bed

#get counts full list of variants
bam-readcount -f $fasta \
 ${tum_name}.bam \
-l $muts_all > $out_put_file_all
