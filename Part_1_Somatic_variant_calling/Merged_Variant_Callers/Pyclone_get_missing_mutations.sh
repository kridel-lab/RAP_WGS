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
cd /cluster/projects/kridelgroup/RAP_ANALYSIS

#pwd
names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}
control_file=$(ls -d ${patient_name}_Ctl*)
str="LY_RAP_0003"

if [ "$patient_name" == "$str" ]; then
  control_file=$(ls $control_file/gatk/*.bam)
else
  control_file=$(ls $control_file/gatk/*.cram)
fi
fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

#file with mutations that need to be identified from each sample (bed_file)
muts_small=/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_small_subset_${patient_name}_pyclone_bam_readcount_input.bed
muts_all=/cluster/projects/kridelgroup/RAP_ANALYSIS/data/pyclone_full_${patient_name}_pyclone_bam_readcount_input.bed

#sample BAM file
bam_file=${names[${SLURM_ARRAY_TASK_ID}]}

#output file to save
out_put_file_small=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone/${tum_name}_missing_muts_small.bed
out_put_file_all=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Pyclone/${tum_name}_missing_muts_all.bed

samtools view -T $fasta -bh ${names[${SLURM_ARRAY_TASK_ID}]} > ${tum_name}.bam
samtools index ${tum_name}.bam

#get counts small list of variants
bam-readcount -f $fasta \
${tum_name}.bam \
-l $muts_small > $out_put_file_small

#get counts full list of variants
bam-readcount -f $fasta \
 ${tum_name}.bam \
-l $muts_all > $out_put_file_all

rm ${tum_name}.bam
