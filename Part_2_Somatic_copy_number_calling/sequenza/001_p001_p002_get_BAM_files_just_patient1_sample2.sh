#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J merging_vars

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

MYVAR=LY_RAP_0001_Aut_FzT_02_files/gatk/LY_RAP_0001_Aut_FzT_02.sorted.dup.recal.cram
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

#sample CRAM file
bam_file=LY_RAP_0001_Aut_FzT_02_files/gatk/LY_RAP_0001_Aut_FzT_02.sorted.dup.recal.cram
samtools view -T $fasta -bh $bam_file > ${tum_name}.bam
samtools index ${tum_name}.bam

#then BAM files were moved into /cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files
