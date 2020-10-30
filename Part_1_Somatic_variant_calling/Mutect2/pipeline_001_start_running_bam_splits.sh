#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=51440M
#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-19 # job array index

#need to run Mutect2 on all bam files from tumour samples
#author: Karin Isaev
#date started: June 25, 2019

module load java/8  #8
module load samtools
module load python
module load gatk
module load annovar

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

names=($(cat all_control_samples.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}
MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal*}
patient_name=${MYVAR%_*_*_*}
output=/cluster/projects/kridelgroup/RAP_ANALYSIS/MUTECT2_WORKDIR/CHR_split/

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

#for i in {1..22} {X,Y,MT}; do echo "* "$i; mkdir -p chr/chr${i}; samtools view -bS <( samtools view -h ${names[${SLURM_ARRAY_TASK_ID}]} $i ) > chr/chr${i}/${tum}.out.${i}.bam; done
filename=`echo ${tum}`

for chrom in `seq 1 22` X Y

do
	samtools view -T $fasta -bh ${names[${SLURM_ARRAY_TASK_ID}]} ${chrom} > ${output}/${patient_name}_${chrom}.bam
	samtools index ${output}/${patient_name}_${chrom}.bam

done
