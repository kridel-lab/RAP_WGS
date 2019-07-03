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

# get file name
#find -L . -name "*.cram" > jobs

names=($(cat jobs))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
#echo "${tum}"
export tum

#for i in {1..22} {X,Y,MT}; do echo "* "$i; mkdir -p chr/chr${i}; samtools view -bS <( samtools view -h ${names[${SLURM_ARRAY_TASK_ID}]} $i ) > chr/chr${i}/${tum}.out.${i}.bam; done
filename=`echo ${tum}` 

for chrom in `seq 1 22` X Y

do 
	samtools view -T human_g1k_v37_decoy.fasta -bh ${names[${SLURM_ARRAY_TASK_ID}]} ${chrom} > chr/${filename}_${chrom}.cram 
	samtools index chr/${filename}_${chrom}.cram 

done

