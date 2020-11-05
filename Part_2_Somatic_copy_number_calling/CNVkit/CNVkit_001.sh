#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J cnvkit
#SBATCH -c 8
#SBATCH --array=0-26 # job array index because 27 total tumour samples

#needs to run on bam files

module load strelka/2.9.10
module load python
module load manta/1.6.0
module load CNVkit
module load samtools

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}
control_file=/cluster/projects/burst2/CNVkit_WORKDIR/Normal_Samples/${patient_name}.bam

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

mkdir /cluster/projects/burst2/CNVkit_WORKDIR/${tum_name}
output=/cluster/projects/burst2/CNVkit_WORKDIR/${tum_name}

samtools view -b  -T $fasta -o ${output}/${tum_name}.bam ${names[${SLURM_ARRAY_TASK_ID}]}

cd $output
cnvkit.py batch ${output}/${tum_name}.bam -n ${control_file} \
        -m wgs -f $fasta --target-avg-size 1000 -p 10
