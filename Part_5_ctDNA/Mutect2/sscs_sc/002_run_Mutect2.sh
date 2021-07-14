#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=21440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J MUTECT2
#SBATCH --array=0-2 # job array index

module load java/8  #8
module load samtools
module load python3
module load gatk/4.0.5.1
module load annovar

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing
#ls *sscs.sc.bam > all_sscs_sc_bam_files_RAP

samples=all_sscs_sc_bam_files_RAP
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

name="${sample%%.bam*}"
echo $name

fasta_file=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta
out_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls

targets_interval_list=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/picard_tools_targets_input.bed
ints=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing/${name}_targets.interval_list

tum=($(samtools view -H ${names[${SLURM_ARRAY_TASK_ID}]} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq))
echo "${tum}"
export tum

cont_sample="${tum##*_}"
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files
normal_file=$(ls LY_RAP_${cont_sample}_Ctl*.bam)
echo $normal_file

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing

gatk Mutect2 \
-R $fasta_file \
-I ${names[${SLURM_ARRAY_TASK_ID}]} \
-tumor ${tum} \
-L $ints \
-I /cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/${normal_file} \
-O $out_folder/Mutect2_VCF_output/${tum}.sscs.sc.vcf.gz \
--germline-resource /cluster/projects/kridelgroup/RAP_ANALYSIS/af-only-gnomad.raw.sites.b37.vcf.gz
