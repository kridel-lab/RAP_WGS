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
module load tabix
module load vt

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls/Mutect2_VCF_output
#ls *filtered.vcf.gz > all_vcf_files #make once

samples=all_vcf_files
names=($(cat $samples))
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

fasta_file=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta

#normalize file
input_vcf=$sample

#normalize variants, send to standard out and remove duplicates.
vt normalize $sample -r $fasta_file -o ${sample}.normalized.vcf.gz

#RUN ANNOVAR
anno_input=${sample}.normalized.vcf.gz
out_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls/Mutect2_VCF_output/annovar

#name
name="${sample%%_filtered*}"
echo $name

table_annovar.pl --buildver hg19 ${anno_input} /cluster/tools/software/annovar/humandb \
--protocol ensGene,gnomad211_genome,cosmic68,avsnp142 --operation g,f,f,f \
--outfile ${out_folder}/${name}_annovar.vcf.gz --vcfinput
