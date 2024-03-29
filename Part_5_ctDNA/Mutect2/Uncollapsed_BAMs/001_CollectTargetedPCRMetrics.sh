#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH -c 6
#SBATCH --mem=40000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J picard
#SBATCH --array=0-2 # job array index - number of jobs = numb of unique samples with top up runs

module load annovar
module load bam-readcount
module load gcc
module load STAR
module load rsem
module load perl
module load python/2.7
module load gatk/4.0.5.1
module load tabix
module load vcftools

#PCR amplicons and target intervals
#prepare using the bed files made by this script:
#https://github.com/kridel-lab/ctdna/blob/main/Analysis/mutect2_prep/001_Prepare_Amplicons_Targets.R

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload/ConsensusCruncher/bamfiles

#get list of sample bam files
#ls *.sorted.bam  > bam_files_for_mutect_all_samples.txt

samples=bam_files_for_mutect_all_samples.txt
names=($(cat $samples))

#get sample file
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

#extract just the bam file name which includes sample name
name="${sample%%.sorted.bam*}"
echo $name

#prepare files and folders required for analysis and storage of downstream files
fasta_file=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta
#where to store output (summary of coverage and pcr metrics from each bam file)
out_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload/ConsensusCruncher/Uncollapse_BAMs_coverage
#coordinates of probes (coding and non-coding) used as provided by IDT and Robert
amplicon_interval_list=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/picard_tools_amps_input.bed
#coordinates of targets (coding and non-coding) used as provided by IDT and Robert
targets_interval_list=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/picard_tools_targets_input.bed

#first create symbolic link for bam file so that can index it
cd /cluster/projects/kridelgroup/RAP_ANALYSIS
ln -s /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload/ConsensusCruncher/bamfiles/$sample ${out_folder}/$name.bam
cd $out_folder
sample=$name.bam #linked bam file that we can use

#index BAM file
samtools index $sample

#make sample specific interval list files required for picard as it doesn't work
#with bed files
#SD is dictionary which in this case is just the BAM file

gatk BedToIntervalList \
      -I $amplicon_interval_list \
      -O ${name}_amplicon.interval_list \
      -SD $sample

gatk BedToIntervalList \
      -I $targets_interval_list \
      -O ${name}_targets.interval_list \
      -SD $sample

#save new interval lists as variables which will be used as input for final picard function
amps=${name}_amplicon.interval_list
ints=${name}_targets.interval_list

gatk CollectTargetedPcrMetrics \
       -I $sample \
       -O ${name}.output_pcr_metrics.txt \
       -R $fasta_file \
       --PER_TARGET_COVERAGE ${name}.per_target_coverage.txt \
       --AMPLICON_INTERVALS $amps \
       --TARGET_INTERVALS $ints
