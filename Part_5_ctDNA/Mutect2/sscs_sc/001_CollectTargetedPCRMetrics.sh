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

cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/consensus_output/consensus

#get list of sample bam files
#ls */sscs_sc/*.sscs.sc.sorted.bam  > sscs_sc_bam_files_for_mutect_all_samples.txt

samples=sscs_sc_bam_files_for_mutect_all_samples.txt
names=($(cat $samples))

#get sample file
sample=${names[${SLURM_ARRAY_TASK_ID}]}
echo $sample

#extract just the bam file name which includes sample name
name="${sample%%/*}"
echo $name

#prepare files and folders required for analysis and storage of downstream files
fasta_file=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta
#where to store output (summary of coverage and pcr metrics from each bam file)
out_folder=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/mutation_calls
#coordinates of probes (coding and non-coding) used as provided by IDT and Robert
amplicon_interval_list=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/picard_tools_amps_input.bed
#coordinates of targets (coding and non-coding) used as provided by IDT and Robert
targets_interval_list=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/probe_coords/picard_tools_targets_input.bed

#first create symbolic link for bam file so that can index it
cd /cluster/projects/kridelgroup/RAP_ANALYSIS
ln -s /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/consensus_output/consensus/$sample /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing/$name.sscs.sc.bam
cd /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing
sample=$name.sscs.sc.bam #linked bam file that we can use

#index BAM file
samtools index $sample

#make sample specific interval list files required for picard as it doesn't work
#with bed files
#SD is dictionary which in this case is just the BAM file

gatk BedToIntervalList \
      -I $amplicon_interval_list \
      -O /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing/${name}.sscs.sc_amplicon.interval_list \
      -SD $sample

gatk BedToIntervalList \
      -I $targets_interval_list \
      -O /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing/${name}.sscs.sc_targets.interval_list \
      -SD $sample

#save new interval lists as variables which will be used as input for final picard function
amps=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing/${name}.sscs.sc_amplicon.interval_list
ints=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/Mutect2/processing/${name}.sscs.sc_targets.interval_list

gatk CollectTargetedPcrMetrics \
       -I $sample \
       -O ${out_folder}/${name}.sscs.sc.output_pcr_metrics.txt \
       -R $fasta_file \
       --PER_TARGET_COVERAGE ${out_folder}/${name}.sscs.sc.per_target_coverage.txt \
       --AMPLICON_INTERVALS $amps \
       --TARGET_INTERVALS $ints
