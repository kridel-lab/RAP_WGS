import pandas as pd
from pathlib import Path
import subprocess
from os.path import join

configfile: "/cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Config/cluster.yaml"
configfile: "/cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Config/config_all_samples.json"


### Globals ---------------------------------------------------------------------

# A Snakemake regular expression matching fastq files.

SAMPLES, = glob_wildcards(join(Config["fastqDir"], "{sample}_R1."+Config["fastqExtension"]))
print(SAMPLES)

wildcard_constraints:
    sample = "\w+"

### Rules -----------------------------------------------------------------------

# Pipeline output files
rule all:
    input:
        expand(join(Config["consensusDir"], "{sample}/dcs_SC/{sample}.dcs.sc.sorted.bam"), sample=SAMPLES)

#Consensus cruncher to convert FASTQ files to bam files
rule consensus_crunch_fastq2bam:
    input:
        fq1 = join(Config["fastqDir"], "{sample}_R1.fastq.gz"),
        fq2 = join(Config["fastqDir"], "{sample}_R2.fastq.gz")
    params:
        fq1_ID = "{sample}_R1.fastq.gz",
        fq2_ID = "{sample}_R2.fastq.gz"
    output:
        bamsort = join(Config["bamDir"], "{sample}.sort.bam"),
        bamout = join(Config["bamDir"], "{sample}.sorted.bam"),
        bamindex = join(Config["bamDir"], "{sample}.sorted.bam.bai"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Alignment of {params.fq1_ID} and {params.fq2_ID} with BWA and extraction of UMIs."
    shell:
            """
                #load modules
                module purge
                module load bwa/0.7.15
                module load python3/3.7.2
                module load samtools/1.10
                module load snakemake/5.20.1

                export TMPDIR=/cluster/projects/kridelgroup/LIBERATE/temp_files/
                export TEMP=/cluster/projects/kridelgroup/LIBERATE/temp_files/
                export TMP=/cluster/projects/kridelgroup/LIBERATE/temp_files/

                fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta
                output=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher

                #run consensus cruncher
                python /cluster/home/kisaev/ConsensusCruncher/ConsensusCruncher.py fastq2bam --fastq1 {input.fq1}  --fastq2 {input.fq2}  -o ${output} -r ${fasta} -b /cluster/tools/software/bwa/0.7.15/bwa  -s /cluster/tools/software/centos7/samtools/1.10/bin/samtools -l /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/Config/IDT_dual_Index.txt -g /cluster/tools/software/picard/2.10.9/picard.jar
            """

#Consensus Cruncher to get summar statistics
rule consensus_crunch_consensus:
    input:
        bam = rules.consensus_crunch_fastq2bam.output.bamout
    output:
        sscssingle = join(Config["consensusDir"], "{sample}.sorted/sscs/{sample}.singleton.sorted.bam"),
        sscsbad = join(Config["consensusDir"], "{sample}/sscs/{sample}.badReads.bam"),
        sscsbam = join(Config["consensusDir"], "{sample}/sscs/{sample}.sscs.sorted.bam"),
        ssSCbam = join(Config["consensusDir"], "{sample}/sscs_SC/{sample}.sscs.sc.sorted.bam"),
        ssSCrescue = join(Config["consensusDir"], "{sample}/sscs_SC/{sample}.sscs.rescue.sorted.bam"),
        ssSCsingle = join(Config["consensusDir"], "{sample}/sscs_SC/{sample}.singleton.rescue.sorted.bam"),
        ssSCremain = join(Config["consensusDir"], "{sample}/sscs_SC/{sample}.rescue.remaining.sorted.bam"),
        ssSCuniq = join(Config["consensusDir"], "{sample}/sscs_SC/{sample}.all.unique.sscs.sorted.bam"),
        dcssingle = join(Config["consensusDir"], "{sample}/dcs/{sample}.sscs.singleton.sorted.bam"),
        dcsbam = join(Config["consensusDir"], "{sample}/dcs/{sample}.dcs.sorted.bam"),
        dcsSCbam = join(Config["consensusDir"], "{sample}/dcs_SC/{sample}.dcs.sc.sorted.bam"),
        dcsSCsingle = join(Config["consensusDir"], "{sample}/dcs_SC/{sample}.sscs.sc.singleton.sorted.bam"),
        dcsSCuniq = join(Config["consensusDir"], "{sample}/dcs_SC/{sample}.all.unique.dcs.sorted.bam"),
    wildcard_constraints:
        sample = "\w+"
    message:
        "Fixing bam files"
    shell:
            """
                #load modules
                module purge
                module load bwa/0.7.15
                module load python3/3.7.2
                module load samtools/1.10
                module load snakemake/5.20.1

                output=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/consensus

                #run consensus cruncher
                python /cluster/home/kisaev/ConsensusCruncher/ConsensusCruncher.py consensus -i {input.bam} -o $output -s /cluster/tools/software/centos7/samtools/1.10/bin/samtools -b False
            """
