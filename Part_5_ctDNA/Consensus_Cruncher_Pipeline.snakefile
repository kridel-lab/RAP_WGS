import pandas as pd
from pathlib import Path
import subprocess
from os.path import join

configfile: "/cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/config/cluster.yaml"
configfile: "/cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/config/config_all_samples.json"

### Globals ---------------------------------------------------------------------

# A Snakemake regular expression matching fastq files.

SAMPLES, = glob_wildcards(join(config["fastqDir"], "{sample}_R1."+config["fastqExtension"]))
print(SAMPLES)

wildcard_constraints:
    sample = "\w+"

### Rules -----------------------------------------------------------------------

# Pipeline output files
rule all:
    input:
        expand(join(config["consensusDir"], "{sample}/dcs_sc/{sample}.dcs.sc.sorted.bam"), sample=SAMPLES)

#Consensus cruncher to convert FASTQ files to bam files
rule consensus_crunch_fastq2bam:
    input:
        fq1 = join(config["fastqDir"], "{sample}_R1.fastq.gz"),
        fq2 = join(config["fastqDir"], "{sample}_R2.fastq.gz")
    params:
        fq1_ID = "{sample}_R1.fastq.gz",
        fq2_ID = "{sample}_R2.fastq.gz"
    output:
        bamsort = join(config["bamDir"], "{sample}.sort.bam"),
        bamout = join(config["bamDir"], "{sample}.sorted.bam"),
        bamindex = join(config["bamDir"], "{sample}.sorted.bam.bai"),
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

                export TMPDIR=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
                export TEMP=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
                export TMP=/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/

                #run consensus cruncher
                python /cluster/home/kisaev/ConsensusCruncher/ConsensusCruncher.py fastq2bam --fastq1 {input.fq1}  --fastq2 {input.fq2}  -o /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload/ConsensusCruncher -r /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta -b /cluster/tools/software/bwa/0.7.15/bwa  -s /cluster/tools/software/centos7/samtools/1.10/bin/samtools -l /cluster/home/kisaev/RAP_WGS/Part_5_ctDNA/config/IDT_dual_Index.txt -g /cluster/tools/software/picard/2.10.9/picard.jar -tmp /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
            """

#Consensus Cruncher to get summar statistics
rule consensus_crunch_consensus:
    input:
        bam = rules.consensus_crunch_fastq2bam.output.bamout
    output:
        ssSCbam = join(config["consensusDir"], "{sample}/sscs_sc/{sample}.sscs.sc.sorted.bam"),
        dcsSCbam = join(config["consensusDir"], "{sample}/dcs_sc/{sample}.dcs.sc.sorted.bam"),
        dcsSCuniq = join(config["consensusDir"], "{sample}/dcs_sc/{sample}.all.unique.dcs.sorted.bam"),
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

                #run consensus cruncher
                python /cluster/home/kisaev/ConsensusCruncher/ConsensusCruncher.py consensus --bedfile /cluster/home/kisaev/ConsensusCruncher/ConsensusCruncher/hg19_cytoBand.txt -i {input.bam} -o /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/June30OICRupload/ConsensusCruncher/consensus -s /cluster/tools/software/centos7/samtools/1.10/bin/samtools -tmp /cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/ConsensusCruncher/temp_files/
            """
