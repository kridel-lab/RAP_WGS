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
module load R/4.0.0

cd /cluster/projects/kridelgroup/RAP_ANALYSIS

export PURECN="/cluster/home/kisaev/R/x86_64-pc-linux-gnu-library/4.0/PureCN/extdata"
#Rscript $PURECN/PureCN.R --help

names=($(cat all_bam_files_raw.txt))
echo ${names[${SLURM_ARRAY_TASK_ID}]}

MYVAR=${names[${SLURM_ARRAY_TASK_ID}]}
tum_loc=${MYVAR%/*}
MYVAR=${MYVAR##*/}
tum_name=${MYVAR%.sorted.dup.recal.cram*}
patient_name=${MYVAR%_*_*_*}
control_file=/cluster/projects/burst2/CNVkit_WORKDIR/Normal_Samples/${patient_name}.bam

fasta=/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta #from gatk resource bundle

mkdir /cluster/projects/burst2/PureCN_WORKDIR/${tum_name}
OUT=/cluster/projects/burst2/PureCN_WORKDIR/${tum_name}
IN=/cluster/projects/burst2/CNVkit_WORKDIR/${tum_name}

#VCF file input from Mutect2
vcf_file=/cluster/projects/burst2/MUTECT2_selected_VCFs/${tum_name}.selected.normalized.vcf.gz
anno_input=merged_MUTECT2_STRELKA/merged_variants_vcfs/${tum_name}_merged_variants.vcf
SAMPLEID=${tum_name}

# Export the segmentation in DNAcopy format
cnvkit.py export seg $IN/${SAMPLEID}.cns --enumerate-chroms \
      -o $OUT/${SAMPLEID}.seg

# Run PureCN by providing the *.cnr and *.seg files
Rscript $PURECN/PureCN.R --out $OUT  \
           --sampleid $SAMPLEID \
           --tumor $IN/${SAMPLEID}.cnr \
           --segfile $OUT/${SAMPLEID}.seg \
           --vcf ${vcf_file} \
           --genome hg19 \
           --force --postoptimize --seed 123
