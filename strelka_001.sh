#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=61440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J strelka
#SBATCH -c 8
#SBATCH --array=0-19 # job array index
#SBATCH -j 20
	
module load strelka/2.9.10  
module load python
module load manta/1.6.0 

STRELKA_INSTALL_PATH=/cluster/tools/software/centos7/strelka/2.9.10
STRELKA_ANALYSIS_PATH=/cluster/projects/kridelgroup/RAP_ANALYSIS/STRELKA_WORKDIR

#This script configures the strelka workflow for somatic variant calling
#    on matched tumor-normal BAM files. The configuration process will
#    produce an analysis makefile and directory structure. The makefile can
#    be used to run the analysis on a workstation or compute cluster via
#    make/qmake or other makefile compatible process.

# configuration
${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam \
    --tumorBam LY_RAP_0003_Aut_FzT_01_files/gatk/LY_RAP_0003_Aut_FzT_01.sorted.dup.recal.cram \
    --referenceFasta /cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta \
    --callRegions \
    --runDir ${STRELKA_ANALYSIS_PATH}
# execution on a single local machine with 20 parallel jobs
${STRELKA_ANALYSIS_PATH}/runWorkflow.py 

