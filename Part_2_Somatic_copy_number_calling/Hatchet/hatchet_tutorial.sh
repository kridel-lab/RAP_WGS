#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=41440M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J hatchet_tut


#
#Preliminaries
#This preliminary part of the script contains all the preliminary
#information that are required to execute the full pipeline.
#

export GUROBI_HOME="/cluster/home/kisaev/gurobi911/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
export GRB_LICENSE_FILE="/cluster/home/kisaev/gurobi.lic"

#one time during gurobi installation
#cd "${GUROBI_HOME}"
#cd src/build/
#make
#cp libgurobi_c++.a ../../lib

conda activate hatchet

REF="/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta"
SAM="/cluster/tools/software/centos7/samtools/1.10/bin/"
BCF="/cluster/tools/software/centos7/samtools/1.10/bin/"

XDIR="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Hatchet/p001/"
NORMAL="/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam"
BAMS="/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Aut_FzT_01_files/gatk/LY_RAP_0003_Aut_FzT_01.sorted.dup.recal.cram.bam /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Aut_FzT_02_files/gatk/LY_RAP_0003_Aut_FzT_02.sorted.dup.recal.cram.bam"
ALLNAMES="Normal Aut_FzT_01 Aut_FzT_02"
NAMES="Aut_FzT_01 Aut_FzT_02"
J=22

set -e
set -o xtrace
PS4='\''[\t]'\'
export PATH=$PATH:${SAM}
export PATH=$PATH:${BCF}
#source /path/to/virtualenv-python2.7/bin/activate

#
#Setting up running directory
#

BIN=${XDIR}bin/
mkdir -p ${BIN}
BAF=${XDIR}baf/
mkdir -p ${BAF}
BB=${XDIR}bb/
mkdir -p ${BB}
BBC=${XDIR}bbc/
mkdir -p ${BBC}
ANA=${XDIR}analysis/
mkdir -p ${ANA}
RES=${XDIR}results/
mkdir -p ${RES}
EVA=${XDIR}evaluation/
mkdir -p ${EVA}

cd ${XDIR}

#
#binBAM
#

\time -v python2 -m hatchet binBAM -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} \
                                   -b 50kb -g ${REF} -j ${J} \
                                   -q 20 -O ${BIN}normal.bin -o ${BIN}bulk.bin -v &> ${BIN}bins.log
