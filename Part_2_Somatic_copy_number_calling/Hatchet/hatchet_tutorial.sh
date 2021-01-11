#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -p himem
#SBATCH --mem=60000M
#SBATCH -t 5-00:00 # Runtime in D-HH:MM
#SBATCH -J hatchet_tut

#
#Preliminaries
#This preliminary part of the script contains all the preliminary
#information that are required to execute the full pipeline.
#

module load python/2.7
export GUROBI_HOME="/cluster/home/kisaev/gurobi911/linux64"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
export GRB_LICENSE_FILE="/cluster/home/kisaev/gurobi.lic"
export PATH=/cluster/home/kisaev/bin:$PATH

#one time during gurobi installation
#cd "${GUROBI_HOME}"
#cd src/build/
#make
#cp libgurobi_c++.a ../../lib

REF="/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta"
SAM="/cluster/home/kisaev/bin/samtools"
BCF="/cluster/home/kisaev/bin/bcftools"

XDIR="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Hatchet/p001/"
NORMAL="/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Ctl_FzG_01_files/gatk/LY_RAP_0003_Ctl_FzG_01.sorted.dup.recal.bam"
BAMS="/cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Aut_FzT_01_files/gatk/LY_RAP_0003_Aut_FzT_01.sorted.dup.recal.cram.bam /cluster/projects/kridelgroup/RAP_ANALYSIS/LY_RAP_0003_Aut_FzT_02_files/gatk/LY_RAP_0003_Aut_FzT_02.sorted.dup.recal.cram.bam"
ALLNAMES="Normal Aut_FzT_01 Aut_FzT_02"
NAMES="Aut_FzT_01 Aut_FzT_02"
J=22

set -e
set -o xtrace
PS4='\''[\t]'\'
#export PATH=$PATH:${SAM}
#export PATH=$PATH:${BCF}
source /cluster/home/kisaev/python2env/bin/activate

#
#Setting up running directory
#

BIN=${XDIR}bin/
#mkdir -p ${BIN}
BAF=${XDIR}baf/
#mkdir -p ${BAF}
BB=${XDIR}bb/
#mkdir -p ${BB}
BBC=${XDIR}bbc/
#mkdir -p ${BBC}
ANA=${XDIR}analysis/
#mkdir -p ${ANA}
RES=${XDIR}results/
#mkdir -p ${RES}
EVA=${XDIR}evaluation/
#mkdir -p ${EVA}

cd ${XDIR}

#
#binBAM
#

#\time -v python2 -m hatchet binBAM -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} \
#                                   -b 50kb -g ${REF} -j ${J} \
#                                   -q 20 -O ${BIN}normal.bin -o ${BIN}bulk.bin -v &> ${BIN}bins.log

#
#deBAF
#

#\time -v python2 -m hatchet deBAF -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} \
#            -r ${REF} -j ${J} -q 10 -Q 10 -U 10 -c 4 \
#            -C 300 -O ${BAF}normal.baf -o ${BAF}bulk.baf -v \
#                          &> ${BAF}bafs.log

#
#comBBo
#

#\time -v python2 -m hatchet comBBo -c ${BIN}normal.bin -C ${BIN}bulk.bin -B ${BAF}bulk.baf -m MIRROR -e 12 > ${BB}bulk.bb

#
#cluBB-cluBB globally clusters genomic bins based on RDR and BAF jointly along the genome and across
#all tumor samples, specified in a BB file ${BB}bulk.bb. The home directory of BNPY is specified through ${BNPY}
#to perform a Dirichelt-process clustering.
#

#BNPY=/cluster/home/kisaev/bnpy/bnpy/
#\time -v python2 -m hatchet cluBB ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e ${RANDOM} -tB 0.04 -tR 0.15 -d 0.08

python2 -m hatchet cluBB $bbDir"/bulk.bb" -o $bbcDir"/bulk.seg" -O $bbcDir"/bulk.bbc" \
-e 12 -tB 0.02 -tR 0.8 -d 0.08 -R 40 -sf 0.001 -K 50 -v
