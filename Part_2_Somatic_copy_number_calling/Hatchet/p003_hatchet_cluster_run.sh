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
export PATH=/cluster/home/kisaev/bin:$PATH #this needs to

#one time during gurobi installation
#cd "${GUROBI_HOME}"
#cd src/build/
#make
#cp libgurobi_c++.a ../../lib

REF="/cluster/projects/kridelgroup/RAP_ANALYSIS/human_g1k_v37_decoy.fasta"
SAM="/cluster/home/kisaev/bin/samtools"
BCF="/cluster/home/kisaev/bin/bcftools"

XDIR="/cluster/projects/kridelgroup/RAP_ANALYSIS/ANALYSIS/Hatchet/p003/"

NORMAL="/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Ctl_FzG_01.bam"

BAMS="/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Dia_FoT_01.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_13.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Dia_FoT_05.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_05.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_12.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_07.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_14.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_18.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_04.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_11.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_16.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_15.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_17.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_06.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_02.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Dia_FoT_03.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_03.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_01.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_09.bam \
/cluster/projects/kridelgroup/RAP_ANALYSIS/CRAM_to_BAM_converted_files/LY_RAP_0003_Aut_FzT_10.bam"

ALLNAMES="Normal Dia_FoT_01 Aut_FzT_13 Dia_FoT_05 Aut_FzT_05 Aut_FzT_12 Aut_FzT_07 Aut_FzT_14 Aut_FzT_18 Aut_FzT_04 Aut_FzT_11 \
Aut_FzT_16 Aut_FzT_15 Aut_FzT_17 Aut_FzT_06 Aut_FzT_02 Dia_FoT_03 Aut_FzT_03 Aut_FzT_01 Aut_FzT_09 Aut_FzT_10"
NAMES="Dia_FoT_01 Aut_FzT_13 Dia_FoT_05 Aut_FzT_05 Aut_FzT_12 Aut_FzT_07 Aut_FzT_14 Aut_FzT_18 Aut_FzT_04 Aut_FzT_11 \
Aut_FzT_16 Aut_FzT_15 Aut_FzT_17 Aut_FzT_06 Aut_FzT_02 Dia_FoT_03 Aut_FzT_03 Aut_FzT_01 Aut_FzT_09 Aut_FzT_10"
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
#                                   -b 150kb -g ${REF} -j ${J} \
#                                   -q 20 -O ${BIN}normal.bin -o ${BIN}bulk.bin -v &> ${BIN}bins.log

#
#deBAF
#

#\time -v python2 -m hatchet deBAF -N ${NORMAL} -T ${BAMS} -S ${ALLNAMES} \
#            --snps /cluster/projects/kridelgroup/RAP_ANALYSIS/deBAF_input_gnomad_snps.txt \
#            -r ${REF} -j ${J} -q 20 -Q 20 -U 20 -c 10 \
#            -C 300 -O ${BAF}normal.baf -o ${BAF}bulk.baf -v \
#                          &> ${BAF}bafs.log

#
#comBBo
#

\time -v python2 -m hatchet comBBo -c ${BIN}normal.bin -C ${BIN}bulk.bin -B ${BAF}bulk.baf -m MIRROR -e 12 > ${BB}bulk.bb

#
#cluBB-cluBB globally clusters genomic bins based on RDR and BAF jointly along the genome and across
#all tumor samples, specified in a BB file ${BB}bulk.bb. The home directory of BNPY is specified through ${BNPY}
#to perform a Dirichelt-process clustering.
#

\time -v python2 -m hatchet cluBB ${BB}bulk.bb -o ${BBC}bulk.seg -O ${BBC}bulk.bbc -e ${RANDOM} -sf 0.1 -tB 0.02 -tR 0.4 -d 0.08 -d 0.07

cd ${ANA}
\time -v python2 -m hatchet BBot -c RD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot -c CRD --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot  -c BAF --figsize 6,3 ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot  -c BB ${BBC}bulk.bbc &
\time -v python2 -m hatchet BBot  -c CBB ${BBC}bulk.bbc -tS 0.03 &
wait

#cd ${RES}

#\time -v python2 -m hatchet solve -i ${BBC}bulk -n1,5 -p 400 -v 3 -u 0.03 -r ${RANDOM} -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.3 &> >(tee >(grep -v Progress > hatchet.log))

## Increase -l to 0.6 to decrease the sensitivity in high-variance or noisy samples, and decrease it to -l 0.3 in low-variance samples to increase the sensitivity and explore multiple solutions with more clones.
## Increase -u if solutions have clone proportions equal to the minimum threshold -u
## Decrease the number of restarts to 200 or 100 for fast runs, as well as user can decrease the number of clones to -n 2,6 when appropriate or when previous runs suggest fewer clones.
## Increase the single-clone confidence to `-c 0.6` to increase the confidence in the presence of a single tumor clone and further increase this value when interested in a single clone.

#cd ${EVA}
#\time -v python2 -m hatchet BBeval ${RES}/best.bbc.ucn
