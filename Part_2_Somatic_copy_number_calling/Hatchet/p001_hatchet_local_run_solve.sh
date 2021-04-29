#
#Preliminaries
#This preliminary part of the script contains all the preliminary
#information that are required to execute the full pipeline.
#

export GUROBI_HOME="/Library/gurobi911/mac64"
export GRB_LICENSE_FILE="/Users/kisaev/gurobi.lic"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

XDIR="/Users/kisaev/Documents/Hatchet_analysis/p001"

#set -e
#set -o xtrace
#PS4='\''[\t]'\'
conda activate hatchet

#
#Setting up running directory
#

BIN=${XDIR}/bin
#mkdir -p ${BIN}
BAF=${XDIR}/baf
#mkdir -p ${BAF}
BB=${XDIR}/bb
#mkdir -p ${BB}
BBC=${XDIR}/bbc
#mkdir -p ${BBC}
ANA=${XDIR}/analysis
#mkdir -p ${ANA}
RES=${XDIR}/results
#mkdir -p ${RES}
EVA=${XDIR}/evaluation
#mkdir -p ${EVA}

cd ${XDIR}

#RUN HATCHET SOLVE
cd ${RES}
J=22
\time python2 -m hatchet solve -i ${BBC}/bulk -n2,6 -p 100 -v 3 -u 0.05 -r 12 -j ${J} -eD 6 -eT 12 -g 0.35 -l 0.6 &> >(tee >(grep -v Progress > hatchet.log))

## Increase -l to 0.6 to decrease the sensitivity in high-variance or noisy samples, and decrease it to -l 0.3 in low-variance samples to increase the sensitivity and explore multiple solutions with more clones.
## Increase -u if solutions have clone proportions equal to the minimum threshold -u
## Decrease the number of restarts to 200 or 100 for fast runs, as well as user can decrease the number of clones to -n 2,6 when appropriate or when previous runs suggest fewer clones.
## Increase the single-clone confidence to `-c 0.6` to increase the confidence in the presence of a single tumor clone and further increase this value when interested in a single clone.

cd ${EVA}
\time python2 -m hatchet BBeval ${RES}/best.bbc.ucn
