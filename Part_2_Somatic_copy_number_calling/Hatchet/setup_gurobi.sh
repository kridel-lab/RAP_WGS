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
