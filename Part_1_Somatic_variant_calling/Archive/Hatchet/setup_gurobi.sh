#on the cluster

#export GUROBI_HOME="/cluster/home/kisaev/gurobi911/linux64"
#export PATH="${PATH}:${GUROBI_HOME}/bin"
#export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
#export GRB_LICENSE_FILE="/cluster/home/kisaev/gurobi.lic"

#one time during gurobi installation on the cluster
#cd "${GUROBI_HOME}"
#cd src/build/
#make
#cp libgurobi_c++.a ../../lib

#locally on mac

export GUROBI_HOME="/Library/gurobi911/mac64"
export GRB_LICENSE_FILE="/Users/kisaev/gurobi.lic"
export PATH="${PATH}:${GUROBI_HOME}/bin"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

#one time during gurobi installation on the cluster
cd "${GUROBI_HOME}"
cd src/build/
make
cp libgurobi_c++.a ../../lib
