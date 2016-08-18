#!/bin/bash

#Setup variable
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CPL_PATH=${DIR}/../../
export CPL_PATH

#Build cpp codes
cd ${DIR}/cpp
rm -rf md cfd
./build.sh

#Build fortran code
cd ${DIR}/fortran
rm -rf md cfd
./build.sh

#Setup python paths
source ${CPL_PATH}/SOURCEME.sh

cd ${DIR}
#mpiexec -n 4 ./fortran/md : -n 1 python ./python/CFD_recv_and_plot.py
#mpiexec -n 4 ./cpp/md : -n 1 python ./python/CFD_recv_and_plot.py
#mpiexec -n 4 ./fortran/md : -n 1 python ./python/CFD_recv_and_plot_grid.py
mpiexec -n 4 ./cpp/md : -n 1 python ./python/CFD_recv_and_plot_grid_interactive.py
