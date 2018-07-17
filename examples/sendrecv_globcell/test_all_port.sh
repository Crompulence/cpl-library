#!/bin/bash
MYPWD=${PWD}

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

#Try all permutations of codes
cd ${DIR}
mpiexec -n 16 ./fortran/md PID=$! & mpiexec -n 4 ./fortran/cfd ; wait $PID
mpiexec -n 16 ./fortran/md PID=$! & mpiexec -n 4 ./cpp/cfd ; wait $PID
mpiexec -n 16 ./fortran/md PID=$! & mpiexec -n 4 python ./python/cfd_send_cells.py ; wait $PID
mpiexec -n 16 ./cpp/md PID=$! & mpiexec -n 4 ./fortran/cfd ; wait $PID
mpiexec -n 16 ./cpp/md PID=$! & mpiexec -n 4 ./cpp/cfd ; wait $PID
mpiexec -n 16 ./cpp/md PID=$! & mpiexec -n 4 python ./python/cfd_send_cells.py ; wait $PID
mpiexec -n 16 python ./python/md_recv_cells.py PID=$! & mpiexec -n 4 ./fortran/cfd ; wait $PID
mpiexec -n 16 python ./python/md_recv_cells.py PID=$! & mpiexec -n 4 ./cpp/cfd ; wait $PID
mpiexec -n 16 python ./python/md_recv_cells.py PID=$! & mpiexec -n 4 python ./python/cfd_send_cells.py ; wait $PID
cd $MYPWD
