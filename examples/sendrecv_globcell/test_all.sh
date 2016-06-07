#!/bin/bash

#Setup variable
PWD=$(pwd)
echo $PWD
CPL_PATH=${PWD}/../../
export CPL_PATH

#Build cpp codes
cd ./cpp
rm -rf md cfd
./build.sh
cd ./../

#Build fortran code
cd ./fortran
rm -rf md cfd
./build.sh
cd ./../

#Setup python paths
source ${CPL_PATH}/SOURCEME.sh

#Try all permutations of codes
mpiexec -n 16 ./fortran/md : -n 4 ./fortran/cfd
mpiexec -n 16 ./fortran/md : -n 4 ./cpp/cfd
mpiexec -n 16 ./fortran/md : -n 4 python ./python/cfd_send_cells.py
mpiexec -n 16 ./cpp/md : -n 4 ./fortran/cfd
mpiexec -n 16 ./cpp/md : -n 4 ./cpp/cfd
mpiexec -n 16 ./cpp/md : -n 4 python ./python/cfd_send_cells.py
mpiexec -n 16 python ./python/md_recv_cells.py : -n 4 ./fortran/cfd
mpiexec -n 16 python ./python/md_recv_cells.py : -n 4 ./cpp/cfd
mpiexec -n 16 python ./python/md_recv_cells.py : -n 4 python ./python/cfd_send_cells.py
