#!/bin/bash

#Setup variable
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CPL_PATH=${DIR}/../../
export CPL_PATH

#Build fortran code
cd ${DIR}
rm -rf md cfd
./build.sh

#Try all permutations of codes
cd ${DIR}
rm -f vg_*
mpiexec -n 4 valgrind --leak-check=full --log-file='vg_md.%q{PMI_RANK}' ./md : -n 2 valgrind --leak-check=full --log-file='vg_cfd.%q{PMI_RANK}' ./cfd
cat vg_* | grep 'definitely'
