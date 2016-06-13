#!/bin/bash

#Setup variable
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CPL_PATH=${DIR}/../../
export CPL_PATH

#Setup python paths
source ${CPL_PATH}/SOURCEME.sh

cd ${DIR}
mpiexec -n 1 python ./cfd_cpl.py : -n 1 python ./md_cpl.py
