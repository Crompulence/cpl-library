#!/bin/bash

# this batch script runs on Archer2 and assumes a single shared MPI_Comm_world communicator

#SBATCH --job-name=my_cpl_demo
#SBATCH --time=0:10:0
#SBATCH --exclusive
#SBATCH --export=none
#SBATCH --account=ecseaf01
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=2

# single thread export overriders any declaration in srun
export OMP_NUM_THREADS=1

module load openfoam/com/v2106
module load lammps/13_Jun_2022
module load cray-python
module load xthi

source /work/ecseaf01/ecseaf01/gavboi/cpl-library/SOURCEME.sh

SHARED_ARGS="--distribution=block:block --hint=nomultithread"

srun ${SHARED_ARGS} --het-group=0 --nodes=1 --tasks-per-node=1 xthi : --het-group=1 --nodes=1 --tasks-per-node=1 xthi

srun ${SHARED_ARGS} --het-group=0 --nodes=1 --tasks-per-node=1 f_MD : --het-group=1 --nodes=1 --tasks-per-node=1 f_CFD



