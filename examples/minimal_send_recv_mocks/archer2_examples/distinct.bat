#!/bin/bash

# this batch script runs on Archer2 and assumes two distinct MPI_Comm_world communicators
# at the time of writing, this only runs with two 1-core jobs: Cray are fixing this bug

#SBATCH --job-name=my_cpl_demo
#SBATCH --time=0:02:00
#SBATCH --export=none
#SBATCH --account=ecseaf01
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

#SBATCH hetjob

#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

# single thread export overriders any declaration in srun
export OMP_NUM_THREADS=1

module load openfoam/com/v2106
module load lammps/13_Jun_2022
module load cray-python
module load xthi
cd /work/ecseaf01/ecseaf01/gavboi/cpl-library
source SOURCEME.sh
cd examples/minimal_send_recv_mocks

SHARED_ARGS="--distribution=block:block --hint=nomultithread"

srun ${SHARED_ARGS} --het-group=0 xthi &
srun ${SHARED_ARGS} --het-group=1 xthi &
wait

srun ${SHARED_ARGS} --het-group=0 f_CFD &
time srun ${SHARED_ARGS} --het-group=1 f_MD &
wait


