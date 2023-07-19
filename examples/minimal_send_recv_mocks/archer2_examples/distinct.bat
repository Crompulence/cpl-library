#!/bin/bash

#SBATCH --job-name=my_cpl_demo
#SBATCH --time=0:02:00
#SBATCH --export=none
#SBATCH --account=y23

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

# following two lines required to enable distinct communicators 
export PMI_UNIVERSE_SIZE=3
export MPICH_SINGLE_HOST_ENABLED=0

module load openfoam/com/v2106
module load gcc/10.3.0
module load cray-python
module load xthi
source /work/y23/y23/gavincpl/cpl-library/SOURCEME.sh
SHARED_ARGS="--distribution=block:block --hint=nomultithread"

srun ${SHARED_ARGS} --het-group=0 xthi &
srun ${SHARED_ARGS} --het-group=1 xthi &
wait
srun ${SHARED_ARGS} --het-group=0 CFD &
srun ${SHARED_ARGS} --het-group=1 MD &
wait



