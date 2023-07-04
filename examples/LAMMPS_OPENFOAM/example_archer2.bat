#!/bin/bash
#SBATCH --job-name=my_cpl_demo
#SBATCH --time=0:30:0
#SBATCH --exclusive
#SBATCH --export=none
#SBATCH --account=y23
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --nodes=2

# single thread export overriders any declaration in srun
export OMP_NUM_THREADS=1

module load other-software
module load cpl-openfoam
source $FOAM_CPL_APP/SOURCEME.sh
module load cpl-lammps

# using your own installation: remove the previous four lines and use these five 'module' lines and three 'source' lines instead
# remmeber to update the path to the three SOURCEME.sh files
#module load openfoam/com/v2106
#module switch gcc gcc/10.3.0
#module load cray-fftw
#module load cray-python
#module load matplotlib
#source /work/y23/y23/gavincpl/cpl-library/SOURCEME.sh
#source /work/y23/y23/gavincpl/CPL_APP_OPENFOAM/SOURCEME.sh
#source /work/y23/y23/gavincpl/CPL_APP_LAMMPS-DEV/SOURCEME.sh

cd openfoam
python clean.py -f
blockMesh
decomposePar
cd ..

SHARED_ARGS="--distribution=block:block --hint=nomultithread"

srun ${SHARED_ARGS} --het-group=0 --nodes=1 --tasks-per-node=2  CPLIcoFoam -case ./openfoam -parallel : --het-group=1 --nodes=1 --tasks-per-node=2 lmp_cpl -i lammps.in
