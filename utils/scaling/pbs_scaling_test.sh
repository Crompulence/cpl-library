#!/bin/bash --login
#PBS -N qscript
#PBS -l select=420
#PBS -l walltime=01:00:00
#PBS -A ecse0803

function change_cfd() {
	sed -i 4s/.*/$1/ input
	sed -i 5s/.*/$2/ input
	sed -i 6s/.*/$3/ input
}

function change_md() {
	sed -i 7s/.*/$1/ input
	sed -i 8s/.*/$2/ input
	sed -i 9s/.*/$3/ input
}


CPL_PATH=/work/ecse0803/ecse0803/es205/cpl-library/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CPL_PATH}/lib

module swap PrgEnv-cray PrgEnv-gnu

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

rm -f qscript.*
rm -f output

echo "Starting: `date`"

MAX=5040
for ((i=24;i<=$MAX;i=i*2)); do
    echo aprun -n $i ./cfd : -n $i ./md
done

change_cfd 4 3 2
change_md 4 3 2
aprun -n 24 ./cfd : -n 24 ./md

change_cfd 4 3 4
change_md 4 3 4
aprun -n 48 ./cfd : -n 48 ./md

change_cfd 4 6 4
change_md 4 6 4
aprun -n 96 ./cfd : -n 96 ./md

change_cfd 4 6 4
change_md 4 6 4
aprun -n 192 ./cfd : -n 192 ./md



