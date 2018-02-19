MYPWD=${PWD}
rm -f ./lmp_cpl
ROOTDIR=/home/es205/codes/cpl_granlammmps/OpenFOAM-3.0.1_LAMMPS-dev/
rm -f ./lmp_cpl
cd ${ROOTDIR}/LAMMPS-dev_coupled/CPL_APP_LAMMPS-DEV/
make 
cp ./bin/lmp_cpl ${MYPWD}
cd ${MYPWD}

