MYPWD=${PWD}
rm -f ./lmp_cpl
cd /home/es205/codes/cpl_granlammmps/OpenFOAM-3.0.1_LAMMPS-dev/LAMMPS-dev_coupled/GranLAMMPS/src
#cd ../../../../GranLAMMPS/src/
make yes-user-cpl
make package-update
make cpl
cp ./lmp_cpl ${MYPWD}/
cd ${MYPWD}
