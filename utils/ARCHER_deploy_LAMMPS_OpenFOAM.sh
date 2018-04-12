###################################################
#       COMPILERS AND MODULES
###################################################

#Swap Cray for GCC compilers
module unload PrgEnv-cray
module load PrgEnv-gnu
#Use GNU compiler 5.1 that was used for OpenFOAM - IMPORTANT TO ENSURE CONSISTENT RESULTS
module swap gcc gcc/5.1.0
#Need version of git
module load git/2.16.2_build1
#Need to import python from compute node for subprocess
module load python-compute/2.7.6
export CRAYPE_LINK_TYPE=dynamic

###################################################
#       Get GranLAMMPS 
###################################################

#Get GranLAMMPS first as this requires password
mkdir LAMMPS
cd ./LAMMPS
git clone https://edwardsmith999@bitbucket.org/granlammps/gitlammps.git gitlammps
##########PASSWORD NEEDED HERE####################
#Switch to Granular branch -- common
cd ./gitlammps
git checkout common
cd ../

###################################################
#       Build CPL library
###################################################

#Get CPL library
git clone https://github.com/Crompulence/cpl-library.git 
cd cpl-library
make PLATFORM=ARCHER
#Next we need to build a version of mpi4py from this century
cd ./3rd-party
python ./build_mpi4py_archer.py
cd ../
#Source all files
source SOURCEME.sh
cd ../

###################################################
#       Install/Copy system OpenFOAM
###################################################

#Get OpenFOAM by copying installed version
mkdir OpenFOAM
cd OpenFOAM 
#We need to copy key third party files here, basically scotch for decomposition of parallel domain
rsync -avP /work/y07/y07/cse/OpenFOAM/ThirdParty-3.0.1/scotch_6.0.3 ./ThirdParty-3.0.1
rsync -avP /work/y07/y07/cse/OpenFOAM/ThirdParty-3.0.1/platforms ./ThirdParty-3.0.1

#Next we copy OpenFOAM itself so it can be patched
#rsync -avP /work/y07/y07/cse/OpenFOAM/OpenFOAM-3.0.1 ./
#Try minimal set of required files:
mkdir -p ./OpenFOAM-3.0.1/platforms/linux64GccDPOpt
rsync -avP /work/y07/y07/cse/OpenFOAM/OpenFOAM-3.0.1/platforms/linux64GccDPOpt ./OpenFOAM-3.0.1/platforms
rsync -avP /work/y07/y07/cse/OpenFOAM/OpenFOAM-3.0.1/etc ./OpenFOAM-3.0.1
rsync -avP /work/y07/y07/cse/OpenFOAM/OpenFOAM-3.0.1/wmake ./OpenFOAM-3.0.1
rsync -avP /work/y07/y07/cse/OpenFOAM/OpenFOAM-3.0.1/src ./OpenFOAM-3.0.1

#Download CPL APP for OpenFOAM and apply patch
OpenFOAM_APP_DIR=./CPL_APP_OPENFOAM-3.0.1
git clone https://github.com/Crompulence/CPL_APP_OPENFOAM-3.0.1.git $OpenFOAM_APP_DIR
pwd > $OpenFOAM_APP_DIR/CODE_INST_DIR
sed -i -e 's/export WM_COMPILER=Gcc/export WM_COMPILER=CC/g' ./config/prefs.sh
cd $OpenFOAM_APP_DIR
source SOURCEME.sh
make sedifoam

###################################################
#       Build GranLAMMPS
###################################################

#Download CPL APP for LAMMPS and add package USER-CPL
LAMMPS_APP_DIR=./CPL_APP_LAMMPS-DEV
git clone https://github.com/Crompulence/CPL_APP_LAMMPS-DEV.git $LAMMPS_APP_DIR
cd gitlammps
pwd > ../$LAMMPS_APP_DIR/CODE_INST_DIR
cd ../$LAMMPS_APP_DIR
echo granular >> config/lammps_packages.in
sed -i -e 's/mpicxx/CC/g' ./config/Makefile.cpl
# We need to patch main.cpp to allow MPMD mode as ARCHER mpich does not
# support using MPI_ports
make patch-lammps-Oct17
make -j 16
cd ../

###################################################
#      Make Run directory
###################################################

mkdir run
cp ./OpenFOAM/CPL_APP_OPENFOAM-3.0.1/bin/CPLSediFOAM ./
cp ./LAMMPS/CPL_APP_LAMMPS-DEV/bin/lmp_cpl ./



